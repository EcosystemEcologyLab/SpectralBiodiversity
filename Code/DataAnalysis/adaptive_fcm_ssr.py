"""
adaptive_fcm_ssr.py

Self-adaptive fuzzy c-means (FCM) Spectral Species Richness (SSR), per
Wu et al. (2026, Computers and Electronics in Agriculture 252:112108).
Bridged into R via reticulate::source_python(), the same source_python()
bridging pattern ComputeSpecBiodiv.R already uses for pyGNDiv (that
script's Section 6) -- reused here NOT because SSR itself is already
implemented in Python anywhere in this repo (it isn't:
AnnualSpectralDiversity.R's compute_spectral_species_richness() is pure R,
RF-proximity + k-means + nearest-centroid). This is a narrowly-scoped
bridge, used only by adaptive-FCM SSR scripts, because Wu et al. themselves
used the `fcmeans` (fuzzy-c-means PyPI package) implementation and
reproducing their method means using the same fitting code, not a
different FCM implementation that could converge differently.

Shared by two R callers as of this writing:
  - CompareSSR_AdaptiveFCM_vs_KMeans.R -- real ABBY/SRER comparison against
    the cached/recomputed existing nearest-centroid SSR.
  - FCM_test.R -- companion diagnostic to KMeans_test.R's n_clusters sweep,
    using this module's `hard_assignment` return value (see
    adaptive_fcm_ssr()'s docstring below) for a Shannon's H' comparable in
    FORMULA to the k-means version's, though not in population scope (see
    FCM_test.R's own header for that caveat).
Both call adaptive_fcm_ssr() unchanged; nothing here is specific to either
caller.

Validity function being maximized over candidate c (from the task prompt
this module was written against, matching Wu et al.'s definition):

    L(c) = [ sum_i sum_j u_ij^w ||v_i - xbar||^2 / (c - 1) ]
         / [ sum_i sum_j u_ij^w ||x_j - v_i||^2 / (m - c) ]

  w      = fuzzy weighting exponent (w > 1).
  c      = candidate cluster count.
  m      = number of pixels in the window.
  u_ij   = fuzzy membership of pixel j in cluster i.
  v_i    = cluster i's prototype (fcmeans' `.centers`).
  xbar   = mean of all pixels in the window.

NAME COLLISION WARNING: the paper's `w` (fuzziness exponent) is passed to
the `fcmeans.FCM` constructor as its `m` argument, and the paper's `m`
(pixel count) is a completely different quantity. To avoid ever writing a
bare `m` that could mean either one, this module always spells them out:
`fuzziness` for the paper's w / package's m, and `n_pixels` for the
paper's m (window pixel count).
"""

import numpy as np
from fcmeans import FCM


def _fit_one(X, c, fuzziness, random_state):
    model = FCM(n_clusters=c, m=fuzziness, random_state=random_state, verbose=False)
    model.fit(X)
    return model.u, model.centers


def fcm_validity_L(X, u, centers, fuzziness):
    """Wu et al. validity function L(c) for one already-fit FCM result.

    X: (n_pixels, n_dims) array, the SAME data the FCM was fit on.
    u: (n_pixels, c) fuzzy membership matrix (fcmeans' `.u`).
    centers: (c, n_dims) cluster prototypes (fcmeans' `.centers`).
    fuzziness: the paper's w (== fcmeans' constructor `m`).
    """
    X = np.asarray(X, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    centers = np.asarray(centers, dtype=np.float64)
    n_pixels, c = u.shape

    xbar = X.mean(axis=0)
    u_w = u ** fuzziness  # (n_pixels, c)

    # between-cluster term: sum_i (sum_j u_ij^w) * ||v_i - xbar||^2
    dist_v_xbar_sq = np.sum((centers - xbar) ** 2, axis=1)  # (c,)
    between = float(np.sum(u_w.sum(axis=0) * dist_v_xbar_sq))

    # within-cluster term: sum_i sum_j u_ij^w * ||x_j - v_i||^2
    # ||x - v||^2 expansion (same vectorization technique
    # AnnualSpectralDiversity.R's assign_nearest_centroid() uses in R, done
    # here in numpy instead of a broadcasted (n, c, d) difference array so
    # memory stays O(n*c) rather than O(n*c*d) for wide-band inputs).
    x_sq = np.sum(X ** 2, axis=1)  # (n_pixels,)
    c_sq = np.sum(centers ** 2, axis=1)  # (c,)
    dist_x_v_sq = x_sq[:, None] + c_sq[None, :] - 2.0 * (X @ centers.T)  # (n_pixels, c)
    within = float(np.sum(u_w * dist_x_v_sq))

    numerator = between / (c - 1)
    denominator = within / (n_pixels - c)
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def adaptive_fcm_ssr(X, c_min=2, c_max=20, fuzziness=2.0, random_state=None):
    """
    Fit FCM for every candidate c in [c_min, c_max], compute L(c) for each,
    and take argmax_c L(c) as this window's adaptive-FCM SSR.

    Returns a dict:
      winning_c: int, the argmax_c L(c).
      L_values: dict[int, float], c -> L(c) for every candidate tried.
      n_pixels: int.
      hit_cap: bool, True iff winning_c == c_max -- the caller should log
        this across windows; per the task this script was written for, a
        winning_c that bunches against c_max would mean the cap is too low
        and needs raising, not that c_max itself is the right answer.
      hard_assignment: list[int], length n_pixels, each pixel's HARD
        cluster label (argmax of its fuzzy membership row) from the
        WINNING c's fit only -- i.e. np.argmax(u, axis=1) for that one
        fit, not a re-fit. Provided for callers (e.g. FCM_test.R) that want
        a Shannon's H'-style diversity index comparable to a hard-assignment
        k-means result; this module itself never uses it (the validity
        function L(c) always uses the full fuzzy membership u, never this
        hardened version).
    """
    X = np.asarray(X, dtype=np.float64)
    n_pixels = X.shape[0]
    if n_pixels <= c_max:
        raise ValueError(
            f"adaptive_fcm_ssr: n_pixels ({n_pixels}) must exceed c_max "
            f"({c_max}) -- the validity function's (m - c) denominator "
            f"term requires strictly more pixels than the largest "
            f"candidate c."
        )

    L_values = {}
    u_by_c = {}
    for c in range(c_min, c_max + 1):
        u, centers = _fit_one(X, c, fuzziness, random_state)
        L_values[c] = fcm_validity_L(X, u, centers, fuzziness)
        u_by_c[c] = u

    valid = {c: L for c, L in L_values.items() if np.isfinite(L)}
    if not valid:
        raise RuntimeError(
            "adaptive_fcm_ssr: L(c) was non-finite for every candidate c "
            "in [c_min, c_max] -- investigate the input data/fit, do not "
            "guess a winner."
        )
    winning_c = max(valid, key=valid.get)
    hard_assignment = np.argmax(u_by_c[winning_c], axis=1)

    return {
        "winning_c": int(winning_c),
        "L_values": {int(c): float(L) for c, L in L_values.items()},
        "n_pixels": int(n_pixels),
        "hit_cap": bool(winning_c == c_max),
        "hard_assignment": [int(v) for v in hard_assignment],
    }
