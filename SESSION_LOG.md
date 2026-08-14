# Session Log

Committed, dated record of work performed in this repository. Reverse
chronological order, newest entry at the top. See CLAUDE.md for the
convention this file follows.

## 2026-08-14 18:57 UTC — Rao's Q in AnnualSpectralDiversity.R: pyGNDiv/reticulate to rasterdiv

Completed an in-progress edit to `Code/AnnualSpectralDiversity.R` that had
been interrupted mid-way: switched Rao's Q computation from
pyGNDiv/reticulate to the R-native `rasterdiv` package (CRAN).

Rationale: this script compares Rao's Q across YEARS at the same site using
the same variable set (NDVI-only, NIRv-only, or all-bands) each time. It
never needs cross-variable-count comparability (e.g. a 1-band Q vs. an
all-bands Q), which was the reason pyGNDiv's generalizable normalization was
used in the first place (see `ComputeSpecBiodiv.R`). Dropping it also removes
the Python/reticulate dependency, which had been a recurring source of
environment friction (missing `pyGNDiv.py`, no real-data access in this
sandbox to validate against, synthetic-only test harnesses standing in for
real validation).

Work done:
- Verified the actual current `rasterdiv` CRAN API (v0.3.8) rather than
  assuming the 2017 Rocchini et al. paper's older `spectralrao()` signature
  still applies. Confirmed `paRao()` takes a `SpatRaster` directly (no
  conversion needed): `method = "classic"` for single-layer inputs (NDVI,
  NIRv), `method = "multidimension"` for all-bands, which requires a *list*
  of single-layer SpatRasters rather than one multi-band object.
- Found and worked around two footguns in `paRao()` during validation:
  - Its default `simplify = 0` rounds input data to whole numbers before
    computing distances, which silently zeroes out Rao's Q entirely for
    0-1-range reflectance/index data. Fixed with `simplify = 4` (matches
    NEON's native ~4-decimal reflectance scale factor), confirmed to restore
    non-trivial, non-zero Q on a synthetic raster.
  - An internal warning claims NA pixels are "treated as 0s". Verified this
    is stale for the code path actually used (window aggregation sums with
    `na.rm = TRUE`, so NA pixels are excluded, not substituted) using a
    homogeneous synthetic raster with an NA hole in the middle, which
    returned near-zero Q throughout rather than the large spurious Q that
    substituting 0 would produce at the hole's edge.
- Replaced `compute_rao_q_raster_pygndiv()` (and its `get_pca_max_dist()` /
  `rao_q_window_pygndiv()` helpers) with `compute_rao_q_raster_rasterdiv()`,
  using classic (non-normalized) Rao's Q, Euclidean distance, alpha = 1
  (arithmetic mean), and the existing `raoq_window` (3x3) parameter. Removed
  `library(reticulate)`, `pygndiv_path`, and `raoq_stride` (no longer
  meaningful -- `paRao()` computes a native per-pixel moving window rather
  than the coarse non-overlapping tiling pyGNDiv needed for per-window
  reticulate-call runtime reasons). Everything else (CV, CHV, spectral
  species richness, 500m buffer, per-site-year loop, incremental CSV
  writing, output column names) is unchanged.
- Validated with a synthetic test harness (no real `Data/` access in this
  sandbox): extracted and ran the actual `compute_rao_q_raster_rasterdiv()`
  function from the script against synthetic single-band and multi-band
  (with NA holes) rasters. All checks passed -- finite, non-negative Q in
  every case, and a two-region synthetic raster showed clearly elevated Q at
  the region boundary (~0.166) vs. a homogeneous interior (~0.005).

Committed as `969dfc4` (not pushed). Not covered this session: CLAUDE.md
itself is still a mostly-unfilled template (PI, collaborators, repository
URL, SCIENCE_PRINCIPLES.md version/commit, data-source rules, etc. all have
`<FILL>` placeholders) -- worth a dedicated pass before treating its governing
rules as fully adopted for this project.
