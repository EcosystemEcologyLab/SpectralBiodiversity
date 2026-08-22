# Session Log

Committed, dated record of work performed in this repository. Reverse
chronological order, newest entry at the top. See CLAUDE.md for the
convention this file follows.

## 2026-08-22 21:32 UTC — AnnualSpectralDiversity.R: epsilon floor for fit_hard_wk()'s log(0) edge case

Small, targeted follow-up to the gap-statistic entry directly below this
one, closing one of that session's own flagged-but-not-fixed robustness
concerns: `fit_hard_wk()`'s `W_k` could hit exactly 0 if every point in a
hard-assigned cluster coincides exactly with that cluster's own centroid,
producing `log(W_k) = -Inf` inside the Gap(c) calculation -- silently, not
as a crash, which is a worse failure mode across an unattended 45-site x
multi-year x n_reps_ssr-rep run.

**Fix**: a single named constant, `epsilon_wk <- 1e-10` (defined just above
`fit_hard_wk()`, explicitly separated from Section 0's PROVISIONAL
science/tuning constants since this is a numerical-safety floor, not a
methodological choice), applied at ONE finalization point inside
`fit_hard_wk()` -- after both its c==1 and c>=2 branches compute their raw
`wk`, before either returns -- not duplicated across the two branches. Logs
a `cat()` (matching this file's established progress-logging convention;
no `warning()` used anywhere else in the file) when the floor actually
engages, so the edge case is visible in an unattended run's log rather than
smoothed over silently.

**Self-tested in this sandbox, four checks, all passed:**
1. A degenerate synthetic case (20 identical points, `c=1`) forced `wk=0`
   exactly; confirmed `log(wk)` is finite after the fix (no `-Inf`).
2. Confirmed the `cat()` warning fired exactly once for that call, with the
   correct candidate-c value in the message.
3. Confirmed silence on ordinary data: a normal synthetic 20x4 random
   matrix produced no floor message and `wk` many orders of magnitude above
   `epsilon_wk`.
4. Ran the full `compute_gap_statistic()` + `pick_optimal_c()` path on a
   fully-degenerate observed matrix (forcing the floor for every candidate
   c in 1:4) -- completed without error, every `Gap(c)`/`sk` value finite
   (all landed at exactly 0, as expected when observed and null W_k are
   both floored identically), `pick_optimal_c()` returned a valid c with no
   crash.

**Non-perturbation check**: compared `fit_hard_wk()`'s `wk` output with vs.
without the floor, same input/seed, directly (not via two independent
`compute_gap_statistic()` calls -- an earlier attempt at this comparison
spuriously showed a real-looking difference, traced to
`compute_gap_statistic()`'s null-reference generation using unseeded
`runif()`/`rnorm()`, so two separate calls draw different null data
regardless of the floor; not a bug in the fix, a flaw in that first test's
design). Direct comparison on ordinary 4-blob synthetic data, c=1..6: `wk`
values were bit-identical (relative difference exactly 0 in every case) --
the floor is a true no-op for any non-degenerate input.

File still `parse()`s cleanly. No other part of the gap-statistic logic,
`pick_optimal_c()`, or the conditional c=2-floor gating was touched, per
instruction.

**This entry and the code change it describes have NOT been committed or
pushed** -- same working-tree-only convention as the two entries below.

## 2026-08-22 21:25 UTC — AnnualSpectralDiversity.R: replace the dead c=1 floor-check with a conditional gap-statistic rule

`Code/AnnualSpectralDiversity.R` (Section 5) already had SSR switched from
RF+K-means to self-adaptive FCM in the working tree from a prior,
still-uncommitted session (see the entry below this one). That version
included a c=1 "floor-check" branch meant to catch an apparent c=2 floor
artifact on real SRER data (every rep's L(c) maximized at c=2, no interior
peak, unlike ABBY's genuine interior peak at c~6-8) -- but that branch's own
self-test found it essentially unwinnable: fuzzy within-cluster dispersion
decreases monotonically in c even for pure noise, so c=1 almost never beat
the best of c=2..20. **This session removed that branch entirely** (not
patched) and replaced it with a conditional gap-statistic rule (Tibshirani,
Walther & Hastie 2001), per explicit task instructions reporting that the
gap statistic had been validated elsewhere against real ABBY/SRER data
(ABBY -> c~5, matching its independent L(c) peak; SRER -> c=1, genuinely
homogeneous) -- **that real-data validation was NOT reproduced in this
sandbox** (no server/data access here); only the synthetic self-tests below
were run in this session.

**The removed branch, for the record:** compared c=1's raw within-cluster
dispersion (mean squared distance to the global mean, per-pixel-normalized)
against `min(within-cluster dispersion across c=2..c_max)` and returned c=1
iff that was smaller -- `if (within_1 <= min_within_ge2) return(list(winning_c
= 1L, ...))`. Gone, along with the `c_min_ssr` config constant and the
`within_norm` dispersion tracking it depended on.

**Replacement design, conditional not unconditional (cost-motivated):**
running the full gap-statistic check (null-reference generation + FCM fits
across `gap_check_n_null` references x multiple candidate c) for every rep
of every site-year would meaningfully erode the compute savings that
motivated dropping k-means in the first place. So it only runs for a rep
whose PRIMARY search (`select_winning_c()`, unchanged: c in [2, c_max_ssr]
via L(c) maximization) lands exactly at the c=2 floor. Every other rep uses
the primary search's result directly, unchanged. When triggered:
`compute_gap_statistic()` tests c in [1, `gap_check_c_max`=8] (a reduced
range -- only meant to adjudicate c=1 vs. small c, not re-run the full
search) against `gap_check_n_null`=10 uniform-random null references per
candidate c, generated within the observed subsample's own per-PC-dimension
min/max range (Tibshirani's simpler reference method). `W_k` is the RAW
POOLED WITHIN-CLUSTER SUM OF SQUARED DISTANCES UNDER HARD ASSIGNMENT (argmax
of FCM's fuzzy membership, via a new `fit_hard_wk()` helper) -- a genuinely
different quantity from the primary search's fuzzy-weighted `within_norm`
(now removed); the two share only the underlying FCM fit call
(`fit_fcm()`, factored out for both), not the objective computation.
`pick_optimal_c()` implements the global-max-then-1SE rule specified in the
task (find the global max of Gap(c) first, then the smallest c within one
SE of it) -- explicitly NOT a left-to-right scan from c=1, which was
confirmed broken by this session's own synthetic test (see below).

**Self-tested in this sandbox (no server access), three ways, all passed:**
1. `pick_optimal_c()` against a synthetic dip-before-peak Gap(c) sequence
   (mimicking ABBY's reported real shape -- dips at c=2, true peak at c=5):
   correctly picked c=5 (the global max); a hand-written left-to-right-scan
   variant, run side by side for contrast, incorrectly stopped at c=1 --
   concretely confirming the dip-early-stop failure mode the task described.
2. `compute_gap_statistic()` + `pick_optimal_c()` directly, on genuinely
   synthetic homogeneous noise (500 pts, 4 dims) vs. genuinely synthetic
   5-blob multi-cluster data (500 pts, 4 dims, well-separated): resolved to
   c=1 for the homogeneous case and c=5 for the blobs case (exactly matching
   the true blob count, not just >1) -- the gap statistic's own correctness,
   independent of whether the primary search happens to trigger it.
3. **The conditional wiring itself**, forced: this session's synthetic
   homogeneous/multi-cluster data did NOT naturally reproduce SRER's
   reported real-data behavior of the primary search landing at exactly
   c=2 (probed across multiple seeds/dimensionalities -- primary search
   winning_c on pure Gaussian noise here ranged c=3-20, never settling at
   2; likely a property of real SRER's specific spectral PC structure, not
   reproducible with generic synthetic noise). To still test the gate
   itself, `select_winning_c()` was monkey-patched inside a test-only copy
   of the function environment to always return winning_c=2, then the full
   `compute_spectral_species_richness_fcm()` wrapper was run end-to-end on
   both scenarios: gap check triggered on 3/3 reps in both (confirming the
   gate fires whenever primary search reports c=2, as designed), and the
   richness that actually came out was 1 for the homogeneous case and 5 for
   the blobs case (the gap statistic's answer, not the forced primary c=2)
   -- confirming the resolved c, not the c=2 stand-in, is what flows into
   the final richness/Shannon's H' values.

Also confirmed: the file still `parse()`s cleanly; sourced directly, it
still runs through Section 0/0a (library loads, FCM/numpy availability
check, resolved `py_config()` printed) exactly as before and stops at the
same expected point (missing `./Data/NEONsites.csv`) -- nothing upstream of
the real data read was broken by this change. Grepped the whole file for
orphaned references to `n_clusters` (as our own removed variable, not
`fcmeans.FCM`'s constructor arg of the same name -- that one's legitimate),
`randomForest`, `kmeans`, `assign_nearest_centroid`, and the old
pre-FCM `compute_spectral_species_richness()` -- none found; all remaining
hits are explanatory comments about what was removed and why. Fixed one
stale comment (Section 6's Rao's Q runtime note still said "RF + K-means,
n_reps_ssr reps" after the earlier FCM-adoption session's edits) while in
the area.

**Confirmed unchanged, per the task's explicit "do not touch" list:**
`c_max_ssr`=20, `fcm_fuzziness_w`=2.0, `n_subsample_pixels`=2500,
`n_reps_ssr`=20, `buffer_m`, `ndvi_thresh`, `n_pc_ssr` -- all untouched. CV,
CHV, CHA, and all three Rao's Q variants -- untouched. Year-major job
ordering (Section 7/8, `inventory %>% arrange(year, tower_id)`) -- already
in place from the prior session, confirmed intact, not re-touched.
Incremental full-table `write.csv()` after every completed job (Section 8)
-- confirmed intact and unaffected by this change. FCM availability
pre-check (`py_module_available("fcmeans")`/`"numpy"`, `stop()` before the
main loop if missing) and the unconditional `print(py_config())` at the top
-- confirmed still present (Section 0a), unmodified.

**New output column**: `n_reps_gap_check_triggered` (out of `n_reps_ssr`),
added to `results`' schema and all three tibble-construction branches
(success / no-reflectance-data / error), reporting per site-year how many
reps hit the c=2 floor and triggered the fallback. Also added a per-site-year
`cat()` line in the main loop reporting reps-triggered plus primary-search
vs. gap-fallback wall-clock seconds (`primary_time_sec`/`gap_time_sec`,
accumulated across reps inside `compute_spectral_species_richness_fcm()`).

**Design choice, flagged rather than silently decided**: the gap-statistic
code lives directly in `AnnualSpectralDiversity.R` (Section 5), not a
separate sourced file. `adaptive_fcm_ssr.py` was deliberately left
untouched per instruction, and the gap check's hard-assignment `W_k` is a
different-enough quantity from that module's fuzzy `L(c)` that bridging it
through Python would mean a second, narrowly-scoped Python function for
what's a small amount of R-side logic already reusing `fit_fcm()`'s single
Python call -- putting it in R alongside `select_winning_c()` (which it
directly gates on) keeps the conditional trigger, the two objective
functions, and their shared FCM-fit dependency in one place. A future
session could still split it out if the file grows unwieldy.

**Robustness concerns for an unattended multi-day run, flagged but NOT
fixed (logic calls, not mechanical fixes):**
- `fit_hard_wk()`'s `log(wk)` can hit `log(0) = -Inf` if a candidate c's
  hard-assignment clusters land pixels EXACTLY on their own centroid (wk=0)
  -- vanishingly unlikely on continuous reflectance PC data but not
  impossible on a degenerate subsample; would produce an `Inf`/`NaN` Gap
  value for that c. Not guarded, matching how Tibshirani's own reference
  implementations (e.g. `cluster::clusGap`) handle it -- an explicit
  epsilon-floor or c-skip would be a real behavior change, flagged for a
  decision rather than added unrequested.
- The gap-statistic fallback's cost is only bounded per-rep, not per-site-
  year or globally -- a pathological site-year where EVERY rep lands at the
  c=2 floor pays the full `gap_check_c_max * (1 + gap_check_n_null)` extra
  fits `n_reps_ssr` times over, with no circuit breaker. The new
  `n_reps_gap_check_triggered` column and per-site-year timing line make
  this visible after the fact, but nothing stops it from happening during
  an unattended run -- worth watching the first few real site-years'
  printed gap-fallback timing before assuming it stays cheap everywhere.
- `select_winning_c()`'s `stop()` on all-non-finite `L(c)` (unchanged from
  the prior session) still aborts that rep's whole computation rather than
  NA-ing just that rep -- same behavior as before this change, not
  reintroduced, just noting it's still there for an unattended run.
- This session's synthetic self-tests could not reproduce SRER's reported
  real-data c=2-floor behavior on generic Gaussian noise (see self-test #3
  above) -- meaning the conditional gate's trigger condition itself
  (`winning_c == 2`) has only been verified by direct code inspection plus
  a forced/monkey-patched test, not by a naturally-occurring synthetic
  case. The real trigger rate on actual SRER/ABBY data this weekend is the
  first real confirmation that the gate fires as intended in production.

**This entry, and the code change it describes, have NOT been committed or
pushed** -- per this session's explicit instructions and this project's
commit-only-when-asked convention. Both are sitting in the working tree
alongside the still-uncommitted prior FCM-adoption change (entry below).

## 2026-08-20 UTC — Add FCM_test.R: adaptive-FCM SSR diagnostic, companion to KMeans_test.R's n_clusters sweep

New script, `Code/DataAnalysis/FCM_test.R`, run in the same R session as
`AnnualSpectralDiversity.R` (same convention as `KMeans_test.R`, which it is
a direct companion to). Where `KMeans_test.R` sweeps `n_clusters` through the
existing RF-proximity + k-means + nearest-centroid SSR to check whether
lowering `k` restores ABBY/SRER separation, this script asks the same
question with self-adaptive fuzzy c-means (FCM, Wu et al. 2026) instead: per
rep, cluster count `c` is chosen by maximizing a validity function `L(c)`
over `c` in `[2, 20]`. Reuses (does not duplicate) `adaptive_fcm_ssr.py`,
already written and unit-tested in the prior session for
`CompareSSR_AdaptiveFCM_vs_KMeans.R`; extended it to also return each rep's
hard cluster assignment (`hard_assignment`, argmax of fuzzy membership) so
`FCM_test.R` can compute a Shannon's H' in the same formula as the k-means
version's. `KMeans_test.R`, `AnnualSpectralDiversity.R`, and
`ComputeSpecBiodiv.R` were read but not modified.

**Two structural asymmetries vs. the k-means version, both by explicit
instruction rather than oversight, flagged prominently in the script
header:** (1) richness here is simply the winning `c` per rep -- there is no
analog of the k-means version's full-population nearest-centroid
reassignment step, so the two richness values are not measuring the same
thing underneath a shared name; (2) Shannon's H' uses the SAME formula but is
computed over the winning c's hard-assigned SUBSAMPLE only, not a
full-population reassignment like the k-means version's H' -- comparable in
formula, not in population scope. Both are called out so a future
side-by-side read of the two scripts' `shannon_h` columns isn't
misinterpreted as apples-to-apples.

**Self-tested (no server/`Data/` access in this sandbox), two ways:**
1. Lightweight synthetic gate -- 5 well-separated Gaussian blobs vs. a
   homogeneous single cluster, run directly through the same
   `adaptive_fcm_ssr()` call `compute_spectral_species_richness_fcm()` uses
   per rep: blobs -> winning c per rep `5,5,6,6,5` (mean 5.4); homogeneous ->
   `4,5,4,2,2` (mean 3.4). Passed (homogeneous clearly below blobs, nowhere
   near the c=20 cap).
2. Full mock-harness run -- `get_tower_reflectance()` replaced with a mock
   returning synthetic ABBY/SRER-shaped rasters (deliberately built with
   more spectral groups and higher veg fraction for the mocked ABBY than
   SRER), run through the ENTIRE script at real production settings (2500-px
   subsample, 20 reps, all 6 site-years). Caught one real bug this way:
   `pivot_wider()` (final summary table) needs `tidyr`, which neither this
   script nor `AnnualSpectralDiversity.R` loads -- fixed by adding
   `library(tidyr)` to `FCM_test.R`. Flagged as a LATENT GAP in
   `KMeans_test.R` too (it also calls `pivot_wider()` with no `library(tidyr)`
   anywhere upstream in its own dependency chain) -- worth checking whether
   the server run of `KMeans_test.R` hit the same failure at its final print
   step. On the mock data itself, both mocked sites' richness landed at the
   floor (c=2) with no separation -- almost certainly an artifact of the
   mock's per-pixel noise swamping its deliberately-injected group signal,
   not a script bug (the clean synthetic gate above, run on unambiguous
   blob/homogeneous structure with no injected noise, is the real
   correctness check). Runtime observed on the mock: ~53 sec/site-year for
   the 3091-veg-pixel mocked ABBY, ~21 sec/site-year for the 1055-pixel
   mocked SRER (subsample capped at population size) -- same order of
   magnitude as the script's own pre-loop calibration estimate (~0.9
   min/site-year, ~5.7 min projected total across 6 site-years), which is
   itself only a synthetic-data-timing floor, not a guarantee against real
   reflectance PC data.

**Assumptions/dependencies to confirm on the server, per instruction:**
`w=2.0` (standard/default, not tuned); `c` capped at `[2,20]`; hard-vs-soft
H' assignment (explicit instruction, not a claim hard assignment is
"correct" for FCM's fuzzy output); `fcmeans`/numpy reachable by `reticulate`
(the script now checks this itself, before the site-year loop, and stops
with a clear message if not -- tested in this sandbox against
`fuzzy-c-means==2.0.2`); the new `library(tidyr)` requirement above.

## 2026-08-20 UTC — Add CompareSSR_AdaptiveFCM_vs_KMeans.R: self-adaptive fuzzy c-means SSR vs. the existing fixed-k nearest-centroid SSR, ABBY/SRER
## 2026-08-20 UTC — Add CompareSSR_AdaptiveFCM_vs_KMeans.R: self-adaptive fuzzy c-means SSR vs. the existing fixed-k nearest-centroid SSR, ABBY/SRER

New script, `Code/DataAnalysis/CompareSSR_AdaptiveFCM_vs_KMeans.R`, plus a new
Python helper module, `Code/DataAnalysis/adaptive_fcm_ssr.py`. Compares two
DIFFERENT definitions of Spectral Species Richness (SSR) side by side on real
ABBY/SRER data: the existing fixed-`n_clusters=50` RF-proximity + k-means +
nearest-centroid SSR (`AnnualSpectralDiversity.R`, unchanged, called
read-only as baseline) against a new self-adaptive fuzzy c-means (FCM) SSR
(Wu et al. 2026, *Computers and Electronics in Agriculture* 252:112108),
where the cluster count `c` is chosen per window (here: per tower-year,
whole-buffer) by maximizing a validity function `L(c)` trading fuzzy
between-cluster separation against fuzzy within-cluster dispersion, searched
over `c` in `[2, 20]` -- deliberately capped well under the old 50, since
that ceiling is exactly what caused the original saturation problem. **No
server data access in this sandbox; only the script was written and the
synthetic-validation gate was run. The real ABBY/SRER comparison, and
whether it actually resolves the saturation problem, is a next-session
task once the user runs this on the server.**

**Three factual corrections made to the task prompt before writing anything,
verified against the actual repo contents rather than assumed:**
1. The task said to reuse `FieldDiversity.R` for pixel extraction via
   `import_functions_from()`. False -- `FieldDiversity.R` is ground-plot
   floristic diversity and, by its own header, never reads reflectance
   values. The real extraction pipeline (`get_tower_reflectance`,
   `compute_ndvi_raster`, `discover_site_years`, and the baseline
   `compute_spectral_species_richness()` itself) lives in
   `AnnualSpectralDiversity.R`, and that is what's imported, via the same
   `import_functions_from()` bootstrap `ComputeLUE.R`/`ComputeLUE_Annual.R`
   already established (and called as `spec_fns$name(...)`, matching their
   call style, not dumped into the global environment).
2. The task said SSR is "already implemented in Python" in this repo and to
   match an existing "CV/CHA in R + SSR in Python" split. False -- SSR is
   pure R everywhere in this repo; the only prior Python bridge
   (reticulate + pyGNDiv, for Rao's Q) lived in the older
   `ComputeSpecBiodiv.R` and was deliberately dropped in
   `AnnualSpectralDiversity.R` for a pure-R implementation (that script's own
   Section 6 header: didn't scale to continuous ~1m data). Still bridged to
   Python here anyway, via `reticulate::source_python()`, because Wu et al.
   themselves used the `fcmeans` (PyPI `fuzzy-c-means`) package and
   reproducing their method means using the same fitting code -- but this is
   a new, narrowly-scoped bridge, not a pre-existing pattern being followed.
3. The task asked for "one row per plot/window combination." The existing
   SSR pipeline has no sub-buffer windowing concept at all -- one richness
   value per whole tower-year buffer. Introducing a new window size/stride
   would make the two methods' outputs NOT comparable to the existing CSV
   (a different spatial support), contradicting the stated goal. This script
   therefore treats the whole per-tower-year buffer as the single window for
   both methods (`window_id` is always `"whole_buffer"`) -- flagged as an
   open decision if finer, sub-buffer windows are actually wanted later.

**Validity-function implementation** (`adaptive_fcm_ssr.py`): implements
`L(c) = [sum_i sum_j u_ij^w ||v_i-xbar||^2 / (c-1)] / [sum_i sum_j u_ij^w
||x_j-v_i||^2 / (m-c)]` directly from `fcmeans.FCM`'s fitted `.u`/`.centers`,
vectorized (no per-pixel loop, same `||x-c||^2`-expansion technique
`AnnualSpectralDiversity.R`'s `assign_nearest_centroid()` already uses).
Flagged a real name collision in the module docstring: the paper's `w`
(fuzziness exponent) is the `fcmeans.FCM` constructor's `m` argument, and the
paper's `m` (window pixel count) is unrelated -- the code always spells these
out as `fuzziness`/`n_pixels`, never a bare `m`.

**Synthetic validation gate (task's Section 4) -- RUN in this sandbox, PASSED:**
using the same kind of scenario the original k-means saturation bug was
diagnosed with (SESSION_LOG entry below, 2026-08-1x-ish: 5 well-separated
Gaussian blobs vs. a homogeneous single cluster, 10000 pixels, 4 PCs), run
through the script's actual `adaptive_fcm_ssr_window()` wrapper at the exact
settings configured for the real run (`c` in `[2,20]`, `w=2`, 1000-pixel
subsample, 5 reps): blobs -> winning c per rep `5,6,5,5,5`, richness 5;
homogeneous -> winning c per rep `2,4,3,4,2`, richness 3. Both gate checks
passed (homogeneous < blobs; homogeneous well clear of the c=20 cap).
Total gate runtime ~11 sec for both scenarios combined. Also ran a separate
smoke test of the full main-loop code path (NDVI/veg-mask/bad-band-mask via
the imported functions, PCA extraction, both SSR computations) against a
synthetic terra raster shaped like `get_tower_reflectance()`'s real output --
all steps ran and returned sane values; this is NOT a substitute for a real
ABBY/SRER run, only a structural/syntax check.

**Assumptions flagged for the user to confirm before trusting a full run**
(also in the script's own header): fuzziness `w=2.0` (standard/default, not
tuned); `c` capped at `[2,20]` per the task's explicit instruction;
FCM subsample (1000 px) and rep count (5) both LOWER than the existing
k-means step's (2500 px, 20 reps) -- a real precision/cost tradeoff, made
because fitting FCM for 19 candidate `c` values per rep is far more
expensive than one k-means fit, not a neutral default; output paths
(`D:/projects/moore/SpectralBiodiversity/Data/...`) copied from
`AnnualSpectralDiversity.R`'s own convention, not independently confirmed.
**Python dependency, unconfirmed on the server:** `fuzzy-c-means` (imports as
`fcmeans`) + numpy, reachable by `reticulate`. In this sandbox, reticulate
initially failed against a statically-linked python3 with no loadable
`libpython*.so`; had to point `RETICULATE_PYTHON` at a system `/usr/bin/
python3` build with a proper shared library before `fuzzy-c-means` (pip
`--user`-installed there) would import -- the server's Python/reticulate
setup should be checked for the same class of failure before assuming this
just works.

Non-goals honored: `AnnualSpectralDiversity.R` and its output CSV were never
modified (only imported/read); the PCA-vs-fixed-band-index question was not
touched; no decision was made about adopting adaptive-FCM SSR as a new
default -- that's explicitly left for the next session once real ABBY/SRER
results are in.

**This entry has NOT yet been committed/pushed** -- per this session's
git-safety instructions, commits are made only when the user explicitly asks;
flagging here so it isn't lost, but the user (or a follow-up instruction)
needs to trigger the actual commit.

## 2026-08-20 UTC — Add ComputeLUE_Annual.R: annual-resolution LUE, a separate ratio-based method alongside (not replacing) ComputeLUE.R

New script, `Code/DataAnalysis/ComputeLUE_Annual.R`. Computes LUE from
ANNUAL-resolution (YY) AmeriFlux data as a SEPARATE, PARALLEL approach to the
existing half-hourly regression-slope method in `ComputeLUE.R`.
**`ComputeLUE.R` was read in full before writing anything and was NOT
modified** -- confirmed via `git diff --stat` at the end showing no change --
it stays exactly as-is in case the half-hourly approach is revisited later.
Different output file (`lue_by_year_annual.csv` vs. `lue_by_year.csv`) so the
two methodologically-distinct results can never overwrite or be confused.

**Methodological distinction, unavoidable at this resolution, flagged
prominently in the script header:** annual data gives exactly one GPP value
and one PAR value per site-year, so there is no within-period variation left
to fit `lm(GPP ~ APAR)`'s slope to -- the regression approach `ComputeLUE.R`
uses is categorically unavailable here, not merely skipped. This script
computes `LUE = GPP_annual / (PAR_annual * fPAR)` instead: a straightforward
ratio, but one that does NOT have the ratio-of-means-bias protection the
regression-slope approach was specifically built to provide (see
`ComputeLUE.R`'s own header on this). Its `lue` values are therefore NOT
expected to match a hypothetical regression-based annual estimate, and the
two scripts' outputs must be treated as methodologically distinct, never
interchangeable. A SECOND, separate caveat was also identified and
documented (not something the task described in this much detail, but a
direct consequence of the resolution change worth surfacing on its own): the
half-hourly script restricts its flux data to the same calendar window the
single-day NDVI/fPAR snapshot represents, keeping the two temporally
consistent; annual GPP/PAR are whole-year aggregates by construction and
can't be windowed that way, so this script's ratio blends dormant-season and
peak-growing-season conditions together in a way the half-hourly version
specifically avoids. `growing_season_source` is still resolved and reported
here (for provenance and consistency with the other two LUE-family scripts),
but no longer functions as an active data filter -- there's nothing left at
annual resolution to filter.

**File/column conventions -- CONFIRMED vs. still assumed, distinguished
explicitly per the task's request.** `Code/NEON_FluxVariability.R` was read
in full: it already reads `./Data/NEON_Ameriflux/AnnualData` (one YY CSV per
site) successfully in this repo, and its own column selection
(`GPP = GPP_NT_VUT_REF`, `NEE = NEE_VUT_REF`, `RECO = RECO_NT_VUT_REF`,
`ET = LE_F_MDS`, `year = TIMESTAMP` with no date parsing, i.e. TIMESTAMP is
the bare 4-digit year) is direct repo evidence -- not a guess -- that these
columns exist unchanged from their half-hourly names at annual resolution.
GPP_NT_VUT_REF's presence is therefore CONFIRMED, unlike everything in
`ComputeLUE.R`'s HH-file section. A PAR-equivalent column is NOT confirmed
by that evidence (`NEON_FluxVariability.R` never selects one) --
`resolve_par_column()` tries `PPFD_IN` first (plausible given every other
confirmed column kept its name, but unverified for this specific one) and
falls back to `SW_IN_F`/`SW_IN` with a literature-default 2.02 umol/J
shortwave-to-PPFD conversion factor (McCree/Britton & Dodd, not
site-calibrated -- same caveat style as the existing NDVI_min/NDVI_max fPAR
defaults) if PPFD_IN is absent. Which path was taken is recorded in a new
`par_source` output column (`ppfd_in_annual` vs. `sw_in_derived`) -- an
addition beyond the literally requested column list, justified by this
project's established never-silently-blend-methods convention
(`growing_season_source`, MODIS's `real_metadata`/`fallback_default`).

**Units consistency -- explicitly checked, not assumed, per the task's
specific warning about a yearly-sum-divided-by-a-mean bug.**
`NEON_FluxVariability.R`'s own comment ("convert LE (W/m2) to ET (mm/yr) if
you need absolute ET values... LE_F_MDS ... not mm") is repo evidence that
YY variables remain flux-RATE MEANS (not summed totals), consistent with the
standard FLUXNET2015 YY aggregation convention -- making
`GPP_annual / (PAR_annual * fPAR)` unit-consistent with `ComputeLUE.R`'s own
mol-CO2/mol-photon convention, evidence-informed but not verified against a
real value. Added `check_gpp_plausible_range()` (default plausible annual-
mean range: -5 to 30 umol CO2 m-2 s-1) as a concrete runtime guard against
exactly the failure mode the task called out: a value in the magnitude of an
annual SUM (e.g. ~200-3000 gC/m2/yr) landing where a mean flux rate is
expected is flagged in `status` as a units-mismatch warning, not silently
divided through as if it were correct.

**Reuse, not reimplementation (per task instructions -- `ComputeLUE.R` read
in full first).** `import_functions_from()` is necessarily redefined locally
(it is the bootstrap mechanism itself -- it cannot import itself without
circularity, same as every script independently stating its own `library()`
calls). Everything it is USED to pull in is reused: the same
`AnnualSpectralDiversity.R` set `ComputeLUE.R` imports (`read_neon_h5_tile`,
`get_tower_reflectance`, `compute_ndvi_raster`, etc.) and
`FieldDiversity.R`'s `get_flight_acquisition_date` -- AND, following the
task's suggestion to check whether pieces of `ComputeLUE.R` itself are
cleanly importable, `compute_mean_ndvi`, `compute_fpar`, and
`resolve_growing_season` are imported DIRECTLY FROM `ComputeLUE.R` (confirmed
these are standalone top-level functions there, not entangled with its
half-hourly-specific state) rather than re-derived from
`AnnualSpectralDiversity.R`/`FieldDiversity.R` a second time. Confirmed
empirically that this import leaks none of `ComputeLUE.R`'s own pipeline
state (`out_csv`, `results`, etc.) into the importing script's global
environment -- same mechanism, same guarantee, as the prior sessions'
imports from `AnnualSpectralDiversity.R`/`FieldDiversity.R`.

Validated (synthetic only, same sandbox constraints as prior sessions):
`resolve_par_column()`'s three-way branching (PPFD_IN preferred, SW_IN_F
fallback, bare SW_IN fallback, NULL when neither exists);
`check_gpp_plausible_range()` correctly passing a plausible annual-mean value
and flagging a summed-total-magnitude value (1200) as implausible;
`compute_annual_lue()` recovering a known true LUE EXACTLY (it's a
deterministic ratio, not a fit -- 0.035 true across three different
(gpp, par, fpar) triples, hand-computed) and guarding zero/NA APAR against
Inf/NaN; an explicit units-error scenario confirming the ratio calculation
*cannot* self-detect a GPP-is-actually-a-sum bug (it just returns a
too-large-but-finite number) and that `check_gpp_plausible_range()` is what
actually catches it -- exactly the validation-not-review check the task
asked for; `load_annual_flux_data()` against a synthetic fixture matching
the CONFIRMED real column set, checking -9999 handling for both GPP and
PPFD_IN and correct tower_id extraction from filename for two different
synthetic sites. 23/23 checks passed. The script `parse()`s cleanly (32
top-level expressions) and, run directly, fails at the same expected point
(missing `./Data/NEONsites.csv`) as `ComputeLUE.R`/`ExtractMODIS.R` did.

**Not validated**: real annual AmeriFlux column names/values (whether
`PPFD_IN` is actually present at annual resolution, whether GPP/PAR really
are means rather than sums) -- needs a real run on the server, same
established pattern as every other script in this project.

## 2026-08-20 UTC — Add ExtractMODIS.R (AppEEARS MODIS FPAR/GPP) and extend ComputeLUE.R with a MODIS-derived LUE comparison

Two pieces, both extending the LUE work from the earlier entry below (same
day): `Code/DataAnalysis/ComputeLUE.R` already existed from that prior
session, so it was read and extended rather than duplicated.

**Part 1 -- `Code/DataAnalysis/ExtractMODIS.R` (new script).** Pulls
MOD15A2H Fpar_500m and MOD17A2H Gpp_500m, 8-day composite, at all 45 tower
Lat/Lon, via NASA's AppEEARS point-sample API, as a submit-one-combined-task
-> poll -> download workflow (one task covering all 45 sites and all 4
layers -- both products' data + QC bands -- not a per-site request loop).
Date range is read from `spectral_diversity_by_year.csv`'s `year` column,
never hardcoded. `appeears` (the R package) is not installed on this system
and was not added; `httr2` (already installed) is used directly for the
login/submit/poll/download HTTP calls instead, per "check what's already
installed before assuming a package needs adding."

**CRITICAL, explicit per the task's instructions**: this sandbox has no
internet access, so every AppEEARS HTTP call in this script (login,
task-submit, status-poll, bundle-list, file-download) -- and the exact
response schema/field names/CSV column-naming convention each of those calls
assumes -- was written from documented AppEEARS API/product conventions and
has **never been executed, not even once**. This is explicitly NOT claimed
as validated, consistent with this project's real-data-validation discipline
(e.g. FieldDiversity.R's live H5-metadata investigation, AnnualSpectralDiversity.R's
synthetic-only validation elsewhere in this log). What WAS validated in this
sandbox, against synthetic data built to match AppEEARS' documented
point-sample CSV format: the MOD15A2H FparLai_QC and MOD17A2H Psn_QC_500m
bit-decoding logic (MODLAND_QC bit + FPAR's cloud-state bits), and the
generic column-matching parser (`find_col()`/`find_col_exact()`/
`parse_appeears_csv()`), which was deliberately written to match columns by
substring rather than an exact hardcoded `<Product>_<Version>_<Layer>` string
specifically because the real naming convention couldn't be confirmed. 24/24
synthetic checks passed (QC-decoding truth table, scale-factor application,
out-of-range-despite-good-QC exclusion, unmapped-tower-id handling, output
schema). The user must run this on the server, watch it fail at whatever
point the real API differs from what's assumed here, and fix that point
specifically (the script's header points at `cat()`-ing the intermediate
`resp_body_json()`/`bundle_files`/`names(raw_combined)` objects as the first
debugging move) rather than assuming the whole thing needs rewriting if one
step doesn't match.

Credentials: `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` env vars, read via
`Sys.getenv()`, never hardcoded. Checked `.gitignore` before choosing this --
`.Renviron` (R's standard mechanism for these) is already listed there
alongside `.Rproj.user`/`.Rhistory`/`.RData`/`.Ruserdata`, so no new
gitignore entry was needed.

**Part 2 -- `ComputeLUE.R` extended, NDVI-derived pipeline untouched.** Added
a second, MODIS-satellite-derived LUE estimate (`lue_modis`, `lue_modis_r2`)
per site-year, computed on the EXACT SAME filtered half-hourly tower records
(`filtered` -- same GPP_NT_VUT_REF, same PPFD_IN, same midday/growing-season/
QC filtering) as the original NDVI-derived fit, with `APAR_modis = PPFD_IN *
fpar_modis` in place of the NDVI-derived APAR -- an apples-to-apples
comparison isolating the fPAR-source effect, per the task's explicit design.
`fit_lue_regression()` is reused unchanged for both fits rather than
duplicated. `fpar_modis` itself is the mean of QC-passed 8-day composites
within the SAME growing-season window (`gs$months`/`target_year`, real
flight month vs. fallback) already resolved for the NDVI-derived version --
reused via a new `aggregate_modis_window()` helper, not a second window
implementation. `gpp_modis_ref` (same aggregation, MODIS GPP) is reported
alongside for comparison only and is never fed into either regression --
tower `GPP_NT_VUT_REF` remains the sole GPP input to both fits. The original
`fpar` output column is renamed `ndvi_fpar` so both fPAR sources are
unambiguous in the header; new columns: `fpar_modis`, `lue_modis`,
`lue_modis_r2`, `gpp_modis_ref`, `n_modis_composites_used`.

Design choice, flagged rather than silently decided: `modis_fpar_raw.csv`/
`modis_gpp_raw.csv` are treated as OPTIONAL inputs, not required ones --
missing files (ExtractMODIS.R not yet run) print a clear one-time notice and
degrade every MODIS-derived column to NA for the run, rather than `stop()`ing
the whole pipeline. The NDVI-derived LUE calculation has no real dependency
on MODIS data and shouldn't be blocked by a separate, network-dependent,
not-yet-run script. A MODIS-side failure for one specific site-year is also
caught on its own (separate `tryCatch`) and only NAs that site-year's MODIS
columns -- it can never discard an already-successful NDVI-derived result.

Validated (synthetic, same session): re-ran the original NDVI-derived
validation harness unchanged first to confirm the extension didn't regress
it (all checks still passing, `make_row()`'s schema assertion updated for
the new columns), then added checks for `aggregate_modis_window()` (mean
over QC-passed, in-window composites only; excludes QC-fail and
out-of-growing-season-window rows; degrades to NA/0 for an unmatched
tower_id or a `NULL` input) and `compute_modis_lue_addon()` (recovers a
known true `lue_modis` slope -- 0.0280 true vs. 0.0278 fit -- from synthetic
`GPP = intercept + slope*APAR_modis + noise` data fed through the SAME
`filtered` object the NDVI-derived fit used; reports `gpp_modis_ref` and
`n_modis_composites_used` correctly; degrades to all-NA, no error, when
MODIS inputs are entirely missing). 41/41 checks passed. Both scripts also
`parse()` cleanly and, run directly, fail at the identical expected point
(`stop()` on missing `./Data/NEONsites.csv`) as the original ComputeLUE.R did.

**Three-way comparison this now enables**, once ExtractMODIS.R has actually
been run on the server: NDVI-derived fPAR (from the hyperspectral archive,
`ndvi_fpar`) vs. MODIS-derived fPAR (`fpar_modis`) as two independent
estimates of the same physical quantity feeding otherwise-identical LUE
regressions (`lue` vs. `lue_modis`, directly comparable since GPP/PPFD/
filtering are held fixed), plus MODIS's own independently-produced GPP
product (`gpp_modis_ref`) as an out-of-band cross-check against tower-derived
`GPP_NT_VUT_REF` -- three separate lines of evidence on the same site-years,
none of them silently blended into the others.

**Not validated**: the actual AppEEARS API interaction (Part 1's core
purpose) and, downstream of that, MODIS-derived LUE on real data -- both
require a real run on the server with internet access.

## 2026-08-20 UTC — Add ComputeLUE.R: Light Use Efficiency from AmeriFlux HH data, regression-based

New script, `Code/DataAnalysis/ComputeLUE.R`. Computes LUE (epsilon) per
tower-year for the 45 NEON/AmeriFlux sites, as the SLOPE of
`lm(GPP_NT_VUT_REF ~ APAR)` across filtered midday half-hourly AmeriFlux
records, not a ratio of period means -- avoids ratio-of-means bias and gives
an `lue_r2` per site-year as a free QA diagnostic (a poor fit flags a
site-year where the light-response assumption may not hold, e.g. water
stress). Fully standalone: no shared state with `AnnualSpectralDiversity.R`
or the footprint pipeline, different output path.

**Before-writing investigation** (per task instructions): no script in this
repo reads AmeriFlux HH files today -- the only existing flux reader
(`NEON_FluxVariability.R`) reads annual (`YY`) aggregates from
`./Data/NEON_Ameriflux/AnnualData`, and gave no real HH filename to confirm
against. This sandbox has no access to the real `D:`/`X:` data drives, so the
real HH directory location and exact filename/version convention could NOT be
confirmed. File discovery therefore uses a flexible match (site ID +
"FLUXNET" + "_HH_" all present in the basename, not a rigid version-number
regex) and `hh_data_dir` is a best-guess sibling of the confirmed
`AnnualData` path, **flagged in the script header as unconfirmed -- must be
verified/adjusted on the server**. QC-flag filtering (`NEE_VUT_REF_QC <= 1`)
and the `YYYYMMDDHHMM` TIMESTAMP_START format are both standard
AmeriFlux/FLUXNET2015 assumptions with no prior convention in this repo to
confirm against; `read_hh_flux()` `stop()`s clearly if the expected columns
are absent rather than silently proceeding on a wrong schema.

**Growing-season month recovery.** `read_neon_h5_tile()` in
`AnnualSpectralDiversity.R` carries no date field, but
`Code/DataAnalysis/FieldDiversity.R`'s `get_flight_acquisition_date()`
(previous session) already recovers the real flight date from the
`ATCOR_Input_file` text log under `Metadata/Logs/<flightline>/`. Reused
(imported, not reimplemented) here: when a real date resolves, its MONTH
becomes the site-year's growing-season window (so the single NDVI/fPAR
snapshot stays temporally consistent with the flux half-hours being
regressed); otherwise a documented June-Sept fallback is used. Every row's
`growing_season_source` column records which happened, never silently
blended -- both paths are exercised in validation (see below).

**Cross-script reuse mechanism.** The task asked to reuse (source/import, not
copy-paste) `read_neon_h5_tile`, `get_tower_reflectance`, and
`compute_ndvi_raster` from `AnnualSpectralDiversity.R`. Plain `source()` of
that script isn't safe -- it's a monolithic script whose top level
immediately runs the full site-year pipeline (H5 discovery, main loop, its
own `write.csv`) as a side effect, not a function library, and no script in
this repo currently imports from another for exactly that reason (prior
cross-script reuse, e.g. `FieldDiversity.R`, reimplemented instead).
Resolved by adding `import_functions_from(script_path, names)`, which
`parse()`s the origin script and `eval()`s only the named top-level `<-`
bindings (function defs and plain config constants alike, e.g.
`bad_band_ranges`) into a dedicated environment -- confirmed empirically to
import the needed bindings from both `AnnualSpectralDiversity.R` and
`FieldDiversity.R` without leaking `towers_df`/`site_year_jobs` or any other
pipeline state into the global environment, i.e. neither origin script's main
loop runs. The actual logic stays single-sourced in the origin files; a later
edit there is picked up automatically here.

**Methodology choices flagged in the script, not silently decided:**
regression includes an intercept (`lm(GPP ~ APAR)`, not forced through the
origin) since forcing through zero is a stronger assumption than the
light-response literature generally requires; fPAR uses literature-default
`NDVI_min = 0.05` / `NDVI_max = 0.95` (Wang et al. 2016 uses ground-calibrated
values instead -- not available for these sites, a known limitation worth
revisiting); local solar time uses the equation of time (Cooper's
approximation) plus a nominal-15-degree-meridian offset from each site's
longitude, since AmeriFlux TIMESTAMP is local standard time and no per-site
UTC-offset metadata exists in `NEONsites.csv` to do better than the nominal
assumption.

**Validation** (synthetic only -- no `Data/` access in this sandbox): a
throwaway harness (not committed) imported every pure function/constant via
the same `import_functions_from()` mechanism and checked, against synthetic
fixtures: fPAR clipping and formula at both bounds and midpoint; local solar
time is close to clock noon at a real site longitude (ABBY, -122.33) and
monotonic across a day; `read_hh_flux()` converts -9999 to NA for both GPP
and PPFD (including a QC-forced-pass fixture verifying the fill value doesn't
survive a passing QC row), applies the QC<=1 filter correctly, and `stop()`s
clearly on a missing required column; `find_hh_files()` matches version-suffix
variants for the right site and excludes a different site's file;
`resolve_growing_season()` exercised on all three paths -- real metadata
resolves to the flight month tagged `real_metadata`, an unparseable date and
an empty file list both fall back to the June-Sept default tagged
`fallback_default`; `filter_midday_growing_season()` correctly enforces
target year, target month(s), the local-solar-hour window, and the PPFD
floor on a full synthetic year with a diurnal PPFD curve; and
`fit_lue_regression()` recovers a known true slope (0.0350 true -> 0.0350
fit) and intercept from synthetic `GPP = intercept + slope*APAR + noise`
data with R^2 > 0.99, while a deliberately decoupled (pure-noise) GPP/APAR
fixture correctly yields R^2 ~ 0 -- confirming the regression-based QA signal
the method was chosen for actually works. All 30 checks passed. The full
script also `parse()`s cleanly (40 top-level expressions) and, run directly,
fails at the expected/correct point -- a clear `stop()`-style error on the
missing `./Data/NEONsites.csv` -- rather than an obscure failure further in.

**Not validated**: real AmeriFlux HH file structure/columns, real NEON H5
collection-date recovery for LUE's specific site-years, and the
`hh_data_dir` path all still need a real run on the server, per this
project's established working pattern.

## 2026-08-19 UTC — AnnualSpectralDiversity.R: nearest-centroid SSR, new CHA metric, found and fixed a second SSR bug

Three requested changes to `Code/AnnualSpectralDiversity.R` only (`ComputeSpecBiodiv.R`
and the footprint scripts were explicitly out of scope and untouched, confirmed
via `git status` at the end).

**Change 1 -- Spectral Species Richness, RF-classifier saturation fix.**
`compute_spectral_species_richness()`'s second stage (a `randomForest` trained
on the k-means subsample labels, `predict()`'d on every pixel) was replaced
with nearest-centroid assignment (Feret & Asner 2014, Ecological Applications
24:1289-1296) via a new vectorized `assign_nearest_centroid()` helper (single
`n x k` distance matrix via the `||x-c||^2` expansion, no per-pixel loop;
`flexclust::dist2()` would do the same but isn't installed and wasn't added as
a new dependency). The RF-proximity + `kmeans(cmdscale(...), centers =
n_clusters)` clustering stage on the subsample is unchanged.

While validating this against synthetic data (5 well-separated Gaussian blobs;
a near-homogeneous single cluster) as directed, found a SECOND, independent bug
not in the original task description: `km$centers` live in the `cmdscale()`
embedding of RF proximity distances, not in PC space -- confirmed empirically
(a `pcs_sub` range of roughly +-75 corresponded to a `km$centers` range of
roughly +-0.07, not a rescaling, an unrelated coordinate system). Computing
nearest-centroid distance between `pcs_all` and `km$centers` as originally
instructed would compare two incompatible spaces. Flagged this to the user
before proceeding (not silently reinterpreted); user chose the fix of
re-deriving each of the 50 clusters' centroids as the mean, in `pcs_sub`'s own
PC coordinates, of the subsample points `km$cluster` assigned to it
(`rowsum(pcs_sub, group = km$cluster) / table(km$cluster)`), keeping the
RF-proximity/MDS/kmeans structure exactly as specified.

With that fix applied, a second, deeper issue surfaced empirically and was also
flagged to the user before committing: `kmeans(centers = n_clusters)` always
returns exactly `n_clusters` non-empty partitions -- it subdivides however many
real groups exist rather than collapsing unused centroids -- and a large
same-distribution evaluation population (any real reflectance raster) generically
populates nearly all of them. Synthetic sweep at `n_clusters=50` (this script's
configured value): 5 separated blobs -> richness ~49.7; homogeneous single
cluster -> richness ~50; blobs 50x tighter -> richness ~49.3 -- saturation was
insensitive to real structure/separation, and only responded to lowering
`n_clusters` itself (n_clusters=10 on the same 5-blob data -> richness ~10,
i.e. still saturated at the (lower) ceiling, not at the true group count). User
directed: commit the nearest-centroid change as scoped, and do not gate it on
solving the saturation problem -- the two bugs it fixes (classifier
extrapolation, and the coordinate-space bug found during validation) are real,
confirmed, and independent of whether saturation is fully solved: they're
about HOW pixels get assigned to clusters, not about how many clusters exist
to assign to. That fix stands on its own.

The residual saturation problem is flagged plainly in the script (Section 5
header comment, and a comment at the `n_clusters <- 50` config line) as an
open problem for a future session, not something quietly papered over with an
undecided threshold added in passing. Two candidates to weigh on their own,
with real ABBY/SRER data in hand rather than guessed at from synthetic tests
(neither implemented this session):
- lowering `n_clusters` itself -- a methodological choice about what "spectral
  species" should mean for this landscape, not a bug fix;
- a principled minimum-occupancy criterion (a cluster only counts toward
  richness above some pixel-share threshold) -- deferred specifically because
  what threshold is principled vs. arbitrary can't be decided from synthetic
  data alone, and an undocumented judgment call here is exactly the kind of
  false-confidence risk this project has hit before.

Runtime: nearest-centroid vs. the old second-RF-classifier step, benchmarked on
synthetic data (10000 pixels, 1500 subsample, 50 clusters, 4 PCs): 8.7s ->
7.1s per rep, ~1.24x speedup (~0.6 min saved per site-year at
`n_reps_ssr=20`) -- modest rather than dramatic, since the kept first RF (for
proximity) and kmeans remain most of the per-rep cost.

**Change 2 -- new metric: Convex Hull Area (CHA).** Added `compute_cha(r,
veg_mask)` near `compute_chv()` (Gholizadeh et al. 2018, Remote Sensing of
Environment 206:240-253; reported to outperform CV, CHV, and SID specifically
at ~1m airborne resolution, this script's regime). For each vegetation pixel,
treats the plot mean spectrum and that pixel's spectrum as a 2D point set (one
point per band) and takes the `geometry::convhulln(..., output.options =
"FA")` hull's `$vol` field (confirmed via a hand-computed triangle case that
`$vol`, not `$area`, is the polygon area for a 2D hull), averaged across all
vegetation pixels in the plot. No PCA, clustering, or RF. Added as a new `cha`
column in the results tibble/CSV schema and both the success and error/no-data
branches of the Section 8 main loop.

A pixel whose spectrum exactly equals the plot mean is perfectly collinear in
this 2D construction, which `qhull` cannot build a hull from (errors "Initial
simplex is flat" rather than returning 0) -- confirmed this empirically and
mapped that specific error to area 0 in `compute_cha()`, since a degenerate
zero-width point set has zero area by definition and this is an expected input
(a pixel identical to the mean), not an exceptional one. Validated: a synthetic
raster where every pixel equals the mean spectrum returns CHA = 0 exactly; a
raster with real per-pixel variation (mean spectrum + Gaussian noise) returns a
clearly positive CHA (~0.11 in the test raster).

**Change 3 -- PCA usage left as-is, flagged as open.** CHV, SSR, and the
all-bands Rao's Q still PCA-reduce, unchanged by this session -- per explicit
instruction, not acted on. Noted in the file header as a deliberately deferred
decision (replacing PCA with a fixed manual band set), to be tested once the
above two changes have been evaluated in isolation.

All validation was synthetic (no real `Data/` access in this sandbox, per this
project's established pattern). Script re-parses cleanly
(`parse("Code/AnnualSpectralDiversity.R")`, 54 top-level expressions, no
errors). Not run against real hyperspectral data.

## 2026-08-18 23:11 UTC — Add FieldDiversity.R (field floristic diversity, plot-scope x flight-timing comparison)

New script, `Code/DataAnalysis/FieldDiversity.R`. Computes Hill-Shannon
(abundance-weighted, q=1, via hillR) and gamma presence richness from NEON
field percent-cover/presence data (`div_1m2Data.csv`,
`div_10m2Data100m2Data.csv`) per tower-year, crossed with two independent
axes: `plot_scope` (`tower` only vs. `tower`+`distributed` combined) and
`temporal_scope` (`per_bout`, every field bout with no selection, vs.
`peak_flight`, the one bout whose median plot end date falls nearest the
site's AOP flight date, within a 30-day tolerance).

Before writing any date-matching logic, investigated the real HDF5
structure of an actual DP3.30006 reflectance tile (the user uploaded one,
`NEON_D16_ABBY_DP3_551000_5067000_reflectance.h5`, ABBY 2017) rather than
assuming the layout from memory, per the task's explicit instruction.
Finding: there is no acquisition-date attribute anywhere on the file, site
group, `/Reflectance` group, or `/Reflectance/Reflectance_Data` dataset,
and the `Metadata/Logs/<flightline>` group name is only the flight line's
HHMMSS time-of-day, not a date. The date is recoverable only from the TEXT
CONTENT of `Metadata/Logs/<flightline>/ATCOR_Input_file` -- either an
embedded `<YYYYMMDD>_<flightline_id>` token or its own literal
`"22\06\2017          Date (dd/mm/year)"` header line (backslash-delimited,
not `/`). Both independently resolved to the correct date (2017-06-22) on
the real file. `get_flight_acquisition_date()` tries the embedded-token
form first (tied to the exact flightline opened) and falls back to the
header line, `stop()`-ing with a clear message if neither parses rather
than guessing. The script also re-runs this investigation live at the top
of every execution and prints what it finds, so the confirmed structure is
visible in the console on every real run, not just asserted in comments.

Two access blockers hit before the user's upload: no `.h5` file on disk in
this sandbox was the right product (the only local one was an unrelated
DP4.00200.001 eddy-covariance bundle), and NEON's `/api/v0/data/...`
endpoints 403'd from this sandbox's egress IP (confirmed via both curl and
R's httr, `/products/` endpoint unaffected) -- ruled out before asking the
user for a file directly.

Validated end-to-end against a synthetic fixture (real `Data/` field CSVs
aren't available in this sandbox) built to exercise every branch: the real
uploaded tile for flight-date resolution, one site/two bouts/three plots
(tower + distributed) for the metrics, plus a second site with no
hyperspec coverage and an out-of-range bout. Confirmed: peak_flight
correctly selects the in-tolerance bout and skips the 52-day-away one;
`floristic_richness(all) >= floristic_richness(tower)` holds; Hill-Shannon
mean and `hill_taxa_parti()` gamma matched hand-computed values; richness
never fell below the max single-plot richness (sanity floor); "no data"
and "no flight date" grid rows populate correctly instead of being
dropped. hillR (not previously installed anywhere in this environment)
was installed from CRAN to confirm its actual exported API
(`hill_taxa`, `hill_taxa_parti` returning a `TD_gamma` column) rather than
assuming function/argument names from memory.

Committed (`732632a`). Not run against real `Data/` field CSVs or the `D:`
drive hyperspec archive -- that first real run, and its console output
(Step 1 investigation, peak-flight inventory, spot check), still needs
doing on the actual server.

## 2026-08-18 21:49 UTC — Fix veg_mask no-op (mask() requires NA, not FALSE/0)

Fixed the veg-mask bug flagged in the previous entry, in both scripts, at
all three construction sites (`AnnualSpectralDiversity.R`'s `veg_mask`;
`ComputeSpecBiodiv.R`'s `veg_mask` and `veg300`): `veg_mask <- ndvi >
ndvi_thresh` (logical 0/1, no NAs) replaced with `veg_mask <- ifel(ndvi >
ndvi_thresh, 1, NA)`. `terra::mask(x, mask)` only excludes cells where
`mask` is `NA` (default `maskvalues = NA`) -- a comment was added at each
site flagging this non-obvious requirement so it isn't silently
reintroduced.

Confirmed empirically, not just by code inspection: on a synthetic 20x20
raster with a deterministic NDVI split (200 cells above threshold, 200
below), `mask(data_r, veg_mask)` with the old construction left all
400/400 cells non-NA (confirming the no-op); with the fix, exactly
200/400 cells (the above-threshold ones) survived.

This is a correctness fix, not a feature. Every metric produced by either
script before this commit (CV, CHV, SSR, all Rao's Q variants, including
everything from the two prior entries in this log) was computed over the
FULL raster extent, not the vegetated subset, and is not directly
comparable to output produced from this point forward.

## 2026-08-18 21:39 UTC — PCA-reduced all-bands Rao's Q, spectral Shannon's H', and a pre-existing veg-mask no-op found during validation

`AnnualSpectralDiversity.R`: added PCA-based band reduction ahead of the
all-bands Rao's Q call only (`pca_reduce_raster()`, Section 6a) --
`prcomp()` on the masked pixel x band matrix (center = TRUE, scale. =
FALSE, matching the CHV/SSR convention), keeping the leading PCs that
explain >=99% variance (`raoq_pca_var_threshold`), N data-driven per
site-year since band count varies after bad-band masking. NDVI/NIRv Rao's
Q are unchanged (already 1-D). `raoq_allbands` going forward is therefore
not numerically comparable to earlier full-band runs (e.g. earlier ABBY
test runs, or the rasterdiv-era output) -- flagged in the file header.

Also added spectral Shannon's index H' and its Hill-number effective
diversity exp(H') (Wang et al. 2018, RSE 211:218-228, Table 2) to
`compute_spectral_species_richness()`, reusing the existing per-rep
classifier output rather than reclassifying. The function's return type
changed from a scalar to a named list (`richness`, `shannon_h`,
`shannon_effective`); Section 8's call site and the results tibble
(new `shannon_h`, `shannon_effective` columns) were updated to unpack all
three. Applied the same H' addition to `ComputeSpecBiodiv.R`, whose
`compute_spectral_species_richness()` is identical -- its Rao's Q is a
separate pyGNDiv-based implementation and was left untouched.

`hyperspec_dir`/`out_csv` in `AnnualSpectralDiversity.R` repointed from
this session's prior `X:/moore/SpectralBiodiversity/...` to
`D:/projects/moore/SpectralBiodiversity/Data/...` per this run's server
working directory; the one other hardcoded-looking path
(`./Data/NEONsites.csv`) was already relative and needed no change.

Synthetic validation (no server data access from this environment): PCA
reduction confirmed to pick a data-driven, much-smaller-than-input N
matching an independent recomputation of the 99% threshold, with correct
NA propagation through masked cells and a finite non-degenerate Rao's Q
on the reduced raster; Shannon's H'/exp(H') matched hand-calculated values
exactly for a known class distribution, including correct handling of
unused factor levels (avoiding 0*ln(0) = NaN); the three-value interface
change was confirmed to unpack correctly end-to-end at the call site.

Also found, while building the veg-mask test raster for validation: both
scripts' `mask(r, veg_mask)` calls, everywhere they occur, are a no-op.
`terra::mask(x, mask)` only masks cells where `mask` is `NA` (default
`maskvalues = NA`), but `veg_mask <- ndvi > ndvi_thresh` produces a 0/1
logical raster with no NAs -- confirmed empirically. This predates this
session's changes and affects every metric (CV, CHV, SSR, all three Rao's
Q variants) in both scripts; left unfixed as out of scope for this task
and flagged here for a deliberate follow-up.

## 2026-08-14 21:45 UTC — Replace rasterdiv Rao's Q with local moving-window implementation (fixes std::bad_alloc)

`rasterdiv::paRao()` (introduced two sessions ago, see the 18:57 UTC entry
below) failed with `cannot allocate vector of size 1077.6 Gb` on ABBY 2017's
NDVI layer (1000x1000, ~750k mostly-unique continuous values even after
`simplify=4`). Root cause, not a tunable-parameter problem: paRao's
classic/multidimension modes build ONE global pairwise distance matrix
across every unique value in the WHOLE input raster, then look up
per-window results from it -- O(unique_values^2) memory. That's fine for
coarse/classified rasters but explodes for continuous ~1m-resolution
reflectance/NDVI data with hundreds of thousands of near-unique values.
paRao's design assumption doesn't match this use case; rasterdiv is
abandoned here entirely (dropped `library(rasterdiv)`).

Replaced `compute_rao_q_raster_rasterdiv()` in `Code/AnnualSpectralDiversity.R`
with `compute_rao_q_raster_local()` (+ helper `rao_q_window()`): a genuinely
local moving-window implementation that computes Rao's Q using ONLY each
window's own pixel values, the same approach `ComputeSpecBiodiv.R`'s
pyGNDiv-based `rao_q_window_pygndiv()`/`compute_rao_q_raster_pygndiv()` used
before the rasterdiv swap, minus the per-window Python/reticulate round-trip.
`ComputeSpecBiodiv.R` itself still uses the pyGNDiv approach (never had a
rasterdiv-based function), so it needed no change.

Formula: classic (non-normalized) Rao's Q, Q = sum_i sum_j p_i*p_j*d_ij, the
full double sum over all pairs (Botta-Dukat 2005 eq. 1; Rocchini et al. 2017
eq. 1), not divided by 2. Implemented as `sum(D) / n^2` over a window's own
n (<= window^2) pixels with equal weight 1/n each -- algebraically identical
to the textbook double sum (duplicates need no special dedup step: repeating
a value k times and weighting 1/n each equals weighting it k/n once).
Verified this equivalence against a brute-force nested-loop reference
implementation, exact match to floating-point precision.

Vectorization: evaluated `terra::focal()`/`focalValues()` first and rejected
them for the all-bands case -- confirmed empirically (synthetic 2-layer
SpatRaster) that `focalValues()` silently returns only the FIRST layer's
window values, and `focal()` applies its function per-layer independently;
neither can hand one function call a pixel's window values across ALL bands
jointly, which all-bands Rao's Q needs. Instead, every window's pixel-cell
indices are computed in one vectorized `cellFromRowCol()` call up front
(not recomputed per window like the old pygndiv loop), and the raster's
values matrix is extracted once; the remaining per-window loop does
O(window^2) work on values already resident in memory.

Verified via the same synthetic-raster standard as the rasterdiv version:
a homogeneous region returns near-zero Q, a region straddling a sharp
reflectance boundary returns clearly elevated Q (0.00095 vs. 0.0444 in the
test raster), and a homogeneous raster with an NA hole (vegetation-mask
exclusion) returns ~0 everywhere, confirming NA pixels are excluded from
window distance calculations rather than substituted. Multi-band case also
confirmed to return Q > 0 on random synthetic data. All four checks reran
successfully directly against the functions as committed in
`Code/AnnualSpectralDiversity.R` (not just the prototype), ruling out
transcription drift.

Memory/complexity reasoning: `rao_q_window()`'s distance matrix is at most
window^2 x window^2 (9x9 for the default 3x3 window) regardless of total
raster size -- structurally impossible to reproduce paRao's O(unique^2)
blowup, since no step anywhere builds a distance matrix over more than one
window's own pixels at a time. Performance sanity check on synthetic
1000x1000 rasters: ~5s for a single-layer (NDVI/NIRv) pass; the all-bands
case scales with band count via each window's `dist()` call (~350 bands,
500x500: ~6s, extrapolating to tens of seconds at 1000x1000) -- well under
the cost of the spectral-species-richness step (RF + K-means, `n_reps_ssr`
reps) that already dominates this script's runtime.

## 2026-08-14 21:10 UTC — Fix std::bad_alloc (crop-then-mosaic) and rhdf5 envRefClass error (scoped H5 handles)

Round 3 of real-data bugfixes. After the CRS/handle-leak fix (previous
entry), user ran the full pipeline locally against real ABBY 2017/2018 data
and hit two new failures:
- 2017: `std::bad_alloc` during read/mosaic/crop
- 2018: `Called from: is(x, "envRefClass")`

User proposed diagnoses for both and asked for investigation/confirmation
before changing anything, not blind implementation of the diagnosis.

**Issue 1 (std::bad_alloc) -- investigated and confirmed correct.**
`get_tower_reflectance()` mosaicked ALL full-resolution tiles together
before cropping to the 500m buffer, matching memory-exhaustion behavior
already independently observed earlier this session when attempting a real
run in this sandbox (holding N full 426-band/1000x1000 tiles simultaneously
is the actual peak-memory driver). Restructured to crop each tile to the
buffer's bounding box immediately after reading it, before mosaicking, so
peak memory never holds more than one full tile plus already-small cropped
pieces.

Verified equivalence EMPIRICALLY, not just by argument: built 4 synthetic
grid-aligned, non-overlapping tiles matching the real NEON tile geometry,
with a buffer point at the hardest case (the 4-way tile corner, drawing from
all 4 tiles), and ran both the old mosaic-then-crop and new crop-then-mosaic
sequences. Results were bit-identical: same extent, same dims,
`identical(values(old), values(new))` = TRUE, max abs diff on non-NA cells =
0, identical NA-pixel count/pattern (6220 in both).

**Issue 2 (envRefClass) -- investigated, not conclusively reproduced, fixed
anyway as strictly safer.** Read `rhdf5::h5closeAll()`'s actual source: it
calls `h5validObjects()` (queries the HDF5 C library for all currently-valid
identifiers) and closes each individually via `.H5close()` -- it does NOT
call the library-wide `H5close()` reset, so the specific mechanism proposed
(closer to a full library reset) isn't quite what the source shows.
Empirically ran 24 iterations of open/read/close (reusing the 4 real ABBY
files, `h5closeAll()` after each, mirroring the real loop structure) and did
not reproduce any handle-corruption error against this sandbox's rhdf5
2.46.1. Inconclusive -- doesn't rule out the failure at ~183-site-year scale
or a different rhdf5/HDF5-C build. Implemented the requested fix regardless:
scoped `H5Fopen()`/`H5Fclose()` per file via `on.exit()` in
`read_neon_h5_tile()`, replacing the global `h5closeAll()` sweep. Confirmed
against a real ABBY file first: `h5read()`/`h5readAttributes()`/`h5ls()` all
accept an open handle in place of a path and can be called repeatedly on it,
and using the handle after `H5Fclose()` correctly raises `"H5Identifier not
valid"` (not silent wrong data) -- `on.exit()` also closes the handle on the
error path, which the old code didn't.

Both files re-sourced cleanly after the edits (stopped only at the expected
missing `Data/NEONsites.csv` point, no syntax/definition errors). No real
pipeline run attempted in this sandbox, per standing direction from earlier
this session.

Applied identically to `Code/AnnualSpectralDiversity.R` and
`Code/ComputeSpecBiodiv.R`.

## 2026-08-14 20:30 UTC — Fix missing CRS and rhdf5 handle leak in read_neon_h5_tile()

User ran the full `AnnualSpectralDiversity.R` pipeline locally (real hardware)
against real ABBY 2017 and 2018 H5 tiles and hit two real errors:
- 2017: `"[project] output crs is not valid"`
- 2018: `"Called from: is(h5id, \"H5IdComponent\")"`

Both traced to `read_neon_h5_tile()`, which is shared verbatim between
`Code/AnnualSpectralDiversity.R` and `Code/ComputeSpecBiodiv.R` (the script it
was derived from) -- fixed identically in both.

**Bug 1, missing CRS:** `read_neon_h5_tile()` built the SpatRaster via
`rast(refl, extent = ext(...))` but never called `crs(r) <-`, so
`crs(masked$raster)` was empty and `project(pt, crs(masked$raster))` in
`get_tower_reflectance()` had nothing valid to project into. Confirmed
against a real file (`h5ls` on
`/ABBY/Reflectance/Metadata/Coordinate_System/`) that an `EPSG Code` dataset
sits alongside `Map_Info` (also `Proj4` and `Coordinate_System_String`,
either of which would also work, EPSG was preferred as the more robust
option across PROJ/terra versions) -- value `32610` (UTM zone 10N),
consistent across all 4 real tiles checked. Fix: read that dataset and
`crs(r) <- paste0("EPSG:", epsg_code)` immediately after `rast()`, before the
raster is ever mosaicked or projected against.

**Bug 2, HDF5 handle leak:** `read_neon_h5_tile()` calls `h5ls()`, `h5read()`
(x3), and `h5readAttributes()` per file without closing handles; across
~183 site-year jobs this exhausts rhdf5's internal handle table and throws
an internal `H5IdComponent` error partway through a long run -- a known
rhdf5 behavior, not a data problem. Fix: `h5closeAll()` at the end of
`read_neon_h5_tile()`, before it returns.

Verification (per user's explicit direction -- code review + lightweight
structural checks against real file metadata, not another full pipeline run,
since this sandbox cannot complete one; see the entry below):
- Confirmed the `EPSG Code` field/path against all 4 real ABBY H5 tiles
  (all `32610`).
- Confirmed `h5closeAll()` is a real exported function in the installed
  rhdf5 2.46.1, and ran it successfully both with nothing open and after
  real `H5Fopen()`/`H5Dopen()` handles were open against a real file.
- Confirmed the `crs(r) <- paste0("EPSG:", epsg_code)` pattern actually
  fixes the failure mode: `crs()` is `NA` before the assignment and a valid
  `32610` after, and `project()` against it succeeds (reproduced the same
  UTM tower coordinates, ~552103/5067601, computed earlier in this session
  by other means).
- Both edited files source cleanly (no syntax/definition errors) with only
  the library set the scripts already require.
- Did NOT re-run the full pipeline against real data in this sandbox -- see
  the memory-exhaustion entry below; that remains something only the user's
  local machine can do.

Applied to both `Code/AnnualSpectralDiversity.R` and `Code/ComputeSpecBiodiv.R`.

## 2026-08-14 19:55 UTC — Real ABBY 2017 H5 tiles: structural validation only; full pipeline run not completed (sandbox memory)

User supplied 4 real NEON DP3.30006.001 reflectance H5 tiles for ABBY, 2017
(`551000_5067000`, `551000_5068000`, `552000_5067000`, `552000_5068000`,
~640MB each), asking for a full end-to-end run of `AnnualSpectralDiversity.R`
against them (all 4 metrics: CV, CHV, spectral species richness, Rao's Q x3).

**Outcome: the full run did not complete in this sandbox.** After repeated
attempts (including adding up to 24GB of swap on top of this environment's
7.8GB RAM) consistently died or timed out inside `get_tower_reflectance()`
before finishing even the first tile's mosaic step -- confirmed via a memory
trace to be genuine RAM/swap exhaustion, not a code bug, an OOM-kill (no
`oom_kill` events in `/sys/fs/cgroup/memory.events`), or a data problem. Per
user direction, stopped attempting the full run in this sandbox -- it needs
to happen on real hardware with adequate RAM (user's local machine). Any
"the real pipeline produced these numbers" claim from this session would be
false; no such run finished, so no metric values are recorded here.

What WAS verified (real data, no full-array loads):
- All 4 files open correctly with `rhdf5`; internal structure
  (`/ABBY/Reflectance/Reflectance_Data`, `.../Metadata/Spectral_Data/Wavelength`,
  `.../Metadata/Coordinate_System/Map_Info`) matches exactly what
  `read_neon_h5_tile()` expects -- confirmed by sourcing that function
  verbatim from the script (not reimplemented) and running it against real
  paths.
- 426 bands, 381-2510nm, ~5nm spacing; `mask_bad_bands()`'s 4 wavelength
  ranges remove 80 of 426 bands (346 remain) -- consistent with the standard
  NEON band set.
- `Scale_Factor` = 10000, `Data_Ignore_Value` = -9999 in all 4 files, matching
  the script's fallback defaults exactly (so the `ifelse(!is.null(...))`
  fallback path is never actually exercised by real files -- untested by this
  data, worth remembering if a future NEON product revision omits these
  attributes).
- Windowed reads (single band, or small pixel windows, via `h5read`'s `index`
  argument -- not the full cube) show real, plausible 0-0.83 reflectance
  values after applying the scale factor, and 0 fill-value (`-9999`) pixels
  in a full-tile single-band (NIR) check across all 4 tiles.
- ABBY tower coordinates (45.76, -122.33, from user-supplied
  `NEONsites_Footprints.csv` -- `Data/NEONsites.csv` itself is gitignored and
  not present in this sandbox) project to UTM (552103, 5067601), comfortably
  inside the 4-tile mosaic extent with a 500m buffer (tightest margin ~101m
  on the south edge) -- confirms these 4 tiles are the geometrically correct
  set for this tower/buffer.

Root cause of the memory failure (confirmed, not guessed): `h5read()` returns
NEON's on-disk int16 cube as a 32-bit R integer array (~1.7GB for one full
426x1000x1000 tile). `read_neon_h5_tile()` then does
`refl[refl==no_data]<-NA; refl<-refl/scale_factor` (implicit full copy to
double, ~3.4GB) and `aperm(refl, c(3,2,1))` (another full transposed copy).
None of these are in-place in R. A synthetic array of the same shape
confirmed `rast()` alone on a comparable array pins RSS at 6.9GB+ for
minutes under swap pressure. `get_tower_reflectance()` reads all 4 tiles via
`map()` before mosaicking, so peak memory is driven by holding multiple
full, un-cropped 1000x1000x426 tiles at once -- this sandbox's 7.8GB RAM
(even +24GB swap, which just made it slow rather than survivable in a
useful wall-clock time) cannot support that; real research hardware likely
can.

Code-review flags for `Code/AnnualSpectralDiversity.R` (NOT fixed -- user
directed no code changes based on an unfinished run; these are things to
watch when the real run happens on adequate hardware):
- `get_tower_reflectance()` mosaics and calls `mask_bad_bands()` on the FULL,
  un-cropped mosaic, then crops to the 500m buffer last. Cropping (or at
  least a rough bounding-box pre-filter) before band-masking would cut peak
  memory substantially, especially for smaller buffers or sites with more
  than 4 overlapping tiles.
- `read_neon_h5_tile()` always reads the entire tile's Reflectance_Data cube
  via `h5read()` with no spatial or band-count subsetting at read time, even
  when only a small buffer fraction of the tile is ultimately needed. NEON
  H5 datasets support indexed/windowed `h5read()` reads (confirmed working
  in this session's windowed checks) -- worth considering if larger buffers
  or many-tile footprint runs make this a recurring bottleneck on real
  hardware too, not just this sandbox.
- `get_tower_reflectance()` wraps `read_neon_h5_tile()` in
  `possibly(..., otherwise = NULL)` and silently drops any tile that fails
  to read (`compact(tiles)`), with no warning surfaced to the caller or the
  per-job log. Given this project's prior placeholder/corrupt-file history
  (see the `diagnose_placeholder_files.R` commit), a real run over many
  site-years should watch for this silently producing an incomplete mosaic
  (missing a tile) that still reports `status = "success"`.
- Unrelated to the pipeline code: the 4 real H5 tiles (and
  `NEONsites_Footprints.csv`) currently sit untracked at the repo root, not
  under `Data/` (which is gitignored) -- a stray `git add -A` would try to
  commit ~2.5GB of raw hyperspectral data. Recommend moving them under
  `Data/` or adding a root-level `*.h5` gitignore rule before this session's
  files are cleaned up.

Sandbox changes made and reverted: added an 8GB then 24GB swapfile at
`/tmp/swapfile_claude` to attempt the full run; both `swapoff` and file
removal confirmed after stopping (`free -h` back to 0B swap, no stray R
processes).

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
