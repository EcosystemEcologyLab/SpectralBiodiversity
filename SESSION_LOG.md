# Session Log

Committed, dated record of work performed in this repository. Reverse
chronological order, newest entry at the top. See CLAUDE.md for the
convention this file follows.

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
