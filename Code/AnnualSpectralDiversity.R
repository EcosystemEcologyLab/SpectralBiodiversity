# AnnualSpectralDiversity.R
#
# Interannual companion to ComputeSpecBiodiv.R. Separate research question
# from the flux-footprint pipeline (NEON_ExtractFootprints.R /
# AnnualFootprintExtent.R) and fully independent from it: different script,
# different output path, no shared state. This does NOT depend on footprint
# size -- it is not about matching a flux measurement to its source area.
#
# Question: does spectral diversity vary meaningfully year-to-year at a
# site, or can it be treated as a roughly constant per-site value? To answer
# that, each metric is computed separately for EVERY YEAR of hyperspectral
# data available at each site (one row per tower-YEAR), instead of one row
# per tower for the most recent year only.
#
# Metrics (identical definitions to ComputeSpecBiodiv.R, except Rao's Q --
# see note below):
#   1. Coefficient of Variation (CV)
#   2. Convex Hull Volume (CHV), first 3 PCs
#   3. Spectral Species Richness -- self-adaptive fuzzy c-means (FCM), Wu et
#      al. 2026 (Computers and Electronics in Agriculture 252:112108). The
#      previous RF-proximity + K-means (50 clusters) + nearest-centroid
#      approach (Feret & Asner 2014) was REMOVED entirely, not merely
#      replaced in this loop: weekend diagnostics on real ABBY/SRER data
#      (Code/KMeans_test.R, all k in {10..50} swept) found richness landing
#      at EXACTLY k every time (100% of k) for every site-year, with no
#      ABBY/SRER separation in Shannon's H' either -- confirming the
#      saturation problem flagged as an open problem in a prior version of
#      this file (see git history / SESSION_LOG.md for that post-mortem) was
#      never resolved by the nearest-centroid fix, and that the RF+K-means
#      approach carries no usable diversity signal on real data, at real
#      expense (~1500-1700 sec/site-year at k=50). Adaptive-FCM SSR showed
#      real, stable separation across 3 years in the same diagnostics
#      (Code/DataAnalysis/FCM_test.R, CompareSSR_AdaptiveFCM_vs_KMeans.R,
#      adaptive_fcm_ssr.py) at comparable-or-lower runtime, and is adopted
#      here as the production method. See Section 5 for the implementation
#      and its PROVISIONAL parameters (Section 0's c_max_ssr/fcm_fuzziness_w).
#      A prior version of this section's SSR used a c=1 "floor-check" branch
#      to try to catch an apparent c=2 floor artifact seen on real SRER
#      data (every rep's L(c) maximized at c=2 with no interior peak, unlike
#      ABBY's genuine interior peak at c~6-8); that branch was found dead by
#      its own self-test (fuzzy within-cluster dispersion decreases
#      monotonically in c even for pure noise, so c=1 essentially never won)
#      and has been REMOVED, replaced by a conditional gap-statistic check
#      (Tibshirani, Walther & Hastie 2001) that only runs for a rep whose
#      primary search lands exactly at the c=2 floor -- see Section 5 for
#      the full rule and its own provisional-parameter caveats.
#   4. Convex Hull Area (CHA), band-selection-based (Gholizadeh et al. 2018,
#      Remote Sensing of Environment 206:240-253) -- found to outperform CV,
#      CHV, and SID specifically at ~1m airborne resolution (this script's
#      resolution regime; CV degraded badly at this scale in that paper).
#      Distinct from CHV: CHV takes pixels as points in PC space and hulls
#      the pixel cloud; CHA takes BANDS as points in a 2D (mean-spectrum,
#      pixel-spectrum) space for one pixel at a time, and is averaged across
#      all vegetation pixels in the plot. No PCA, clustering, or RF -- see
#      Section 4a.
#   5. Rao's Q (NDVI, NIRv, all-bands), classic (non-normalized) form, via a
#      local per-window implementation (Section 6) -- NOT via
#      pyGNDiv/reticulate, and NOT via rasterdiv::paRao() (tried first; see
#      Section 6 for why it doesn't work at this resolution). Comparability
#      across variable-set sizes (e.g. a 1-band NDVI Q vs. an all-bands Q,
#      which is what pyGNDiv's generalizable normalization was for) is not
#      needed here: every comparison in this analysis is across YEARS at the
#      same site, using the same variable set each time. A hypothetical
#      future analysis that compared Rao's Q *across index types* directly
#      would need that normalization back; this one doesn't.
#      The all-bands variant is additionally PCA-reduced before Rao's Q is
#      computed (Section 6a) -- raoq_allbands measures Rao's Q on the
#      dominant ~N axes of spectral variation (99% variance, N data-driven
#      per site-year), not literally every band, consistent with how
#      CHV/SSR already reduce dimensionality via PCA. Earlier full-band
#      Rao's Q numbers (e.g. earlier ABBY test runs, or the rasterdiv-era
#      output) are NOT numerically comparable to this version's
#      raoq_allbands.
#   6. Spectral Shannon's index (H') and its Hill-number effective diversity
#      (exp(H'), Wang et al. 2018, Remote Sensing of Environment 211:218-228,
#      Table 2), computed from the same spectral-species cluster assignments
#      spectral species richness already produces (Section 5) -- no separate
#      classification step.
#
# NOTE on PCA usage (open decision, not acted on): CHV, SSR, and the
# all-bands Rao's Q still reduce dimensionality via PCA, unchanged by this
# session's edits. Whether to replace PCA with a fixed manual band set
# instead is an open decision, deliberately deferred so the nearest-centroid
# SSR fix and the new CHA metric (which uses no PCA at all) can be tested in
# isolation first.
#
# Differences from ComputeSpecBiodiv.R:
#   - a single fixed 500m buffer is used for ALL METRICS. The original
#     script's two different buffers (530m for CV/CHV/SSR, 300m for Rao Q)
#     were tied to footprint-related design choices that don't apply to this
#     question -- every metric for every site-year uses the same 500m buffer
#     here, and the reflectance raster is read once per site-year and reused
#     for all four metrics (rather than reading it twice at two buffer sizes).
#   - one row per tower-YEAR: loops over whatever years of hyperspectral
#     tiles are actually present on disk for each site under
#     NEON_Hyperspec/<tower_id>/, the same way AnnualFootprintExtent.R
#     discovers available years per site rather than assuming a fixed list.
#   - results are written to CSV incrementally, after each site-year, so an
#     interrupted run doesn't lose completed work (matches
#     AnnualFootprintExtent.R).
#   - CHV is still z-score standardized, but now across all site-YEAR
#     observations (not across sites only), since the point is to compare a
#     site's own values across years as well as against other sites.
#
# All H5-reading, spectral-index, and metric-computation helpers below are
# reused verbatim from ComputeSpecBiodiv.R -- this script only restructures
# the loop and buffer logic around them -- EXCEPT Rao's Q (Section 6), which
# uses a local per-window R implementation instead of ComputeSpecBiodiv.R's
# pyGNDiv approach, and EXCEPT Spectral Species Richness (Section 5), which
# uses adaptive-FCM (reticulate + adaptive_fcm_ssr.py) instead of
# ComputeSpecBiodiv.R's still-unchanged RF+K-means approach; see the notes
# above and in Sections 5/6 for why.
#
# Runtime: dominated by spectral species richness (adaptive-FCM, n_reps_ssr
# reps per site-year, ~19 candidate-c fits per rep for the primary search's
# c in [2, c_max_ssr], PLUS a conditional gap-statistic fallback -- up to
# gap_check_c_max*(1 + gap_check_n_null) additional fits -- on whichever
# reps land at the c=2 floor; see Section 5's header for the trigger rule)
# and the three local Rao's Q moving-window passes (one per variant, no
# per-window Python round-trip and no raster-wide distance matrix).
# Diagnostics found adaptive-FCM SSR runtime comparable to or cheaper than
# the single-k-means-fit approach it replaced when the gap fallback isn't
# triggered; see Section 5's header for the diagnostic history and the
# per-site-year gap-fallback timing note the main loop now prints. CHA
# (Section 4a) adds a
# fast per-pixel geometric pass with no PCA/clustering/RF, negligible next to
# SSR and Rao's Q. The 500m
# buffer here is smaller than ComputeSpecBiodiv.R's CV/CHV/SSR buffer (530m,
# so that part of the cost is about the same) but larger than its Rao Q
# buffer (300m) -- area scales as radius^2, so each Rao Q variant here covers
# roughly (500/300)^2 ~= 2.8x as many pixels/windows as the original script's
# Rao Q step. Total run time scales with the number of site-YEARS, not the
# number of sites: a site with 4 years on disk costs ~4x a single-year site.
# Rather than guess an absolute number, the loop below times and prints each
# completed site-year and a running average/ETA, so you can extrapolate for
# your actual job count.

Sys.setenv(PROJ_LIB = "C:/Program Files/R/R-4.4.1/library/terra/proj")

library(terra)
library(rhdf5)          # low-level NEON H5 reading
library(geometry)        # convhulln for CHV
library(reticulate)      # adaptive-FCM SSR bridge (Section 5), adaptive_fcm_ssr.py
library(cluster)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# ============================================================================
# 0. Setup
# ============================================================================
hyperspec_dir <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
out_csv       <- "D:/projects/moore/SpectralBiodiversity/Data/spectral_diversity_by_year.csv"

buffer_m           <- 500   # single buffer for ALL metrics (CV, CHV, SSR, RaoQ x3)
ndvi_thresh         <- 0.4
n_subsample_pixels   <- 2500
n_pc_ssr             <- 4
n_pc_chv             <- 3
n_reps_ssr           <- 20
raoq_window          <- 3   # side of the square moving window for local Rao's Q (Section 6)
raoq_pca_var_threshold <- 0.99   # variance fraction to retain when PCA-reducing all-bands Rao's Q (Section 6a)

# ---- Adaptive-FCM spectral species richness config (Section 5) ----
# ALL VALUES BELOW ARE PROVISIONAL -- see Section 5's header. Kept as named
# constants (not buried inline) so a follow-up session can change them
# without hunting through the function body.
c_max_ssr       <- 20    # unchanged from diagnostics; primary search always
                          # runs c in [2, c_max_ssr] (see select_winning_c()).
fcm_fuzziness_w <- 2.0   # PROVISIONAL -- standard/default FCM fuzziness
                          # exponent (Wu et al.'s w), not tuned for this data.
gap_check_c_max <- 8     # PROVISIONAL -- c range [1, gap_check_c_max] tested
                          # by the gap-statistic fallback (Section 5) when a
                          # rep's primary search lands at the c=2 floor. Only
                          # meant to adjudicate c=1 vs. small c, not to re-run
                          # the full [2, c_max_ssr] search -- see Section 5.
gap_check_n_null <- 10   # null (uniform-random) references per candidate c
                          # in the gap-statistic fallback.
adaptive_fcm_py <- "./Code/DataAnalysis/adaptive_fcm_ssr.py"

# ============================================================================
# 0a. Python/FCM environment check -- printed unconditionally at the VERY TOP
#     of execution, before towers_df or anything else is read (not only on
#     failure), so an unattended weekend run's log clearly shows which
#     Python/reticulate environment was actually used, for later debugging
#     if results look off. Fails HERE, before any site-year work starts,
#     rather than deep into a multi-hour loop.
# ============================================================================
if (!file.exists(adaptive_fcm_py)) {
  stop("Required input not found: ", adaptive_fcm_py)
}

cat("==== reticulate Python configuration ====\n")
print(py_config())
cat("==========================================\n\n")

if (!py_module_available("numpy")) {
  stop("Python module 'numpy' is not importable via reticulate. Check ",
       "RETICULATE_PYTHON / reticulate::py_config() above before running ",
       "the real loop -- do not proceed on a guess.")
}
if (!py_module_available("fcmeans")) {
  stop("Python module 'fcmeans' (PyPI package 'fuzzy-c-means') is not ",
       "importable via reticulate. Install it in whatever Python ",
       "reticulate is pointed at (e.g. `pip install fuzzy-c-means`), or ",
       "point RETICULATE_PYTHON at an environment that already has it, ",
       "before running the real loop -- do not proceed on a guess.")
}
source_python(adaptive_fcm_py)   # exposes adaptive_fcm_ssr(), fcm_validity_L(), FCM
                                   # (the fcmeans.FCM class) -- NOT _fit_one(), which
                                   # source_python() does not bind (leading underscore);
                                   # see select_winning_c()'s header (Section 5).
cat("  numpy: available\n  fcmeans (fuzzy-c-means): available\n\n")

towers_df <- read.csv("./Data/NEONsites.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites in NEONsites.csv.\n")

# ============================================================================
# 1. NEON reflectance reading helpers (reused from ComputeSpecBiodiv.R)
# ============================================================================
# NEON DP3.30006.00x H5 tiles store reflectance as int16 scaled by a factor
# (typically 10000) with a fill/no-data value (typically -9999 or -100).
# Wavelengths (nm) are stored as a dataset alongside the reflectance cube.
# Bad bands are the standard water-vapor absorption + edge regions widely
# used in NEON hyperspectral tutorials -- verify against your product
# version's metadata before trusting these ranges blindly.
bad_band_ranges <- list(c(300, 400), c(1340, 1445), c(1790, 1955), c(2400, 2600))

read_neon_h5_tile <- function(h5_path) {
  # Open this file's handle once and pass it into every h5*() call below,
  # instead of the path (which rhdf5 would open+close internally per call).
  # Scoped to this one file/call via on.exit -- never touches any other
  # file's handles or session-wide HDF5 state. (An earlier version used
  # h5closeAll() here, which sweeps every currently-open identifier in the
  # whole R session; that was suspected of leaving rhdf5 in a state that
  # broke a later, unrelated file's read. Not conclusively reproduced in
  # testing, but scoped per-file open/close is strictly more precise
  # regardless and is rhdf5's documented pattern for read loops.)
  fid <- H5Fopen(h5_path)
  on.exit(H5Fclose(fid), add = TRUE)

  refl_path <- h5ls(fid) %>% filter(str_detect(name, "Reflectance$")) %>% pull(group) %>% unique()
  site_group <- str_split(refl_path, "/")[[1]][2]

  wavelengths <- h5read(fid, paste0("/", site_group, "/Reflectance/Metadata/Spectral_Data/Wavelength"))
  refl <- h5read(fid, paste0("/", site_group, "/Reflectance/Reflectance_Data"))
  attrs <- h5readAttributes(fid, paste0("/", site_group, "/Reflectance/Reflectance_Data"))
  scale_factor <- ifelse(!is.null(attrs$Scale_Factor), attrs$Scale_Factor, 10000)
  no_data <- ifelse(!is.null(attrs$Data_Ignore_Value), attrs$Data_Ignore_Value, -9999)

  map_info <- h5read(fid, paste0("/", site_group, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
  map_info <- str_split(map_info, ",")[[1]]
  px_size  <- as.numeric(map_info[6])
  x_min    <- as.numeric(map_info[4])
  y_max    <- as.numeric(map_info[5])
  epsg_code <- h5read(fid, paste0("/", site_group, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))

  refl[refl == no_data] <- NA
  refl <- refl / scale_factor   # -> 0-1 reflectance

  # H5 array dims are (band, col, row); reorder to (row, col, band) for terra
  refl <- aperm(refl, c(3, 2, 1))
  nrow_r <- dim(refl)[1]; ncol_r <- dim(refl)[2]

  r <- rast(refl, extent = ext(x_min, x_min + ncol_r * px_size,
                               y_max - nrow_r * px_size, y_max))
  crs(r) <- paste0("EPSG:", epsg_code)
  names(r) <- paste0("b", round(wavelengths))

  list(raster = r, wavelengths = wavelengths)
}

# drop water-absorption / edge bands
mask_bad_bands <- function(r, wavelengths) {
  bad <- rep(FALSE, length(wavelengths))
  for (rng in bad_band_ranges) bad <- bad | (wavelengths >= rng[1] & wavelengths <= rng[2])
  list(raster = r[[!bad]], wavelengths = wavelengths[!bad])
}

# Read + mosaic a GIVEN set of h5 tiles (already filtered to one tower-YEAR),
# crop to a single circular buffer around the tower. Unlike
# ComputeSpecBiodiv.R's get_tower_reflectance() -- which globs every h5 file
# under the tower directory and reads it at two different buffer sizes --
# this takes explicit file paths (one tower-year's worth) and one buffer, and
# is read once per site-year for reuse across all four metrics.
get_tower_reflectance <- function(tower_id, h5_files, buffer_m) {
  if (length(h5_files) == 0) return(NULL)

  tower_row <- towers_df %>% filter(Site.ID == tower_id)
  pt <- vect(cbind(tower_row$Lon, tower_row$Lat), crs = "EPSG:4326")

  # Crop each tile to the tower's buffer immediately after reading it, BEFORE
  # mosaicking, rather than mosaicking all full-resolution tiles together and
  # cropping once at the end. Holding N full, un-cropped tiles in memory
  # simultaneously (mosaic-then-crop) is what exhausted memory (std::bad_alloc)
  # on real multi-tile site-years; crop-then-mosaic never holds more than one
  # full tile plus a set of already-small, buffer-sized pieces.
  #
  # This does not change the result. mask_bad_bands() only selects LAYERS
  # (by wavelength), never touches spatial extent, so it commutes with
  # crop/mosaic order regardless of where it runs. And NEON tiles are
  # grid-aligned and non-overlapping (fixed 1km tiles on a common 1m-pixel
  # UTM grid), so mosaic() assigns each output cell the value from whichever
  # single tile covers it -- there is no blending to get wrong. Cropping a
  # tile to the buffer's bounding box before mosaicking only discards cells
  # that the old mosaic-then-crop order would also discard in its final
  # crop(); it cannot change the value of any cell that survives either way.
  crop_tile_to_buffer <- function(h5_path) {
    tile <- read_neon_h5_tile(h5_path)
    pt_proj <- project(pt, crs(tile$raster))
    buf <- buffer(pt_proj, buffer_m)
    tile$raster <- crop(tile$raster, buf)
    tile
  }

  tiles <- map(h5_files, possibly(crop_tile_to_buffer, otherwise = NULL))
  tiles <- compact(tiles)
  if (length(tiles) == 0) return(NULL)

  wavelengths <- tiles[[1]]$wavelengths
  rasters <- map(tiles, "raster")
  mosaic_r <- if (length(rasters) > 1) do.call(terra::mosaic, rasters) else rasters[[1]]

  masked <- mask_bad_bands(mosaic_r, wavelengths)

  # final circular mask + crop -- the pieces are already buffer-bbox-sized,
  # so this crop is now cheap/near-no-op; the mask() to the circle (rather
  # than just the bbox) still has to happen once on the assembled mosaic.
  pt_proj <- project(pt, crs(masked$raster))
  buf <- buffer(pt_proj, buffer_m)
  cropped <- crop(mask(masked$raster, buf), buf)
  list(raster = cropped, wavelengths = masked$wavelengths)
}

# ============================================================================
# 2. Spectral index helpers (reused from ComputeSpecBiodiv.R)
# ============================================================================
nearest_band <- function(wavelengths, target_nm) which.min(abs(wavelengths - target_nm))

compute_ndvi_raster <- function(r, wavelengths) {
  red_i <- nearest_band(wavelengths, 670)
  nir_i <- nearest_band(wavelengths, 800)
  (r[[nir_i]] - r[[red_i]]) / (r[[nir_i]] + r[[red_i]])
}

compute_nirv_raster <- function(r, wavelengths) {
  ndvi <- compute_ndvi_raster(r, wavelengths)
  nir_i <- nearest_band(wavelengths, 800)
  ndvi * r[[nir_i]]
}

# ============================================================================
# 3. Metric 1: Coefficient of Variation (reused from ComputeSpecBiodiv.R)
# ============================================================================
compute_cv <- function(r, veg_mask) {
  r_masked <- mask(r, veg_mask)
  band_cv <- global(r_masked, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  mean(band_cv[, 1], na.rm = TRUE)
}

# ============================================================================
# 4. Metric 2: Convex Hull Volume (reused from ComputeSpecBiodiv.R)
# ============================================================================
compute_chv <- function(r, veg_mask, n_pc = 3) {
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  if (nrow(vals) < n_pc + 1) return(NA_real_)
  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pcs <- pca$x[, 1:n_pc]
  ch <- tryCatch(convhulln(pcs, output.options = "FA"), error = function(e) NULL)
  if (is.null(ch)) return(NA_real_)
  ch$vol
}

# ============================================================================
# 4a. Metric: Convex Hull Area (CHA), band-selection-based (Gholizadeh et al.
#     2018, Remote Sensing of Environment 206:240-253) -- found to outperform
#     CV, CHV, and SID specifically at ~1m airborne resolution (this script's
#     resolution regime; CV degraded badly at this scale in that paper).
# ============================================================================
# Distinct from CHV above: CHV treats PIXELS as points in PC space and hulls
# the pixel cloud (one hull per plot). CHA treats BANDS as points in a 2D
# space -- for ONE pixel, each band contributes one point (plot mean
# reflectance at that band, this pixel's reflectance at that band) -- and is
# averaged across all vegetation pixels in the plot. No PCA, clustering, or
# RF -- a fast per-pixel geometric calculation using the same
# geometry::convhulln(..., output.options = "FA") approach as CHV; for a 2D
# point set $vol is the qhull "volume" field, which for 2 dimensions is the
# polygon's area (not $area, which is the perimeter -- confirmed against a
# hand-computed triangle area during validation).
#
# A pixel whose spectrum equals the plot mean spectrum exactly produces a
# perfectly collinear point set (every point falls on the y = x line), which
# qhull cannot build a hull from -- it errors ("Initial simplex is flat")
# rather than returning 0. That collinear-input error is mapped to area 0
# here (not NA): a degenerate, zero-width point set has zero enclosed area by
# definition, and this is an expected, not exceptional, case (a pixel
# identical to the plot mean).
compute_cha <- function(r, veg_mask) {
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)   # n_pixels x n_bands
  if (nrow(vals) == 0) return(NA_real_)
  mean_spectrum <- colMeans(vals, na.rm = TRUE)

  cha_vals <- apply(vals, 1, function(px) {
    pts <- cbind(mean_spectrum, px)
    ch <- tryCatch(convhulln(pts, output.options = "FA"), error = function(e) NULL)
    if (is.null(ch)) return(0)
    ch$vol
  })
  mean(cha_vals, na.rm = TRUE)
}

# ============================================================================
# 5. Metric 3: Spectral Species Richness -- self-adaptive fuzzy c-means
#    (FCM), Wu et al. 2026 (Computers and Electronics in Agriculture
#    252:112108), plus spectral Shannon's index H' and its Hill-number
#    effective diversity exp(H'), reusing the SAME per-rep FCM output -- NOT
#    a separate classification pass. See the file header for why the prior
#    RF-proximity + K-means + nearest-centroid approach was removed entirely
#    (it is NOT preserved below, not even as a fallback).
# ============================================================================
# Returns a named list: list(richness, shannon_h, shannon_effective,
# mean_winning_c, max_winning_c, n_reps_gap_check_triggered, primary_time_sec,
# gap_time_sec). Every call site MUST unpack all eight (see Section 8) --
# unpack explicitly, don't rely on partial-match coercion.
#
# ALL PARAMETERS THIS FUNCTION IS CALLED WITH (Section 0's c_max_ssr,
# fcm_fuzziness_w, gap_check_c_max, gap_check_n_null) ARE PROVISIONAL --
# flagged there, not repeated here.
#
# SHANNON'S H' POPULATION-SCOPE CAVEAT: unlike the removed k-means version's
# H' (computed over a FULL-population nearest-centroid reassignment), this
# H' is computed over the winning c's HARD assignment (argmax of fuzzy
# membership, or the single group for c=1) of that rep's SUBSAMPLE only --
# there is no full-population reassignment step for FCM here (see
# compute_spectral_species_richness_fcm() below). Same formula as before,
# NOT the same population scope -- a genuine change in what H' measures, not
# just in method, inherited unchanged from the FCM_test.R diagnostic this
# replaces.
#
# select_winning_c() (PRIMARY SEARCH, c in [2, c_max]): reuses
# adaptive_fcm_ssr.py's own, UNMODIFIED fcm_validity_L() (actual
# validity-function arithmetic) AND the same `fcmeans.FCM` class that
# module's own (private) _fit_one() calls internally -- calling
# `FCM(n_clusters=c, m=fuzziness, random_state=..., verbose=FALSE);
# model$fit(X); model$u; model$centers` directly here (via fit_fcm(), shared
# with the gap-statistic fallback below for the fit step only -- see that
# function's header) reproduces _fit_one()'s exact body (confirmed by
# reading adaptive_fcm_ssr.py's source), NOT a different or reimplemented
# fitting routine, and uses the identical FCM class object that module
# imports. _fit_one() itself could not be called directly here -- confirmed
# empirically that reticulate::source_python() does not bind
# underscore-prefixed Python names into R at all (only fcm_validity_L,
# adaptive_fcm_ssr, and FCM appear after sourcing the module; a leading
# underscore is Python's own "not public API" convention, and reticulate
# respects it), so calling it was never actually an option regardless of
# preference -- this is the closest available reuse of the module's real
# fitting code without editing the file.
#
# THE C=1 FLOOR-CHECK BRANCH THAT USED TO LIVE HERE HAS BEEN REMOVED, NOT
# PATCHED. It compared c=1's raw within-cluster dispersion against
# min(within-cluster dispersion across c=2..c_max) and returned c=1 iff that
# was smaller. Self-testing on synthetic homogeneous data (pure Gaussian
# noise, no real substructure) found this branch essentially unwinnable: it
# never won across 5 reps, because fuzzy within-cluster dispersion is
# MONOTONICALLY DECREASING in c even for pure noise (more fuzzy prototypes
# mechanically reduce total weighted squared-distance-to-nearest-center
# regardless of whether the extra prototypes fit real signal or just noise
# -- the same reason k-means SSE is monotonically non-increasing in k), so
# the multi-cluster side of that comparison was always anchored near
# c_max's value and c=1 essentially never beat it. That branch is GONE
# (previously: `if (c_min == 1) { within_1 <- ...; if (within_1 <=
# min_within_ge2) return(list(winning_c = 1L, ...)) }`), replaced below by
# the conditional gap-statistic rule.
#
# GAP-STATISTIC FALLBACK (Tibshirani, Walther & Hastie 2001, JRSS-B
# 63:411-423), conditional -- NOT run for every rep. The full gap-statistic
# check (null-reference generation + FCM fits across gap_check_n_null
# references x multiple candidate c) is too expensive to run unconditionally
# for every rep of every site-year -- it would meaningfully undercut the
# compute savings that motivated dropping k-means in the first place. So
# it's gated: it only runs for a rep whose PRIMARY search (select_winning_c,
# c in [2, c_max]) lands exactly at the c=2 floor (argmax_c L(c) == 2). For
# every other rep (winning c == 3, 4, ..., c_max), the primary search's
# result is used directly and the gap check never runs.
#
# When triggered, compute_gap_statistic() tests c in [1, gap_check_c_max]
# (a SMALLER, reduced range than the primary search's [2, c_max] -- the
# point is only to adjudicate c=1 vs. small c, not to re-run the full
# search) using gap_check_n_null uniform-random null references per
# candidate c, generated within the OBSERVED subsample's own per-PC-dimension
# min/max range (Tibshirani's simpler reference method, not a PCA-realigned
# bounding box). W_k for the gap statistic is the RAW POOLED WITHIN-CLUSTER
# SUM OF SQUARED DISTANCES UNDER HARD ASSIGNMENT (argmax of FCM's fuzzy
# membership matrix, or the single group for c=1) -- see fit_hard_wk()
# below. This is a DIFFERENT quantity from select_winning_c()'s fuzzy-
# weighted `within_norm`, computed differently (raw hard-assignment SS to a
# cluster's own mean, vs. u^fuzziness-weighted squared distance to the FCM
# prototype) -- the two are not interchangeable and the code below does not
# try to share anything between them beyond the FCM fit call itself
# (fit_fcm()).
#
# Optimal-c selection (pick_optimal_c(), below): finds the GLOBAL maximum of
# Gap(c) across the tested range first, then picks the SMALLEST c whose Gap
# is within one standard error (sk) of that global max -- NOT a naive
# left-to-right scan from c=1, which breaks on a non-monotonic Gap curve
# (confirmed by self-test below against a synthetic dip-before-peak Gap
# sequence).
#
# If the gap check picks c=1 for a rep: that rep's richness value is 1,
# entered into the site-year's richness average like any other rep's value
# (simple mean across all n_reps reps, unchanged). If it picks c>=2: that c
# is used (the gap statistic's answer, not the primary search's c=2, since
# it's the more rigorous check when it has actually been run -- the two
# usually agree in that case, but if they don't, the gap check's result
# wins).
fit_fcm <- function(X, c, fuzziness, seed) {
  model <- FCM(n_clusters = as.integer(c), m = fuzziness,
               random_state = as.integer(seed), verbose = FALSE)
  model$fit(X)
  model
}

select_winning_c <- function(pcs_sub, c_max, fuzziness, seed) {
  L_vals <- setNames(numeric(c_max - 1), as.character(2:c_max))
  u_by_c <- list()

  for (c in 2:c_max) {
    model <- fit_fcm(pcs_sub, c, fuzziness, seed + c)
    u <- model$u; centers <- model$centers
    L_vals[as.character(c)] <- fcm_validity_L(pcs_sub, u, centers, fuzziness)
    u_by_c[[as.character(c)]] <- u
  }

  finite_L <- L_vals[is.finite(L_vals)]
  if (length(finite_L) == 0) {
    stop("select_winning_c: L(c) was non-finite for every candidate c in ",
         "[2, ", c_max, "] -- investigate the input/fit, do not guess a winner.")
  }
  winning_c <- as.integer(names(which.max(finite_L)))

  list(winning_c = winning_c,
       hard = max.col(u_by_c[[as.character(winning_c)]], ties.method = "first"))
}

# Purely a log(0) guard for fit_hard_wk() below -- NOT a science/tuning
# parameter (unlike Section 0's PROVISIONAL c_max_ssr/fcm_fuzziness_w/
# gap_check_* constants), so kept separate from that block. W_k=0 happens
# iff every point in every hard-assigned cluster coincides EXACTLY with
# that cluster's own centroid -- a degenerate edge case (e.g. an exact
# duplicate point pair collapsing to a singleton "cluster", or a genuinely
# zero-variance PC subsample), vanishingly unlikely on continuous
# reflectance PC data but not impossible across an unattended, multi-day,
# 45-site x multiple-year x n_reps_ssr-rep run. Without this floor,
# log(W_k) -> -Inf and silently propagates into that candidate c's Gap(c),
# which is a worse failure mode than a loud one -- it wouldn't stop the run
# or necessarily be obvious in the output, but could corrupt that rep's
# Gap(c)/pick_optimal_c() result.
epsilon_wk <- 1e-10

# Raw pooled within-cluster sum of squared distances under HARD assignment,
# for one candidate c on one data matrix X -- the gap statistic's W_k
# (Tibshirani et al. 2001), computed via each hard cluster's OWN mean (not
# the FCM prototype), which for squared Euclidean distance is algebraically
# identical to that paper's sum-of-pairwise-distances/(2*n_r) formulation
# (the standard implementation convention, e.g. cluster::clusGap). c=1 is
# handled directly (no FCM fit needed): a single group, dispersion around
# the global mean. Returns list(wk, hard) -- `hard` is kept so the caller
# can reuse the OBSERVED side's cluster assignment for Shannon's H' without
# a second fit, if this c is the one that ultimately wins.
#
# epsilon_wk floor applied at ONE finalization point below (after either
# branch computes its raw wk, before either returns) -- not duplicated
# across the c==1 / c>=2 branches. This only guards the exact-W_k==0 edge
# case; it does not meaningfully perturb wk for any non-degenerate cluster
# (1e-10 vs. any real-data sum-of-squared-distances is negligible). Logs a
# cat() (this file's established progress-logging convention, not a hard
# failure) when it actually engages, so the edge case is visible in an
# unattended run's log rather than silently smoothed over.
fit_hard_wk <- function(X, c, fuzziness, seed) {
  n <- nrow(X)
  if (c == 1) {
    centroid <- colMeans(X)
    wk <- sum(rowSums(sweep(X, 2, centroid) ^ 2))
    hard <- rep(1L, n)
  } else {
    model <- fit_fcm(X, c, fuzziness, seed)
    hard <- max.col(model$u, ties.method = "first")
    wk <- 0
    for (k in seq_len(c)) {
      idx <- which(hard == k)
      if (length(idx) == 0) next
      pts <- X[idx, , drop = FALSE]
      centroid <- if (length(idx) == 1) pts[1, ] else colMeans(pts)
      wk <- wk + sum(rowSums(sweep(pts, 2, centroid) ^ 2))
    }
  }

  if (wk < epsilon_wk) {
    cat("    WARNING: fit_hard_wk() W_k hit the epsilon floor for c=", c,
        " -- cluster(s) collapsed to their centroid exactly (raw wk=", wk, ").\n", sep = "")
    wk <- epsilon_wk
  }

  list(wk = wk, hard = hard)
}

# Gap(c) = mean_b[log(W_k,null_b)] - log(W_k,observed), sk = sd(log(W_k,null))
# * sqrt(1 + 1/n_null), for c in c_range -- Tibshirani et al. 2001's simpler
# (non-PCA-aligned) uniform reference distribution, bounded to the OBSERVED
# subsample's own per-dimension min/max range. Returns list(gap_df, hard_by_c)
# -- hard_by_c holds the OBSERVED side's per-c hard assignment (indexed by
# fit_hard_wk() above), for the caller to pick up once pick_optimal_c() has
# decided the winning c, without a second fit.
compute_gap_statistic <- function(pcs_sub, c_range, fuzziness, seed, n_null = 10) {
  n_pixels <- nrow(pcs_sub); n_dims <- ncol(pcs_sub)
  dim_min <- apply(pcs_sub, 2, min)
  dim_range <- apply(pcs_sub, 2, max) - dim_min

  observed <- lapply(c_range, function(c) fit_hard_wk(pcs_sub, c, fuzziness, seed + c))
  names(observed) <- as.character(c_range)
  wk_observed <- vapply(observed, function(o) o$wk, numeric(1))

  log_wk_null <- matrix(NA_real_, nrow = n_null, ncol = length(c_range))
  for (b in seq_len(n_null)) {
    ref <- matrix(runif(n_pixels * n_dims), nrow = n_pixels, ncol = n_dims)
    ref <- sweep(sweep(ref, 2, dim_range, "*"), 2, dim_min, "+")
    for (ci in seq_along(c_range)) {
      null_seed <- seed + 100000L + b * 1000L + c_range[ci]
      log_wk_null[b, ci] <- log(fit_hard_wk(ref, c_range[ci], fuzziness, null_seed)$wk)
    }
  }

  gap <- colMeans(log_wk_null) - log(wk_observed)
  sk  <- apply(log_wk_null, 2, sd) * sqrt(1 + 1 / n_null)

  list(gap_df = tibble(c = c_range, gap = gap, sk = sk), hard_by_c = observed)
}

# Global-max-then-1SE rule, NOT a left-to-right scan (which stops early on
# any dip before the true peak -- see this function's self-test below).
pick_optimal_c <- function(gap_df) {
  gap_df <- gap_df %>% arrange(c)
  i_max <- which.max(gap_df$gap)
  threshold <- gap_df$gap[i_max] - gap_df$sk[i_max]
  candidates <- gap_df$c[gap_df$gap >= threshold]
  min(candidates)
}

compute_spectral_species_richness_fcm <- function(r, veg_mask, c_max = 20,
                                                    fuzziness = 2.0, n_subsample = 2500,
                                                    n_pc = 4, n_reps = 20, seed_base = 0,
                                                    gap_c_max = 8, gap_n_null = 10) {
  # SAME masking + PCA recipe as the removed k-means version -- identical
  # PC-space inputs, only the clustering/search step differs.
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  this_n_subsample <- n_subsample
  if (nrow(vals) < this_n_subsample) this_n_subsample <- nrow(vals)
  if (nrow(vals) <= c_max) {
    return(list(richness = NA_real_, shannon_h = NA_real_, shannon_effective = NA_real_,
                mean_winning_c = NA_real_, max_winning_c = NA_real_,
                n_reps_gap_check_triggered = NA_integer_,
                primary_time_sec = NA_real_, gap_time_sec = NA_real_))
  }

  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pcs_all <- pca$x[, 1:n_pc]

  richness_reps  <- numeric(n_reps)
  shannon_reps   <- numeric(n_reps)
  winning_c_reps <- integer(n_reps)
  gap_triggered  <- logical(n_reps)
  primary_time_sec <- 0
  gap_time_sec     <- 0

  for (rep_i in seq_len(n_reps)) {
    sub_idx <- sample(seq_len(nrow(pcs_all)), this_n_subsample)
    pcs_sub <- pcs_all[sub_idx, , drop = FALSE]

    t0 <- Sys.time()
    sel <- select_winning_c(pcs_sub, c_max = c_max, fuzziness = fuzziness,
                             seed = seed_base + rep_i * 1000L)
    primary_time_sec <- primary_time_sec + as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (sel$winning_c == 2) {
      gap_triggered[rep_i] <- TRUE
      t0 <- Sys.time()
      gap_result <- compute_gap_statistic(pcs_sub, c_range = 1:min(gap_c_max, c_max),
                                           fuzziness = fuzziness,
                                           seed = seed_base + rep_i * 1000L + 500000L,
                                           n_null = gap_n_null)
      gap_time_sec <- gap_time_sec + as.numeric(difftime(Sys.time(), t0, units = "secs"))

      optimal_c <- pick_optimal_c(gap_result$gap_df)
      final_c    <- optimal_c
      final_hard <- gap_result$hard_by_c[[as.character(optimal_c)]]$hard
    } else {
      final_c    <- sel$winning_c
      final_hard <- sel$hard
    }

    winning_c_reps[rep_i] <- final_c
    # final_c (primary argmax_c L(c), or the gap statistic's answer when the
    # c=2-floor fallback triggered) IS this rep's richness value directly --
    # no full-population reassignment step, unlike the removed k-means
    # version (see this function's header caveat).
    richness_reps[rep_i]  <- final_c

    # Spectral Shannon's index (Wang et al. 2018, RSE 211:218-228, Table 2),
    # same formula as the removed k-means version, different population
    # scope -- see this function's header caveat. For c=1 (single group),
    # p_i = 1 and H' = 0, same formula, no special case needed.
    tbl <- table(final_hard)
    p_i <- as.numeric(tbl[tbl > 0]) / length(final_hard)
    shannon_reps[rep_i] <- -sum(p_i * log(p_i))
  }

  list(richness = mean(richness_reps),
       shannon_h = mean(shannon_reps),
       shannon_effective = mean(exp(shannon_reps)),
       mean_winning_c = mean(winning_c_reps),
       max_winning_c = max(winning_c_reps),
       n_reps_gap_check_triggered = sum(gap_triggered),
       primary_time_sec = primary_time_sec,
       gap_time_sec = gap_time_sec)
}

# ============================================================================
# 6. Metric 4: Rao's Q, local moving-window implementation (no global
#    distance matrix)
# ============================================================================
# This used to call rasterdiv::paRao() (see git history / SESSION_LOG.md for
# that version). paRao()'s "classic" and "multidimension" modes both build
# ONE global pairwise distance matrix across every UNIQUE value in the
# ENTIRE input raster, then look up each window's result from it. That
# scales as O(unique_values^2) in memory -- fine for coarse/classified
# rasters, but for continuous ~1m-resolution reflectance/NDVI data a single
# site-year raster has hundreds of thousands of near-unique values, and
# paRao failed outright on ABBY 2017 NDVI (1000x1000, ~750k mostly-unique
# values even after simplify=4) trying to allocate a >1 TB distance matrix.
# That's not a tunable-parameter problem -- paRao's design assumption
# (few unique classes, looked up globally) doesn't match this data.
#
# Rao's Q is inherently local: Rocchini et al. (2017) and the classic
# Botta-Dukat (2005) formulation both define it per moving window, using
# only that window's own pixel values -- a window's Q never needs distances
# to pixel values from elsewhere in the raster. The implementation below
# computes Q per window directly from that window's own (<= window^2) pixel
# vectors, so memory for the pairwise-distance step is bounded by window^2
# regardless of total raster size. This is the same window-local approach
# ComputeSpecBiodiv.R's pyGNDiv-based rao_q_window_pygndiv() /
# compute_rao_q_raster_pygndiv() used before the (now-reverted) rasterdiv
# swap, just without the per-window Python/reticulate round-trip.
#
# Formula: classic (non-normalized) Rao's Q, alpha = 1 (arithmetic-mean
# distance weighting) -- Q = sum_i sum_j p_i * p_j * d_ij, the FULL double
# sum over all pairs (i, j) including i == j (which contributes 0), NOT
# divided by 2 (Botta-Dukat 2005 eq. 1; Rocchini et al. 2017 eq. 1). d_ij is
# Euclidean distance between pixel spectral vectors (1 dimension for
# NDVI/NIRv, n_bands dimensions for all-bands). rao_q_window() below
# computes this as sum(D) / n^2 over the window's own n (<= window^2)
# pixels with equal weight p_i = 1/n each -- algebraically identical to the
# textbook sum_i sum_j p_i*p_j*d_ij (duplicate pixel values need no special
# handling: repeating an identical value k times and weighting it 1/n each
# is the same sum as giving it weight k/n once).
#
# Windows are non-overlapping (stride = window), matching the tiling the
# pyGNDiv version used, and the per-window Q values are averaged for the
# site-year result, matching how the rasterdiv version aggregated
# (mean(values(...), na.rm = TRUE)). Incomplete windows at the raster's
# bottom/right edge are dropped rather than padded, same as the pyGNDiv
# version's row_centers/col_centers behavior.
#
# Vectorization: terra::focal() and terra::focalValues() were evaluated and
# rejected for the all-bands case -- confirmed empirically against a
# synthetic 2-layer SpatRaster that focalValues() silently returns only the
# FIRST layer's window values, and focal() applies its function per-layer
# independently, so neither can hand one function call a pixel's window
# values across ALL bands together, which all-bands Rao's Q needs (distance
# is computed jointly across bands, not layer-by-layer). Instead, every
# window's pixel-cell indices are computed in ONE vectorized
# cellFromRowCol() call up front (not recomputed per window, unlike the old
# pygndiv version's per-window call inside its loop), and feature_r's values
# matrix is extracted once. The remaining per-window loop does O(window^2)
# work on values already resident in memory -- no per-window disk I/O, no
# per-window Python call, and no raster-wide distance matrix. Benchmarked on
# synthetic 1000x1000 rasters: ~5s for a single-layer (NDVI/NIRv) pass, and
# runtime for the all-bands case scales with band count via each window's
# dist() call (~350 bands, 500x500: ~6s; scales roughly linearly with pixel
# count from there) -- on the order of tens of seconds to low minutes per
# site-year for all three Rao's Q variants combined, well under the cost of
# the spectral-species-richness step (adaptive-FCM, n_reps_ssr reps) above.
# ----------------------------------------------------------------------------
# 6a. PCA-based band reduction, ALL-BANDS Rao's Q ONLY.
# ----------------------------------------------------------------------------
# rao_q_window() below runs a pairwise-distance computation per window; fed
# the full ~346-band raster directly, every 3x3 window's dist() call is over
# ~346 dimensions. Reduce to the leading PCs that explain a fixed variance
# threshold first, same convention already used for CHV/SSR (prcomp(),
# center = TRUE, scale. = FALSE) -- N is data-driven per site-year (band
# count varies slightly after bad-band masking), not a fixed count. NDVI and
# NIRv Rao's Q are NOT reduced -- they're already 1-dimensional.
#
# PCA loadings are fit on complete-case (non-NA, i.e. in-veg-mask) pixels
# only, then applied to every pixel via predict.prcomp()'s matrix multiply,
# which propagates NA rows through untouched -- preserving the input
# raster's spatial NA structure (out-of-mask cells stay NA in every output
# PC layer) without a second masking pass.
pca_reduce_raster <- function(r, veg_mask, var_threshold = 0.99) {
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked)   # ncell x nlyr, NA preserved outside the mask
  complete_idx <- stats::complete.cases(vals)
  if (sum(complete_idx) < 2) return(r_masked)

  pca <- prcomp(vals[complete_idx, , drop = FALSE], center = TRUE, scale. = FALSE)
  cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_pc <- which(cum_var >= var_threshold)[1]
  if (is.na(n_pc)) n_pc <- length(pca$sdev)

  scores_all <- predict(pca, newdata = vals)[, seq_len(n_pc), drop = FALSE]
  pc_r <- rast(r_masked, nlyr = n_pc)
  values(pc_r) <- scores_all
  names(pc_r) <- paste0("PC", seq_len(n_pc))
  pc_r
}

rao_q_window <- function(window_vals) {
  window_vals <- window_vals[stats::complete.cases(window_vals), , drop = FALSE]
  n <- nrow(window_vals)
  if (n < 2) return(NA_real_)
  D <- as.matrix(stats::dist(window_vals, method = "euclidean"))
  sum(D) / (n * n)
}

compute_rao_q_raster_local <- function(feature_r, window = 3) {
  stopifnot(window %% 2 == 1)
  half <- window %/% 2
  nr <- nrow(feature_r); nc <- ncol(feature_r)
  if (nr < window || nc < window) return(NA_real_)

  row_centers <- seq(half + 1, nr - half, by = window)
  col_centers <- seq(half + 1, nc - half, by = window)
  n_row <- length(row_centers); n_col <- length(col_centers)
  if (n_row == 0 || n_col == 0) return(NA_real_)

  # offsets of the window^2 positions within a window, relative to its center
  off <- (-half):half
  row_off <- rep(off, times = window)
  col_off <- rep(off, each  = window)

  # every window's center, repeated once per within-window offset -- lets a
  # single vectorized cellFromRowCol() call resolve cell indices for ALL
  # windows at once instead of once per window
  centers_row <- rep(row_centers, times = n_col)
  centers_col <- rep(col_centers, each  = n_row)
  n_windows <- length(centers_row)

  cell_rows <- outer(row_off, centers_row, "+")   # window^2 x n_windows
  cell_cols <- outer(col_off, centers_col, "+")   # window^2 x n_windows
  cell_idx  <- cellFromRowCol(feature_r, as.vector(cell_rows), as.vector(cell_cols))
  dim(cell_idx) <- c(window * window, n_windows)

  vals_all <- values(feature_r)   # ncell x nlyr, extracted once

  raoq_vals <- vapply(seq_len(n_windows), function(w) {
    rao_q_window(vals_all[cell_idx[, w], , drop = FALSE])
  }, numeric(1))

  mean(raoq_vals, na.rm = TRUE)
}

# ============================================================================
# 7. Discover which years of hyperspectral data are actually on disk, per
#    site -- do not assume every site has the same year range, and do not
#    assume the most recent year is representative (mirrors how
#    AnnualFootprintExtent.R discovers available years per site).
# ============================================================================
# byTileAOP() writes tiles under nested product/year/FullSite/domain/
# <year>_<SITECODE>_<visit>/... directories -- the year is not reliably
# encoded in the h5 filename itself, so it's read off the directory
# structure. Prefer the "<year>_<SITECODE>_<visit>" folder (unambiguous),
# falling back to any standalone 4-digit path component in a plausible NEON
# operating year range.
extract_year_from_path <- function(path) {
  parts <- str_split(path, "[/\\\\]")[[1]]
  m <- str_match(parts, "^(\\d{4})_[A-Za-z]{4}_\\d+$")
  yr <- m[!is.na(m[, 2]), 2]
  if (length(yr) >= 1) return(yr[1])
  yr2 <- parts[str_detect(parts, "^(19|20)\\d{2}$")]
  if (length(yr2) >= 1) return(yr2[1])
  NA_character_
}

discover_site_years <- function(tower_id) {
  tile_dir <- file.path(hyperspec_dir, tower_id)
  if (!dir.exists(tile_dir)) return(tibble(file = character(), year = character()))
  h5_files <- list.files(tile_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) return(tibble(file = character(), year = character()))
  tibble(file = h5_files, year = map_chr(h5_files, extract_year_from_path)) %>%
    filter(!is.na(year))
}

site_years_files <- setNames(map(towers_df$Site.ID, discover_site_years), towers_df$Site.ID)

inventory <- imap(site_years_files, function(df, tid) {
  if (nrow(df) == 0) return(NULL)
  tibble(tower_id = tid, year = sort(unique(df$year)))
}) %>% compact() %>% bind_rows()

cat("\n==== Hyperspectral data inventory (years on disk per site) ====\n")
if (nrow(inventory) == 0) {
  cat("  No hyperspectral h5 files found under", hyperspec_dir, "for any site in NEONsites.csv.\n")
  single_year_tower_ids <- character(0)
} else {
  inventory_summary <- inventory %>%
    group_by(tower_id) %>%
    summarise(years = paste(year, collapse = ", "), n_years = n(), .groups = "drop") %>%
    arrange(tower_id)
  print(inventory_summary, n = Inf)

  single_year_tower_ids <- inventory_summary$tower_id[inventory_summary$n_years == 1]
  if (length(single_year_tower_ids) > 0) {
    cat("\n  Single-year sites (no interannual comparison possible -- included in the\n",
        "  output as a baseline/reference point, flagged in `status`):\n   ",
        paste(single_year_tower_ids, collapse = ", "), "\n", sep = "")
  }
}
cat("=================================================================\n\n")

# YEAR-MAJOR ordering: all available sites for the earliest year first, then
# all sites for the next year, and so on -- NOT the tower-major (all years
# for one site, then the next site) order this loop used previously. Reuses
# `inventory` (Section 7 above; identical tower_id/year pairs to
# site_years_files) sorted by year first, tower second, so an interrupted
# weekend run has touched a spread of sites at each year rather than
# exhausting one site's full year range before starting the next site. Each
# job's list(tower_id=, year=, files=) entry is constructed identically to
# before -- only the iteration order changed.
site_year_jobs <- list()
if (nrow(inventory) > 0) {
  inventory_ordered <- inventory %>% arrange(year, tower_id)
  for (i in seq_len(nrow(inventory_ordered))) {
    tid <- inventory_ordered$tower_id[i]
    yr  <- inventory_ordered$year[i]
    df  <- site_years_files[[tid]]
    site_year_jobs[[length(site_year_jobs) + 1]] <- list(
      tower_id = tid, year = yr, files = df$file[df$year == yr]
    )
  }
}

cat("Total site-year jobs to process:", length(site_year_jobs), "\n")
if (length(site_year_jobs) == 0) {
  stop("No site-years with hyperspectral data found under ", hyperspec_dir, " -- nothing to process.")
}

# ============================================================================
# 8. Main loop: one row per tower-YEAR, single 500m buffer for every metric
# ============================================================================
results <- tibble(
  tower_id = character(), neon_site = character(), year = character(),
  cv = double(), chv = double(), chv_standardized = double(), cha = double(),
  spectral_species_richness = double(), shannon_h = double(), shannon_effective = double(),
  ssr_mean_winning_c = double(), ssr_max_winning_c = double(),
  n_reps_gap_check_triggered = integer(),
  raoq_ndvi = double(), raoq_nirv = double(), raoq_allbands = double(),
  status = character()
)

pb <- txtProgressBar(min = 0, max = length(site_year_jobs), style = 3)
job_minutes <- numeric(0)

for (j in seq_along(site_year_jobs)) {
  job <- site_year_jobs[[j]]
  tower_id  <- job$tower_id
  yr        <- job$year
  neon_site <- towers_df$neon_site[towers_df$Site.ID == tower_id][1]
  job_start <- Sys.time()

  cat("\n==== [", j, "/", length(site_year_jobs), "] ", tower_id,
      " (", neon_site, ") - ", yr, " ====\n", sep = "")

  row <- tryCatch({
    data <- get_tower_reflectance(tower_id, job$files, buffer_m)
    if (is.null(data)) {
      tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
             cv = NA_real_, chv = NA_real_, cha = NA_real_,
             spectral_species_richness = NA_real_, shannon_h = NA_real_, shannon_effective = NA_real_,
             ssr_mean_winning_c = NA_real_, ssr_max_winning_c = NA_real_,
             n_reps_gap_check_triggered = NA_integer_,
             raoq_ndvi = NA_real_, raoq_nirv = NA_real_, raoq_allbands = NA_real_,
             status = "no reflectance data")
    } else {
      r <- data$raster; wl <- data$wavelengths

      ndvi <- compute_ndvi_raster(r, wl)
      nirv <- compute_nirv_raster(r, wl)
      # terra::mask(r, veg_mask) only excludes cells where veg_mask is NA
      # (default maskvalues = NA) -- a logical/0-1 raster (e.g. ndvi >
      # ndvi_thresh) has no NAs and masks nothing. Must be 1/NA, not TRUE/FALSE.
      veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)

      cat("  computing CV...\n");  cv_val  <- compute_cv(r, veg_mask)
      cat("  computing CHV...\n"); chv_val <- compute_chv(r, veg_mask, n_pc_chv)
      cat("  computing CHA...\n"); cha_val <- compute_cha(r, veg_mask)
      cat("  computing spectral species richness (adaptive-FCM) + Shannon's H' (", n_reps_ssr, " reps)...\n", sep = "")
      job_seed_base <- 20260101L + j * 1000L
      ssr_result <- compute_spectral_species_richness_fcm(r, veg_mask, c_max_ssr,
                                                    fcm_fuzziness_w, n_subsample_pixels,
                                                    n_pc_ssr, n_reps_ssr, seed_base = job_seed_base,
                                                    gap_c_max = gap_check_c_max, gap_n_null = gap_check_n_null)
      ssr_val             <- ssr_result$richness
      shannon_h_val        <- ssr_result$shannon_h
      shannon_eff_val      <- ssr_result$shannon_effective
      ssr_mean_wc_val      <- ssr_result$mean_winning_c
      ssr_max_wc_val       <- ssr_result$max_winning_c
      n_gap_triggered_val  <- ssr_result$n_reps_gap_check_triggered
      cat(sprintf("    gap-statistic fallback triggered on %d/%d reps -- primary search %.1fs, gap fallback %.1fs\n",
                   n_gap_triggered_val, n_reps_ssr,
                   ssr_result$primary_time_sec, ssr_result$gap_time_sec))

      cat("  computing Rao Q (NDVI), local moving-window...\n")
      raoq_ndvi <- compute_rao_q_raster_local(mask(ndvi, veg_mask), raoq_window)
      cat("  computing Rao Q (NIRv), local moving-window...\n")
      raoq_nirv <- compute_rao_q_raster_local(mask(nirv, veg_mask), raoq_window)
      cat("  computing Rao Q (all bands), local moving-window (PCA-reduced)...\n")
      r_allbands_pca <- pca_reduce_raster(r, veg_mask, raoq_pca_var_threshold)
      cat("    reduced to", nlyr(r_allbands_pca), "PC(s) (>=",
          raoq_pca_var_threshold * 100, "% variance)\n", sep = "")
      raoq_all  <- compute_rao_q_raster_local(r_allbands_pca, raoq_window)

      status <- if (tower_id %in% single_year_tower_ids) "success (single-year site)" else "success"

      tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
             cv = cv_val, chv = chv_val, cha = cha_val,
             spectral_species_richness = ssr_val, shannon_h = shannon_h_val, shannon_effective = shannon_eff_val,
             ssr_mean_winning_c = ssr_mean_wc_val, ssr_max_winning_c = ssr_max_wc_val,
             n_reps_gap_check_triggered = n_gap_triggered_val,
             raoq_ndvi = raoq_ndvi, raoq_nirv = raoq_nirv, raoq_allbands = raoq_all,
             status = status)
    }
  }, error = function(e) {
    cat("  unexpected error:", conditionMessage(e), "\n")
    tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
           cv = NA_real_, chv = NA_real_, cha = NA_real_,
           spectral_species_richness = NA_real_, shannon_h = NA_real_, shannon_effective = NA_real_,
           ssr_mean_winning_c = NA_real_, ssr_max_winning_c = NA_real_,
           n_reps_gap_check_triggered = NA_integer_,
           raoq_ndvi = NA_real_, raoq_nirv = NA_real_, raoq_allbands = NA_real_,
           status = paste("error:", conditionMessage(e)))
  })

  results <- bind_rows(results, row)

  # CHV z-score across all site-YEAR observations completed SO FAR (not just
  # across sites, and not only computed once at the very end) -- recomputing
  # this on every completed job means an interrupted run still leaves a
  # best-available standardization in the csv, and the final write reflects
  # the full site-year sample once the loop finishes.
  results <- results %>%
    mutate(chv_standardized = (chv - mean(chv, na.rm = TRUE)) / sd(chv, na.rm = TRUE)) %>%
    select(tower_id, neon_site, year, cv, chv, chv_standardized, cha,
           spectral_species_richness, shannon_h, shannon_effective,
           ssr_mean_winning_c, ssr_max_winning_c, n_reps_gap_check_triggered,
           raoq_ndvi, raoq_nirv, raoq_allbands, status)

  # save incrementally so a crash/interruption partway through doesn't lose
  # completed site-years (matches AnnualFootprintExtent.R)
  write.csv(results, out_csv, row.names = FALSE)

  elapsed_min <- as.numeric(difftime(Sys.time(), job_start, units = "mins"))
  job_minutes <- c(job_minutes, elapsed_min)
  avg_min <- mean(job_minutes)
  remaining <- length(site_year_jobs) - j
  cat(sprintf("  done in %.1f min (avg %.1f min/job, ~%.0f min remaining for %d job(s))\n",
              elapsed_min, avg_min, avg_min * remaining, remaining))

  setTxtProgressBar(pb, j)
}

close(pb)
cat("\nDone. Annual spectral diversity written to:", out_csv, "\n")
