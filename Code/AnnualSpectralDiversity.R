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
#   3. Spectral Species Richness (RF + K-means, 50 clusters)
#   4. Rao's Q (NDVI, NIRv, all-bands), via rasterdiv (CRAN), classic
#      (non-normalized) form -- NOT via pyGNDiv/reticulate. Comparability
#      across variable-set sizes (e.g. a 1-band NDVI Q vs. an all-bands Q,
#      which is what pyGNDiv's generalizable normalization was for) is not
#      needed here: every comparison in this analysis is across YEARS at the
#      same site, using the same variable set each time. A hypothetical
#      future analysis that compared Rao's Q *across index types* directly
#      would need that normalization back; this one doesn't.
#
# Differences from ComputeSpecBiodiv.R:
#   - a single fixed 500m buffer is used for ALL FOUR metrics. The original
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
# uses rasterdiv instead of ComputeSpecBiodiv.R's pyGNDiv approach; see the
# note above and in Section 6 for why.
#
# Runtime: dominated by spectral species richness (RF + K-means, n_reps_ssr
# reps per site-year) and the three Rao's Q rasterdiv calls (one native
# moving-window pass per variant, no per-window Python round-trip). The 500m
# buffer here is smaller than ComputeSpecBiodiv.R's CV/CHV/SSR buffer (530m,
# so that part of the cost is about the same) but larger than its Rao Q
# buffer (300m) -- area scales as radius^2, so each Rao Q variant here covers
# roughly (500/300)^2 ~= 2.8x as many pixels/windows as the original script's
# Rao Q step. Total run time scales with the number of site-YEARS, not the
# number of sites: a site with 4 years on disk costs ~4x a single-year site.
# Rather than guess an absolute number, the loop below times and prints each
# completed site-year and a running average/ETA, so you can extrapolate for
# your actual job count.

library(terra)
library(rhdf5)          # low-level NEON H5 reading
library(geometry)        # convhulln for CHV
library(randomForest)
library(cluster)
library(dplyr)
library(stringr)
library(purrr)
library(rasterdiv)
library(tibble)

# ============================================================================
# 0. Setup
# ============================================================================
hyperspec_dir <- "X:/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
out_csv       <- "X:/moore/SpectralBiodiversity/Data/spectral_diversity_by_year.csv"

buffer_m           <- 500   # single buffer for ALL metrics (CV, CHV, SSR, RaoQ x3)
ndvi_thresh         <- 0.4
n_clusters           <- 50
n_subsample_pixels   <- 2500
n_pc_ssr             <- 4
n_pc_chv             <- 3
n_reps_ssr           <- 20
raoq_window          <- 3   # side of the square moving window passed to rasterdiv::paRao()

towers_df <- read.csv("./Data/NEONsites.csv") %>%
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
# 5. Metric 3: Spectral Species Richness (reused from ComputeSpecBiodiv.R)
# ============================================================================
compute_spectral_species_richness <- function(r, veg_mask, n_clusters = 50,
                                              n_subsample = 2500, n_pc = 4,
                                              n_reps = 20) {
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  if (nrow(vals) < n_subsample) n_subsample <- nrow(vals)
  if (nrow(vals) < n_clusters) return(NA_real_)

  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pcs_all <- pca$x[, 1:n_pc]

  richness_reps <- numeric(n_reps)
  for (rep_i in seq_len(n_reps)) {
    sub_idx <- sample(seq_len(nrow(pcs_all)), n_subsample)
    pcs_sub <- pcs_all[sub_idx, ]

    rf <- randomForest(x = pcs_sub, ntree = 500, proximity = TRUE)
    prox_dist <- as.dist(1 - rf$proximity)
    km <- kmeans(cmdscale(prox_dist, k = n_pc), centers = n_clusters, nstart = 10)

    classifier <- randomForest(x = pcs_sub, y = as.factor(km$cluster), ntree = 500)
    pred_all <- predict(classifier, newdata = pcs_all)

    richness_reps[rep_i] <- length(unique(pred_all))
  }
  mean(richness_reps)
}

# ============================================================================
# 6. Metric 4: Rao's Q via rasterdiv (classic, non-normalized)
# ============================================================================
# ComputeSpecBiodiv.R computes Rao's Q via pyGNDiv/reticulate, whose
# "generalizable normalization" exists to make Q comparable across variable
# sets of different sizes (e.g. a 1-band NDVI-derived Q vs. an all-bands Q).
# This script doesn't need that: every comparison here is across YEARS at
# the same site, using the SAME variable set (NDVI-only, NIRv-only, or
# all-bands) each time. So Rao's Q here is computed with rasterdiv (CRAN;
# Rocchini et al. 2017, Ecological Indicators, describes the original
# spectralrao() function this package grew out of) via its current public
# function, paRao(), in the classic (non-normalized) parametric form. This
# drops the Python/reticulate dependency entirely -- verified against the
# installed rasterdiv CRAN package (not assumed from the 2017 paper's older
# spectralrao() signature, which no longer applies).
#
# paRao() takes a SpatRaster (or a list of them) directly, no conversion
# needed. method = "classic" is used for single-layer inputs (NDVI, NIRv);
# method = "multidimension" for the all-bands case, which paRao requires as
# a LIST of single-layer SpatRasters rather than one multi-band SpatRaster.
# alpha = 1 (arithmetic-mean distance weighting) is the classic Rao's Q
# formulation. dist_m = "euclidean" and window = raoq_window (3x3) match
# what the pyGNDiv version used.
#
# simplify controls how many decimal places of the input are kept before
# rasterdiv rounds and casts values to integer internally for its distance
# calculations. Its DEFAULT (simplify = 0) rounds to whole numbers, which
# silently collapses reflectance/NDVI/NIRv values (all roughly in a 0-1
# range) to 0 and zeroes out Rao's Q entirely -- confirmed with a synthetic
# raster during validation. simplify = 4 preserves ~4 decimal places, in
# line with the original NEON int16 reflectance scale factor (1/10000), and
# was confirmed to produce non-trivial, non-zero Q values on the same
# synthetic raster.
#
# NA handling: veg_mask-excluded pixels arrive as NA. rasterdiv emits a
# warning claiming NAs "will be treated as 0s", but that text is stale for
# this code path -- its internal distance function propagates NA per pixel
# pair and the window aggregation sums with na.rm = TRUE, so NA pixels are
# excluded from the distance calculation, not substituted with a value. This
# was verified empirically: a homogeneous synthetic raster with an NA hole
# in the middle returns near-zero Q everywhere (as expected for a
# genuinely homogeneous raster), not the large spurious Q that would appear
# at the hole's edge if NA were actually being treated as 0.
compute_rao_q_raster_rasterdiv <- function(feature_r, window = 3, dist_m = "euclidean", simplify = 4) {
  if (nlyr(feature_r) == 1) {
    res <- paRao(x = feature_r, window = window, alpha = 1, method = "classic",
                dist_m = dist_m, na.tolerance = 1, rasterOut = TRUE,
                simplify = simplify, np = 1, progBar = FALSE)
  } else {
    band_list <- lapply(seq_len(nlyr(feature_r)), function(i) feature_r[[i]])
    res <- paRao(x = band_list, window = window, alpha = 1, method = "multidimension",
                dist_m = dist_m, na.tolerance = 1, rasterOut = TRUE,
                simplify = simplify, np = 1, progBar = FALSE)
  }
  mean(values(res[[1]][[1]]), na.rm = TRUE)
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

site_year_jobs <- list()
for (tid in names(site_years_files)) {
  df <- site_years_files[[tid]]
  if (nrow(df) == 0) next
  for (yr in sort(unique(df$year))) {
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
  cv = double(), chv = double(), chv_standardized = double(),
  spectral_species_richness = double(),
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
             cv = NA_real_, chv = NA_real_, spectral_species_richness = NA_real_,
             raoq_ndvi = NA_real_, raoq_nirv = NA_real_, raoq_allbands = NA_real_,
             status = "no reflectance data")
    } else {
      r <- data$raster; wl <- data$wavelengths

      ndvi <- compute_ndvi_raster(r, wl)
      nirv <- compute_nirv_raster(r, wl)
      veg_mask <- ndvi > ndvi_thresh

      cat("  computing CV...\n");  cv_val  <- compute_cv(r, veg_mask)
      cat("  computing CHV...\n"); chv_val <- compute_chv(r, veg_mask, n_pc_chv)
      cat("  computing spectral species richness (", n_reps_ssr, " reps)...\n", sep = "")
      ssr_val <- compute_spectral_species_richness(r, veg_mask, n_clusters,
                                                    n_subsample_pixels, n_pc_ssr, n_reps_ssr)

      cat("  computing Rao Q (NDVI) via rasterdiv...\n")
      raoq_ndvi <- compute_rao_q_raster_rasterdiv(mask(ndvi, veg_mask), raoq_window)
      cat("  computing Rao Q (NIRv) via rasterdiv...\n")
      raoq_nirv <- compute_rao_q_raster_rasterdiv(mask(nirv, veg_mask), raoq_window)
      cat("  computing Rao Q (all bands) via rasterdiv...\n")
      raoq_all  <- compute_rao_q_raster_rasterdiv(mask(r, veg_mask), raoq_window)

      status <- if (tower_id %in% single_year_tower_ids) "success (single-year site)" else "success"

      tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
             cv = cv_val, chv = chv_val, spectral_species_richness = ssr_val,
             raoq_ndvi = raoq_ndvi, raoq_nirv = raoq_nirv, raoq_allbands = raoq_all,
             status = status)
    }
  }, error = function(e) {
    cat("  unexpected error:", conditionMessage(e), "\n")
    tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
           cv = NA_real_, chv = NA_real_, spectral_species_richness = NA_real_,
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
    select(tower_id, neon_site, year, cv, chv, chv_standardized,
           spectral_species_richness, raoq_ndvi, raoq_nirv, raoq_allbands, status)

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
