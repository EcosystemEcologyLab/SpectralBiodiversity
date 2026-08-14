# analyze_spectral_diversity.R
#
# Computes, for each NEON tower's hyperspectral footprint (tiles already
# downloaded via download_hyperspec_data.R), one value per tower for:
#   1. Coefficient of Variation (CV)                     -- 530m buffer
#   2. Convex Hull Volume (CHV), first 3 PCs              -- 530m buffer
#   3. Spectral Species Richness (RF + K-means, 50 clusters) -- 530m buffer
#   4. Rao's Q (NDVI, NIRv, all-bands), via pyGNDiv        -- 300m buffer
#
# Rao's Q uses the actual pyGNDiv.py functions (mpaRaoS_freq +
# max_eucl_dist_PCA_fun) via reticulate, applying the "generalizable
# normalization" from Pacheco-Labrador et al. (2023), so that RaoQ values
# computed from 1-variable indices (NDVI, NIRv) and from the full band set
# are on a comparable scale.

library(terra)
library(rhdf5)          # low-level NEON H5 reading
library(geometry)        # convhulln for CHV
library(randomForest)
library(cluster)
library(dplyr)
library(stringr)
library(purrr)
library(reticulate)

# ============================================================================
# 0. Setup
# ============================================================================
hyperspec_dir <- "X:/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
out_csv       <- "X:/moore/SpectralBiodiversity/Data/spectral_diversity_summary.csv"
pygndiv_path  <- "./Code/pyGNDiv.py"   # place the uploaded pyGNDiv.py here

buffer_cv_chv_ssr_m <- 530   # matches your hyperspec download buffer
buffer_raoq_m       <- 300   # matches Rossi/Pacheco-Labrador-style radius
ndvi_thresh          <- 0.4
n_clusters            <- 50
n_subsample_pixels    <- 2500
n_pc_ssr              <- 4
n_pc_chv              <- 3
n_reps_ssr             <- 20
raoq_window            <- 3
raoq_stride            <- raoq_window  # non-overlapping by default; see notes

towers_df <- read.csv("./Data/NEONsites.csv") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites to process.\n")

# ============================================================================
# 1. NEON reflectance reading helpers
# ============================================================================
# NEON DP3.30006.00x H5 tiles store reflectance as int16 scaled by a factor
# (typically 10000) with a fill/no-data value (typically -9999 or -100).
# Wavelengths (nm) are stored as a dataset alongside the reflectance cube.
# Bad bands are the standard water-vapor absorption + edge regions widely
# used in NEON hyperspectral tutorials -- verify against your product
# version's metadata before trusting these ranges blindly.
bad_band_ranges <- list(c(300, 400), c(1340, 1445), c(1790, 1955), c(2400, 2600))

read_neon_h5_tile <- function(h5_path) {
  refl_path <- h5ls(h5_path) %>% filter(str_detect(name, "Reflectance$")) %>% pull(group) %>% unique()
  site_group <- str_split(refl_path, "/")[[1]][2]
  
  wavelengths <- h5read(h5_path, paste0("/", site_group, "/Reflectance/Metadata/Spectral_Data/Wavelength"))
  refl <- h5read(h5_path, paste0("/", site_group, "/Reflectance/Reflectance_Data"))
  attrs <- h5readAttributes(h5_path, paste0("/", site_group, "/Reflectance/Reflectance_Data"))
  scale_factor <- ifelse(!is.null(attrs$Scale_Factor), attrs$Scale_Factor, 10000)
  no_data <- ifelse(!is.null(attrs$Data_Ignore_Value), attrs$Data_Ignore_Value, -9999)
  
  map_info <- h5read(h5_path, paste0("/", site_group, "/Reflectance/Metadata/Coordinate_System/Map_Info"))
  map_info <- str_split(map_info, ",")[[1]]
  px_size  <- as.numeric(map_info[6])
  x_min    <- as.numeric(map_info[4])
  y_max    <- as.numeric(map_info[5])
  epsg_code <- h5read(h5_path, paste0("/", site_group, "/Reflectance/Metadata/Coordinate_System/EPSG Code"))

  refl[refl == no_data] <- NA
  refl <- refl / scale_factor   # -> 0-1 reflectance

  # H5 array dims are (band, col, row); reorder to (row, col, band) for terra
  refl <- aperm(refl, c(3, 2, 1))
  nrow_r <- dim(refl)[1]; ncol_r <- dim(refl)[2]

  r <- rast(refl, extent = ext(x_min, x_min + ncol_r * px_size,
                               y_max - nrow_r * px_size, y_max))
  crs(r) <- paste0("EPSG:", epsg_code)
  names(r) <- paste0("b", round(wavelengths))

  # rhdf5 accumulates open file/dataset identifiers across repeated
  # h5read()/h5readAttributes() calls; across many site loop iterations
  # this exhausts HDF5's handle table and throws an internal
  # H5IdComponent error. Close everything this call opened before returning.
  h5closeAll()

  list(raster = r, wavelengths = wavelengths)
}

# drop water-absorption / edge bands
mask_bad_bands <- function(r, wavelengths) {
  bad <- rep(FALSE, length(wavelengths))
  for (rng in bad_band_ranges) bad <- bad | (wavelengths >= rng[1] & wavelengths <= rng[2])
  list(raster = r[[!bad]], wavelengths = wavelengths[!bad])
}

# read + mosaic all tiles for one tower, crop to a circular buffer
get_tower_reflectance <- function(tower_id, buffer_m) {
  tile_dir <- file.path(hyperspec_dir, tower_id)
  h5_files <- list.files(tile_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) return(NULL)
  
  tiles <- map(h5_files, possibly(read_neon_h5_tile, otherwise = NULL))
  tiles <- compact(tiles)
  if (length(tiles) == 0) return(NULL)
  
  wavelengths <- tiles[[1]]$wavelengths
  rasters <- map(tiles, "raster")
  mosaic_r <- if (length(rasters) > 1) do.call(terra::mosaic, rasters) else rasters[[1]]
  
  masked <- mask_bad_bands(mosaic_r, wavelengths)
  
  tower_row <- towers_df %>% filter(Site.ID == tower_id)
  pt <- vect(cbind(tower_row$Lon, tower_row$Lat), crs = "EPSG:4326")
  pt_proj <- project(pt, crs(masked$raster))
  buf <- buffer(pt_proj, buffer_m)
  
  cropped <- crop(mask(masked$raster, buf), buf)
  list(raster = cropped, wavelengths = masked$wavelengths)
}

# ============================================================================
# 2. Spectral index helpers
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
# 3. Metric 1: Coefficient of Variation
# ============================================================================
compute_cv <- function(r, veg_mask) {
  r_masked <- mask(r, veg_mask)
  band_cv <- global(r_masked, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  mean(band_cv[, 1], na.rm = TRUE)
}

# ============================================================================
# 4. Metric 2: Convex Hull Volume
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
# 5. Metric 3: Spectral Species Richness
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
# 6. Metric 4: Rao's Q via pyGNDiv (reticulate)
# ============================================================================
source_python(pygndiv_path)
np <- import("numpy")

# Generalizable-normalization max distance for a variable set of n_vars
# columns (n_vars = 1 for NDVI/NIRv, n_vars = n_bands for all-bands).
get_pca_max_dist <- function(X, n_sigmas = 6, fexp_var = 0.98) {
  n_vars <- ncol(X)
  max_eucl_dist_PCA_fun(as.integer(n_vars), fexp_var = fexp_var, n_sigmas = n_sigmas)
}

# Rao Q for ONE window, via pyGNDiv's mpaRaoS_freq, using only the unique
# pixel vectors within that window (keeps combinations_with_replacement cheap
# -- at most window^2 unique rows, not the whole raster).
rao_q_window_pygndiv <- function(window_vals, max_dist_) {
  window_vals <- window_vals[complete.cases(window_vals), , drop = FALSE]
  n <- nrow(window_vals)
  if (n < 2) return(NA_real_)
  
  uniq <- unique(window_vals)
  freq_row <- sapply(seq_len(nrow(uniq)), function(i) {
    mean(apply(window_vals, 1, function(rw) all(rw == uniq[i, ])))
  })
  freq_mat <- matrix(freq_row, nrow = 1)
  
  X_py    <- np$array(uniq, dtype = "float64")
  freq_py <- np$array(freq_mat, dtype = "float64")
  
  res <- mpaRaoS_freq(X_py, freq_py, alphas_ = list(1L),
                      normalize_dist = TRUE, max_dist_ = max_dist_)
  raoQ <- res[[1]]
  as.numeric(py_to_r(raoQ))
}

# Moving-window Rao Q over a raster (1 or multiple layers). stride = window
# (non-overlapping) is the default for runtime reasons -- a 300m-radius area
# at ~1m resolution is ~280k pixels, and each window is a separate reticulate
# call. Reduce `stride` for denser coverage at the cost of much longer runs.
compute_rao_q_raster_pygndiv <- function(feature_r, window = 3, stride = window) {
  vals_all <- if (nlyr(feature_r) == 1) matrix(values(feature_r), ncol = 1) else values(feature_r)
  max_dist_ <- get_pca_max_dist(vals_all)
  
  nr <- nrow(feature_r); nc <- ncol(feature_r)
  half <- window %/% 2
  row_centers <- seq(half + 1, nr - half, by = stride)
  col_centers <- seq(half + 1, nc - half, by = stride)
  
  raoq_vals <- c()
  for (i in row_centers) {
    for (j in col_centers) {
      cell_idx <- cellFromRowCol(feature_r,
                                 rep((i - half):(i + half), each = window),
                                 rep((j - half):(j + half), times = window))
      w <- vals_all[cell_idx, , drop = FALSE]
      rq <- rao_q_window_pygndiv(w, max_dist_)
      if (!is.na(rq)) raoq_vals <- c(raoq_vals, rq)
    }
  }
  mean(raoq_vals, na.rm = TRUE)
}

# ============================================================================
# 7. Main loop
# ============================================================================
results <- tibble(
  tower_id = character(), neon_site = character(),
  cv = double(), chv = double(), spectral_species_richness = double(),
  raoq_ndvi = double(), raoq_nirv = double(), raoq_allbands = double(),
  status = character()
)

for (i in seq_len(nrow(towers_df))) {
  
  tower_id  <- towers_df$Site.ID[i]
  neon_site <- towers_df$neon_site[i]
  cat("\n==== [", i, "/", nrow(towers_df), "] ", tower_id,
      " (", neon_site, ") ====\n", sep = "")
  
  # ---- 530m buffer: CV, CHV, spectral species richness ----
  data_530 <- tryCatch(get_tower_reflectance(tower_id, buffer_cv_chv_ssr_m),
                       error = function(e) { cat("  read failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(data_530)) {
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site,
                       cv = NA, chv = NA, spectral_species_richness = NA,
                       raoq_ndvi = NA, raoq_nirv = NA, raoq_allbands = NA,
                       status = "no reflectance data")
    next
  }
  r <- data_530$raster; wl <- data_530$wavelengths
  
  ndvi <- compute_ndvi_raster(r, wl)
  veg_mask <- ndvi > ndvi_thresh
  
  cat("  computing CV...\n");  cv_val  <- compute_cv(r, veg_mask)
  cat("  computing CHV...\n"); chv_val <- compute_chv(r, veg_mask, n_pc_chv)
  cat("  computing spectral species richness (", n_reps_ssr, "reps)...\n")
  ssr_val <- compute_spectral_species_richness(r, veg_mask, n_clusters,
                                               n_subsample_pixels, n_pc_ssr, n_reps_ssr)
  
  # ---- 300m buffer: Rao Q (NDVI, NIRv, all-bands) via pyGNDiv ----
  data_300 <- tryCatch(get_tower_reflectance(tower_id, buffer_raoq_m),
                       error = function(e) NULL)
  if (!is.null(data_300)) {
    r300 <- data_300$raster; wl300 <- data_300$wavelengths
    ndvi300 <- compute_ndvi_raster(r300, wl300)
    nirv300 <- compute_nirv_raster(r300, wl300)
    veg300  <- ndvi300 > ndvi_thresh
    
    cat("  computing Rao Q (NDVI) via pyGNDiv...\n")
    raoq_ndvi <- compute_rao_q_raster_pygndiv(mask(ndvi300, veg300), raoq_window, raoq_stride)
    cat("  computing Rao Q (NIRv) via pyGNDiv...\n")
    raoq_nirv <- compute_rao_q_raster_pygndiv(mask(nirv300, veg300), raoq_window, raoq_stride)
    cat("  computing Rao Q (all bands) via pyGNDiv...\n")
    raoq_all  <- compute_rao_q_raster_pygndiv(mask(r300, veg300), raoq_window, raoq_stride)
  } else {
    raoq_ndvi <- NA; raoq_nirv <- NA; raoq_all <- NA
  }
  
  results <- add_row(results, tower_id = tower_id, neon_site = neon_site,
                     cv = cv_val, chv = chv_val, spectral_species_richness = ssr_val,
                     raoq_ndvi = raoq_ndvi, raoq_nirv = raoq_nirv, raoq_allbands = raoq_all,
                     status = "success")
}

# ============================================================================
# 8. Standardize CHV across towers (z-score) and save
# ============================================================================
results <- results %>%
  mutate(chv_standardized = (chv - mean(chv, na.rm = TRUE)) / sd(chv, na.rm = TRUE))

write.csv(results, out_csv, row.names = FALSE)
cat("\nDone. Summary written to:", out_csv, "\n")