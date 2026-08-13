# process_footprints.R
# For each downloaded site: builds a footprint raster stack, keeps only
# daytime (solar elevation > 0) half-hourly intervals, computes the
# footprint extent (radius of the 80% cumulative source-area contour)
# for each interval, and averages those into one value per tower.

library(neonUtilities)
library(terra)
library(suncalc)
library(dplyr)
library(stringr)
library(purrr)

# ---- 1. Setup ----------------------------------------------------------
data_dir  <- "X:/moore/SpectralBiodiversity/Data/NEON_Flux"
out_csv   <- "X:/moore/SpectralBiodiversity/Data/NEON_Flux/avg_daytime_footprint_extent.csv"
contour   <- 0.80   # cumulative footprint contour used to define source area (eg extent)

# NEON force-publishes a monthly HDF5 file with valid structure/timestamps but
# every footprint cell set to NaN when source data isn't received within 30 days
# of the processing window. footRaster() reads these without error, so they must
# be detected by content, not by read failure. A file is treated as a placeholder
# if even its most-complete layer is below this fraction of non-NA cells.
placeholder_frac_threshold <- 0.01

towers_df <- read.csv("./Data/NEONsites.csv") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

# ---- 2. Helper: extent of the X% cumulative footprint contour ---------
# Footprint rasters from footRaster() are tower-centered (tower at x=0, y=0).
# For a given layer, sort cells by footprint value, accumulate until the
# contour threshold is reached, and return the max distance from the tower
# among included cells.
footprint_extent <- function(r, threshold = 0.90) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "val"
  if (nrow(df) == 0 || sum(df$val, na.rm = TRUE) == 0) return(NA_real_)
  df$dist <- sqrt(df$x^2 + df$y^2)
  df <- df[order(-df$val), ]
  df$cum <- cumsum(df$val) / sum(df$val)
  keep <- df[df$cum <= threshold, ]
  if (nrow(keep) == 0) keep <- df[1, ]
  max(keep$dist)
}

# ---- 3. Helper: parse a timestamp out of a footRaster layer name ------
# footRaster layer names are the full internal HDF5 path, e.g.:
#   ".../dp04/data/foot/grid/turb/20250101T003000Z"
# The timestamp is always the last 8-digit-date + "T" + 6-digit-time segment.
parse_layer_time <- function(nm) {
  ts <- str_extract(nm, "\\d{8}T\\d{6}(?=Z?$)")
  if (is.na(ts)) return(as.POSIXct(NA))
  as.POSIXct(ts, format = "%Y%m%dT%H%M%S", tz = "UTC")
}

# ---- 4. Process each site (most recent year only) -----------------------
results <- tibble(
  tower_id = character(), neon_site = character(), year = character(),
  avg_daytime_extent_m = numeric(), n_daytime_intervals = integer(),
  n_files_skipped = integer(), status = character()
)

pb <- txtProgressBar(min = 0, max = nrow(towers_df), style = 3)

for (i in seq_len(nrow(towers_df))) {
  
  tower_id  <- towers_df$Site.ID[i]
  neon_site <- towers_df$neon_site[i]
  lat       <- towers_df$Lat[i]
  lon       <- towers_df$Lon[i]
  site_dir  <- file.path(data_dir, tower_id)
  
  cat("\n==== [", i, "/", nrow(towers_df), "] ", tower_id,
      " (", neon_site, ") ====\n", sep = "")
  
  if (!dir.exists(site_dir) || length(list.files(site_dir, recursive = TRUE)) == 0) {
    cat("  No downloaded data found - skipping.\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = NA,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = NA, status = "no data downloaded")
    setTxtProgressBar(pb, i)
    next
  }
  
  all_files <- list.files(site_dir, pattern = "\\.h5$|\\.zip$",
                          recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) {
    cat("  No .h5 or .zip files found under", site_dir, "- skipping.\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = NA,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = NA, status = "no data files found")
    setTxtProgressBar(pb, i)
    next
  }
  
  file_years <- str_extract(basename(all_files), "(?<=\\.)\\d{4}(?=-\\d{2}\\.)")
  if (all(is.na(file_years))) {
    cat("  Could not parse year from filenames - check naming pattern.\n")
    cat("  Example:", basename(all_files[1]), "\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = NA,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = NA, status = "year parsing failed")
    setTxtProgressBar(pb, i)
    next
  }
  
  latest_year <- max(file_years, na.rm = TRUE)
  data_files <- all_files[!is.na(file_years) & file_years == latest_year]
  cat("  most recent year available:", latest_year,
      "(", length(data_files), "of", length(all_files), "files )\n")
  
  # store each footRaster() result wrapped -- protects against terra's
  # external-pointer invalidation when SpatRasters are stashed in a list
  # across loop iterations (values look fine at creation, then silently
  # vanish later without wrap()/unwrap())
  good_stacks <- list()
  n_skipped <- 0
  n_files <- length(data_files)
  for (f_idx in seq_along(data_files)) {
    df <- data_files[f_idx]
    cat("\r  reading file", f_idx, "of", n_files, "...          ")
    st <- tryCatch({
      raw_st <- footRaster(filepath = df, progress = FALSE)
      frac_valid <- max(terra::global(raw_st, fun = "notNA")$notNA) / terra::ncell(raw_st)
      if (frac_valid < placeholder_frac_threshold) {
        stop("placeholder file (footprint data is NaN-filled)")
      }
      terra::wrap(raw_st)
    }, error = function(e) {
      cat("\n  skipping", basename(df), "-", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(st)) good_stacks[[length(good_stacks) + 1]] <- st
    else n_skipped <- n_skipped + 1
  }
  cat("\r  ", length(good_stacks), "of", n_files,
      "files read successfully (", n_skipped, "skipped)          \n")
  
  if (length(good_stacks) == 0) {
    cat("  All files failed to read - skipping site.\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = latest_year,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = n_skipped, status = "all files failed")
    setTxtProgressBar(pb, i)
    next
  }
  
  # unwrap before combining -- restores real SpatRaster objects from the
  # protected wrapper representation
  r_stack <- rast(terra::unwrap(good_stacks[[1]]))
  if (length(good_stacks) > 1) {
    for (k in 2:length(good_stacks)) r_stack <- c(r_stack, rast(terra::unwrap(good_stacks[[k]])))
  }
  lyr_names <- names(r_stack)
  keep_idx <- which(!grepl("^mean$", lyr_names, ignore.case = TRUE))
  r_stack <- r_stack[[keep_idx]]
  lyr_names <- lyr_names[keep_idx]
  
  times <- map(lyr_names, parse_layer_time) %>% reduce(c)
  valid <- !is.na(times)
  
  if (sum(valid) == 0) {
    cat("  Could not parse timestamps from layer names - check names(r_stack) manually.\n")
    cat("  Example names:", paste(head(lyr_names, 3), collapse = ", "), "\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = latest_year,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = n_skipped, status = "timestamp parsing failed")
    setTxtProgressBar(pb, i)
    next
  }
  
  sun_pos <- getSunlightPosition(date = times[valid], lat = lat, lon = lon)
  is_day  <- !is.na(sun_pos$altitude) & sun_pos$altitude > 0
  
  day_idx <- which(valid)[is_day]
  if (length(day_idx) == 0) {
    cat("  No daytime intervals identified - skipping.\n")
    results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = latest_year,
                       avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                       n_files_skipped = n_skipped, status = "no daytime intervals")
    setTxtProgressBar(pb, i)
    next
  }
  
  cat("  averaging", length(day_idx), "daytime footprints...\n")
  daytime_stack <- r_stack[[day_idx]]
  mean_daytime_footprint <- mean(daytime_stack, na.rm = TRUE)
  avg_extent <- footprint_extent(mean_daytime_footprint, threshold = contour)
  n_intervals <- length(day_idx)
  
  cat("  ", n_intervals, "daytime intervals averaged; ", contour * 100,
      "% footprint extent =", round(avg_extent, 1), "m\n")
  
  results <- add_row(results, tower_id = tower_id, neon_site = neon_site, year = latest_year,
                     avg_daytime_extent_m = avg_extent,
                     n_daytime_intervals = n_intervals,
                     n_files_skipped = n_skipped, status = "success")
  
  setTxtProgressBar(pb, i)
}

close(pb)
# ---- 5. Save summary CSV ------------------------------------------------
write.csv(results, out_csv, row.names = FALSE)
cat("\nDone. Summary written to:", out_csv, "\n")