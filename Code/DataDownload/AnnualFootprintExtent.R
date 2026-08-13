# process_footprints_annual.R
# Same footprint-extent logic as process_footprints.R, but computes one
# value PER SITE PER YEAR instead of restricting to the most recent year.
# Meant to run in the background / overnight, since it processes every
# available year across all 45 sites.

library(neonUtilities)
library(terra)
library(suncalc)
library(dplyr)
library(stringr)
library(purrr)

# ---- 1. Setup ----------------------------------------------------------
data_dir  <- "X:/moore/SpectralBiodiversity/Data/NEON_Flux"
out_csv   <- "X:/moore/SpectralBiodiversity/Data/NEON_Flux/annual_daytime_footprint_extent.csv"
contour   <- 0.80

towers_df <- read.csv("./Data/NEONsites.csv") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

# ---- 2. Helpers (identical to process_footprints.R) ---------------------
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

parse_layer_time <- function(nm) {
  ts <- str_extract(nm, "\\d{8}T\\d{6}(?=Z?$)")
  if (is.na(ts)) return(as.POSIXct(NA))
  as.POSIXct(ts, format = "%Y%m%dT%H%M%S", tz = "UTC")
}

# process one site-year's worth of files -> one result row
# process one site-year's worth of files -> one result row
process_site_year <- function(tower_id, neon_site, lat, lon, yr, data_files) {
  good_stacks <- list()
  n_skipped <- 0
  for (df in data_files) {
    st <- tryCatch({
      raw_st <- footRaster(filepath = df, progress = FALSE)
      terra::wrap(raw_st)
    }, error = function(e) NULL)
    if (!is.null(st)) good_stacks[[length(good_stacks) + 1]] <- st
    else n_skipped <- n_skipped + 1
  }
  
  if (length(good_stacks) == 0) {
    return(tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
                  avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                  n_files_skipped = n_skipped, status = "all files failed"))
  }
  
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
    return(tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
                  avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                  n_files_skipped = n_skipped, status = "timestamp parsing failed"))
  }
  
  sun_pos <- getSunlightPosition(date = times[valid], lat = lat, lon = lon)
  is_day  <- !is.na(sun_pos$altitude) & sun_pos$altitude > 0
  day_idx <- which(valid)[is_day]
  
  if (length(day_idx) == 0) {
    return(tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
                  avg_daytime_extent_m = NA, n_daytime_intervals = 0,
                  n_files_skipped = n_skipped, status = "no daytime intervals"))
  }
  
  daytime_stack <- r_stack[[day_idx]]
  mean_daytime_footprint <- mean(daytime_stack, na.rm = TRUE)
  avg_extent <- footprint_extent(mean_daytime_footprint, threshold = contour)
  
  tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
         avg_daytime_extent_m = avg_extent, n_daytime_intervals = length(day_idx),
         n_files_skipped = n_skipped, status = "success")
}
# ---- 3. Build the full list of site-year jobs ---------------------------
site_year_jobs <- list()
for (i in seq_len(nrow(towers_df))) {
  tower_id <- towers_df$Site.ID[i]
  site_dir <- file.path(data_dir, tower_id)
  if (!dir.exists(site_dir)) next
  
  all_files <- list.files(site_dir, pattern = "\\.h5$|\\.zip$",
                          recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) next
  
  file_years <- str_extract(basename(all_files), "(?<=\\.)\\d{4}(?=-\\d{2}\\.)")
  years_present <- sort(unique(file_years[!is.na(file_years)]))
  
  for (yr in years_present) {
    site_year_jobs[[length(site_year_jobs) + 1]] <- list(
      tower_id = tower_id, yr = yr,
      files = all_files[file_years == yr]
    )
  }
}

cat("Total site-year jobs to process:", length(site_year_jobs), "\n")

# ---- 4. Process every site-year, with a progress bar ---------------------
results <- tibble(
  tower_id = character(), neon_site = character(), year = character(),
  avg_daytime_extent_m = numeric(), n_daytime_intervals = integer(),
  n_files_skipped = integer(), status = character()
)

pb <- txtProgressBar(min = 0, max = length(site_year_jobs), style = 3)

for (j in seq_along(site_year_jobs)) {
  job <- site_year_jobs[[j]]
  tower_row <- towers_df %>% filter(Site.ID == job$tower_id)
  
  cat("\n==== [", j, "/", length(site_year_jobs), "] ", job$tower_id,
      " - ", job$yr, " ====\n", sep = "")
  
  row <- tryCatch(
    process_site_year(job$tower_id, tower_row$neon_site, tower_row$Lat,
                      tower_row$Lon, job$yr, job$files),
    error = function(e) {
      cat("  unexpected error:", conditionMessage(e), "\n")
      tibble(tower_id = job$tower_id, neon_site = tower_row$neon_site, year = job$yr,
             avg_daytime_extent_m = NA, n_daytime_intervals = 0,
             n_files_skipped = NA, status = paste("error:", conditionMessage(e)))
    }
  )
  
  cat("  ", row$status, "- extent:", round(row$avg_daytime_extent_m, 1), "m\n")
  results <- bind_rows(results, row)
  setTxtProgressBar(pb, j)
  
  # save incrementally so a crash/interruption overnight doesn't lose progress
  write.csv(results, out_csv, row.names = FALSE)
}

close(pb)
cat("\nDone. Annual summary written to:", out_csv, "\n")