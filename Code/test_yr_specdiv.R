# ---- Test: process only US-xAB across all its discovered years ------------
test_jobs <- keep(site_year_jobs, ~ .x$tower_id == "US-xAB")
cat("Testing", length(test_jobs), "site-year job(s) for US-xAB:\n")
walk(test_jobs, ~ cat("  ", .x$tower_id, .x$year, "-", length(.x$files), "files\n"))

results <- tibble(
  tower_id = character(), neon_site = character(), year = character(),
  cv = double(), chv = double(),
  spectral_species_richness = double(),
  raoq_ndvi = double(), raoq_nirv = double(), raoq_allbands = double(),
  status = character()
)

for (job in test_jobs) {
  
  t0 <- Sys.time()
  cat("\n==== ", job$tower_id, " - ", job$year, " ====\n", sep = "")
  
  neon_site <- towers_df$neon_site[towers_df$Site.ID == job$tower_id]
  
  data <- tryCatch(
    get_tower_reflectance(job$tower_id, job$files, buffer_m),
    error = function(e) { cat("  read/mosaic/crop failed:", conditionMessage(e), "\n"); NULL }
  )
  
  if (is.null(data)) {
    results <- add_row(results, tower_id = job$tower_id, neon_site = neon_site, year = job$year,
                       cv = NA, chv = NA, spectral_species_richness = NA,
                       raoq_ndvi = NA, raoq_nirv = NA, raoq_allbands = NA,
                       status = "read failed")
    next
  }
  
  r <- data$raster; wl <- data$wavelengths
  cat("  raster dims:", nrow(r), "x", ncol(r), "x", nlyr(r), "layers\n")
  cat("  read+mosaic+crop:", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")
  
  ndvi <- compute_ndvi_raster(r, wl)
  nirv <- compute_nirv_raster(r, wl)
  veg_mask <- ndvi > ndvi_thresh
  cat("  veg pixels:", sum(values(veg_mask), na.rm = TRUE), "/",
      sum(!is.na(values(ndvi))), "\n")
  
  t1 <- Sys.time()
  cv_val <- compute_cv(r, veg_mask)
  cat("  CV:", round(cv_val, 4), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  t1 <- Sys.time()
  chv_val <- compute_chv(r, veg_mask, n_pc_chv)
  cat("  CHV:", round(chv_val, 4), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  t1 <- Sys.time()
  ssr_val <- compute_spectral_species_richness(r, veg_mask, n_clusters,
                                               n_subsample_pixels, n_pc_ssr, n_reps_ssr)
  cat("  SSR:", round(ssr_val, 2), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  t1 <- Sys.time()
  raoq_ndvi <- compute_rao_q_raster_rasterdiv(mask(ndvi, veg_mask), raoq_window)
  cat("  RaoQ (NDVI):", round(raoq_ndvi, 5), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  t1 <- Sys.time()
  raoq_nirv <- compute_rao_q_raster_rasterdiv(mask(nirv, veg_mask), raoq_window)
  cat("  RaoQ (NIRv):", round(raoq_nirv, 5), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  t1 <- Sys.time()
  raoq_all <- compute_rao_q_raster_rasterdiv(mask(r, veg_mask), raoq_window)
  cat("  RaoQ (all-bands):", round(raoq_all, 5), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")
  
  results <- add_row(results, tower_id = job$tower_id, neon_site = neon_site, year = job$year,
                     cv = cv_val, chv = chv_val, spectral_species_richness = ssr_val,
                     raoq_ndvi = raoq_ndvi, raoq_nirv = raoq_nirv, raoq_allbands = raoq_all,
                     status = "success")
  
  cat("  TOTAL for this site-year:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
}

cat("\n==== US-xAB test results ====\n")
print(results, n = Inf)