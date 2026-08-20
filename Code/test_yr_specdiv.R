# ---- Test: US-xAB 2017, using the CURRENT (fixed) pipeline functions ------
# Confirms, on real data: BOM-fixed towers_df, ifel()-based veg_mask,
# compute_rao_q_raster_local() (not rasterdiv), PCA-reduced all-bands Rao's Q,
# and the new list(richness, shannon_h, shannon_effective) SSR return value.

job <- keep(site_year_jobs, ~ .x$tower_id == "US-xSR" & .x$year == "2017")[[1]]
cat("Testing", job$tower_id, job$year, "-", length(job$files), "files\n")

t0 <- Sys.time()
cat("\n====", job$tower_id, "-", job$year, "====\n")

neon_site <- towers_df$neon_site[towers_df$Site.ID == job$tower_id]

data <- tryCatch(
  get_tower_reflectance(job$tower_id, job$files, buffer_m),
  error = function(e) { cat("  read/mosaic/crop failed:", conditionMessage(e), "\n"); NULL }
)
stopifnot(!is.null(data))

r <- data$raster; wl <- data$wavelengths
cat("  raster dims:", nrow(r), "x", ncol(r), "x", nlyr(r), "layers\n")
cat("  read+mosaic+crop:", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")

ndvi <- compute_ndvi_raster(r, wl)
nirv <- compute_nirv_raster(r, wl)

# ifel()-based mask -- confirm it's actually excluding cells this time
veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)
n_veg   <- sum(!is.na(values(veg_mask)))
n_total <- sum(!is.na(values(ndvi)))
cat("  veg pixels (mask):", n_veg, "/", n_total,
    " (", round(100 * n_veg / n_total, 1), "% )\n", sep = "")

t1 <- Sys.time()
cv_val <- compute_cv(r, veg_mask)
cat("  CV:", round(cv_val, 4), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

t1 <- Sys.time()
chv_val <- compute_chv(r, veg_mask, n_pc_chv)
cat("  CHV:", round(chv_val, 4), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

t1 <- Sys.time()
cha_val <- compute_cha(r, veg_mask)
cat("  CHA:", round(cha_val, 4), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

t1 <- Sys.time()
ssr_result <- compute_spectral_species_richness(r, veg_mask, n_clusters,
                                                n_subsample_pixels, n_pc_ssr, n_reps_ssr)
cat("  SSR:", round(ssr_result$richness, 2),
    " | Shannon H':", round(ssr_result$shannon_h, 4),
    " | exp(H'):", round(ssr_result$shannon_effective, 2),
    " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n", sep = "")

t1 <- Sys.time()
raoq_ndvi <- compute_rao_q_raster_local(mask(ndvi, veg_mask), raoq_window)
cat("  RaoQ (NDVI):", round(raoq_ndvi, 5), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

t1 <- Sys.time()
raoq_nirv <- compute_rao_q_raster_local(mask(nirv, veg_mask), raoq_window)
cat("  RaoQ (NIRv):", round(raoq_nirv, 5), " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

t1 <- Sys.time()
r_allbands_pca <- pca_reduce_raster(r, veg_mask, raoq_pca_var_threshold)
cat("  PCA reduction: all-bands ->", nlyr(r_allbands_pca), "PC(s) at",
    raoq_pca_var_threshold * 100, "% variance (",
    round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n", sep = "")

t1 <- Sys.time()
raoq_all <- compute_rao_q_raster_local(r_allbands_pca, raoq_window)
cat("  RaoQ (all-bands, PCA-reduced):", round(raoq_all, 5),
    " (", round(difftime(Sys.time(), t1, units = "secs"), 1), "sec )\n")

cat("\n  TOTAL for this site-year:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

test_result <- tibble(
  tower_id = job$tower_id, neon_site = neon_site, year = job$year,
  cv = cv_val, chv = chv_val, cha = cha_val,
  spectral_species_richness = ssr_result$richness,
  shannon_h = ssr_result$shannon_h, shannon_effective = ssr_result$shannon_effective,
  raoq_ndvi = raoq_ndvi, raoq_nirv = raoq_nirv, raoq_allbands = raoq_all,
  veg_pixel_pct = round(100 * n_veg / n_total, 1)
)
print(test_result)