# ---- Diagnostic: n_clusters sensitivity sweep, ABBY vs SRER, multiple years ----
# Tests candidate (a) from the SSR open-problem note: does lowering n_clusters
# restore separation between a dense-forest site and a sparse-desert site, or
# does the near-ceiling behavior (SSR saturation, Shannon's H' near ln(k))
# persist even at smaller k? Fixed k WITHIN each comparison (never adaptive
# per site-year) -- see Feret & Asner 2014 on why k must stay constant across
# sites being compared. This is diagnostic only; does not change the
# production pipeline or n_clusters config value.

k_values   <- c(10, 20, 30, 40, 50)
test_sites <- list(
  list(tower_id = "US-xAB", years = c("2017", "2018", "2019")),  # add more years as available
  list(tower_id = "US-xSR", years = c("2017", "2018", "2019"))
)

sweep_results <- list()

for (site in test_sites) {
  for (yr in site$years) {
    job <- keep(site_year_jobs, ~ .x$tower_id == site$tower_id & .x$year == yr)
    if (length(job) == 0) {
      cat("  SKIP:", site$tower_id, yr, "-- not found in site_year_jobs\n")
      next
    }
    job <- job[[1]]
    neon_site <- towers_df$neon_site[towers_df$Site.ID == job$tower_id]
    
    cat("\n==== Reading", job$tower_id, "(", neon_site, ")", yr, "====\n")
    t0 <- Sys.time()
    data <- tryCatch(
      get_tower_reflectance(job$tower_id, job$files, buffer_m),
      error = function(e) { cat("  read failed:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(data)) next
    r <- data$raster; wl <- data$wavelengths
    cat("  read+mosaic+crop:", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")
    
    ndvi <- compute_ndvi_raster(r, wl)
    veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)
    n_veg <- sum(!is.na(values(veg_mask)))
    cat("  veg pixels:", n_veg, "\n")
    
    # Read the raster ONCE per site-year, reuse across all k values in the
    # sweep -- this is the expensive step (crop-then-mosaic), no reason to
    # repeat it per k.
    for (k in k_values) {
      cat("  -- k =", k, "--\n")
      t1 <- Sys.time()
      ssr_result <- compute_spectral_species_richness(
        r, veg_mask,
        n_clusters   = k,
        n_subsample  = n_subsample_pixels,
        n_pc         = n_pc_ssr,
        n_reps       = n_reps_ssr
      )
      elapsed <- round(difftime(Sys.time(), t1, units = "secs"), 1)
      
      cat("     richness:", round(ssr_result$richness, 2),
          " (", round(100 * ssr_result$richness / k, 1), "% of k)",
          " | H':", round(ssr_result$shannon_h, 4),
          " (", round(100 * ssr_result$shannon_h / log(k), 1), "% of ln(k) ceiling)",
          " | exp(H'):", round(ssr_result$shannon_effective, 2),
          " (", elapsed, "sec )\n", sep = "")
      
      sweep_results[[length(sweep_results) + 1]] <- tibble(
        tower_id = job$tower_id, neon_site = neon_site, year = yr,
        n_clusters = k, veg_pixels = n_veg,
        richness = ssr_result$richness,
        richness_pct_of_k = 100 * ssr_result$richness / k,
        shannon_h = ssr_result$shannon_h,
        shannon_h_pct_of_ceiling = 100 * ssr_result$shannon_h / log(k),
        shannon_effective = ssr_result$shannon_effective,
        elapsed_sec = elapsed
      )
    }
  }
}

sweep_df <- bind_rows(sweep_results)
cat("\n\n==== Full sweep results ====\n")
print(sweep_df, n = Inf)

# Quick visual check: does the gap between sites widen at lower k?
cat("\n==== ABBY vs SRER gap by k (2017) ====\n")
sweep_df %>%
  filter(year == "2017") %>%
  select(tower_id, n_clusters, richness, shannon_h, shannon_effective) %>%
  pivot_wider(names_from = tower_id,
              values_from = c(richness, shannon_h, shannon_effective)) %>%
  print()