# ---- Diagnostic: self-adaptive FCM SSR, ABBY vs SRER, multiple years ----
# Companion to KMeans_test.R (same repo, same directory). KMeans_test.R sweeps
# n_clusters in {10,20,30,40,50} through the EXISTING compute_spectral_species_
# richness() (RF-proximity -> k-means -> nearest-centroid) to check whether
# lowering n_clusters restores ABBY/SRER separation or whether SSR saturation
# persists regardless of k. This script asks the same question with a
# DIFFERENT clustering method: self-adaptive fuzzy c-means (FCM, Wu et al.
# 2026, Computers and Electronics in Agriculture 252:112108), where the
# cluster count c is chosen PER REP by maximizing a validity function L(c)
# (Section 1 below) over c in [2, 20] -- capped well below the old
# n_clusters=50 ceiling that caused the saturation this whole diagnostic
# lineage (KMeans_test.R included) is chasing.
#
# LIKE KMeans_test.R, this script is NOT standalone -- it assumes
# AnnualSpectralDiversity.R has already been run/sourced in the same R
# session, so site_year_jobs, towers_df, get_tower_reflectance(),
# compute_ndvi_raster(), buffer_m, ndvi_thresh, n_subsample_pixels, n_pc_ssr,
# and n_reps_ssr already exist in the global environment. Diagnostic only --
# does not change the production pipeline, and per instruction does NOT
# modify KMeans_test.R, AnnualSpectralDiversity.R, or ComputeSpecBiodiv.R
# (confirmed untouched -- this script only READS the objects those define).
#
# Output table shape matches KMeans_test.R's sweep_df exactly in spirit (one
# row per site/year), so the two can be joined/compared directly: tower_id,
# neon_site, year, veg_pixels, richness, richness_pct_of_cap (this script's
# analog of KMeans_test.R's richness_pct_of_k, using the FIXED c_max=20 cap
# in place of a swept k), mean_winning_c, max_winning_c (cap-bunching
# check), shannon_h, shannon_effective, elapsed_sec.
#
# ============================================================================
# DESIGN NOTES / ASSUMPTIONS -- flagged per instruction, not hard-coded gospel:
#
# 1. FUZZINESS w = 2.0 (Section 0) -- standard/default FCM fuzziness exponent
#    in the literature. Not tuned or re-derived for this data.
#
# 2. c-RANGE capped at [2, 20] (Section 0), per explicit instruction --
#    deliberately far below the old n_clusters=50 ceiling under diagnosis.
#    mean_winning_c / max_winning_c are logged per site-year specifically so
#    a winning c that bunches against 20 (this script's cap) can be spotted,
#    the same way KMeans_test.R's richness_pct_of_k lets you see richness
#    tracking k rather than true structure.
#
# 3. RICHNESS DEFINITION differs structurally from the k-means version, by
#    the task's own explicit design, not an oversight: KMeans_test.R's
#    compute_spectral_species_richness() clusters a SUBSAMPLE, then
#    reassigns the FULL vegetation-masked population to the nearest
#    resulting centroid (assign_nearest_centroid(), AnnualSpectralDiversity.R
#    Section 5) and counts how many of the n_clusters labels the full
#    population actually uses -- richness can come out BELOW n_clusters.
#    compute_spectral_species_richness_fcm() below does NOT do a full-
#    population reassignment step at all -- each rep's richness is simply
#    the winning c (argmax_c L(c)) from that rep's SUBSAMPLE-only search, as
#    specified. There is no equivalent to k-means's "count of labels actually
#    used after reassigning everyone" for FCM here.
#
# 4. SHANNON'S H' -- SAME FORMULA as the k-means version
#    (-sum(p_i * log(p_i))), but NOT the same POPULATION SCOPE. The k-means
#    version's H' comes from proportions over the FULL reassigned population
#    (same cluster_all as its richness). This script's H' comes from
#    proportions over the WINNING c's HARD assignment (argmax of fuzzy
#    membership, via adaptive_fcm_ssr.py's `hard_assignment`) of that rep's
#    SUBSAMPLE only -- there is no full-population step to draw from (see
#    note 3). Hard assignment itself was also an explicit instruction (for
#    formula comparability with k-means's hard cluster labels), not a claim
#    that hard assignment is the "correct" way to summarize FCM's fuzzy
#    output -- a soft/fuzzy diversity index (e.g. weighted by membership
#    rather than argmax) would be a legitimate alternative for a future
#    session. Net effect: shannon_h / shannon_effective in this script's
#    output and in KMeans_test.R's sweep_df are comparable in FORMULA only,
#    not in what population they were computed over -- do not treat a
#    side-by-side H' comparison as apples-to-apples without accounting for
#    this.
#
# 5. n_subsample = n_subsample_pixels (2500), n_reps = n_reps_ssr (20) --
#    SAME as the k-means sweep, per explicit instruction, even though this
#    makes FCM far more expensive per rep (~19 candidate-c fits vs. 1 k-means
#    fit). Not reduced for speed -- see Section 2's runtime estimate, printed
#    BEFORE the main loop runs, for the actual projected cost at these
#    settings.
#
# 6. FCM fitting via Python's `fcmeans` package (PyPI `fuzzy-c-means`) through
#    reticulate::source_python() on adaptive_fcm_ssr.py -- the SAME module
#    CompareSSR_AdaptiveFCM_vs_KMeans.R uses (not duplicated here), because
#    Wu et al. used this exact package and reproducing their method means
#    reusing the same fitting code, not a different FCM implementation that
#    could converge to different partitions. Section 1 below checks
#    fcmeans/numpy are actually importable BEFORE the site-year loop starts,
#    and stop()s with a clear message (not a deep-loop failure) if not.
#
# 7. NAME COLLISION WARNING (inherited from adaptive_fcm_ssr.py, restated
#    here since it's easy to miss): the paper's `w` (fuzzy weighting
#    exponent) is passed to `fcmeans.FCM` as its `m` argument; the paper's
#    `m` (pixel count) is unrelated. This script always writes `fuzziness`
#    for the former and lets `nrow(pcs_sub)` speak for the latter -- never a
#    bare `m`.
#
# SELF-TEST: done EXTERNALLY before delivery (not an embedded section in
# this file, matching KMeans_test.R's own structure, which has no embedded
# self-test either) -- a synthetic PC-space matrix with 5 well-separated
# Gaussian blobs and a separate homogeneous single-cluster matrix, both run
# through compute_spectral_species_richness_fcm()'s inner per-rep call
# (adaptive_fcm_ssr()) directly (no raster/H5 needed), confirming winning c
# behaves sensibly (low for homogeneous, correctly recovers ~5 for the
# blobs). Also smoke-tested the full site-year loop shape against a mocked
# get_tower_reflectance() returning synthetic rasters. Neither is a
# real-data validation -- only the server run (real site_year_jobs, real H5
# tiles) can confirm this against ABBY/SRER themselves.

library(reticulate)
# tidyr is needed for pivot_wider() (Section 4's final ABBY-vs-SRER table).
# Neither KMeans_test.R nor AnnualSpectralDiversity.R load tidyr anywhere,
# even though KMeans_test.R also calls pivot_wider() at its own end -- that
# is a pre-existing latent gap in this diagnostic-script pair (both only
# work if the interactive session has tidyr loaded from some other, unseen
# step), found while self-testing THIS script in a clean session. Loading
# it explicitly here so FCM_test.R works regardless of what else is loaded;
# worth checking whether KMeans_test.R's run on the server hit the same gap.
library(tidyr)

# ============================================================================
# 0. Config
# ============================================================================
c_min           <- 2
c_max           <- 20     # deliberately << the old n_clusters=50 -- see note 2
fcm_fuzziness_w <- 2.0    # ASSUMPTION -- see note 1
fcm_seed_base   <- 20260101   # arbitrary fixed base, only for run-to-run reproducibility

adaptive_fcm_py <- "./Code/DataAnalysis/adaptive_fcm_ssr.py"

test_sites <- list(
  list(tower_id = "US-xAB", years = c("2017", "2018", "2019")),
  list(tower_id = "US-xSR", years = c("2017", "2018", "2019"))
)

# Guard against being run standalone (outside an AnnualSpectralDiversity.R
# session) -- KMeans_test.R has no equivalent guard and relies on convention
# alone; this one is a strict addition, not a structural deviation, and
# fails clearly instead of erroring confusingly deep in the loop on a
# missing `site_year_jobs`.
required_objects <- c("site_year_jobs", "towers_df", "get_tower_reflectance",
                       "compute_ndvi_raster", "buffer_m", "ndvi_thresh",
                       "n_subsample_pixels", "n_pc_ssr", "n_reps_ssr")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("FCM_test.R must be run in the SAME R session as AnnualSpectralDiversity.R ",
       "(same convention as KMeans_test.R) -- missing: ",
       paste(missing_objects, collapse = ", "))
}
if (!file.exists(adaptive_fcm_py)) {
  stop("Required input not found: ", adaptive_fcm_py)
}

# ============================================================================
# 1. FCM availability check -- fail HERE, clearly, not deep in the site-year
#    loop (per instruction)
# ============================================================================
cat("==== Checking fcmeans/numpy availability via reticulate ====\n")
if (!py_module_available("numpy")) {
  stop("Python module 'numpy' is not importable via reticulate. Check ",
       "RETICULATE_PYTHON / reticulate::py_config() before running the ",
       "real loop -- do not proceed on a guess.")
}
if (!py_module_available("fcmeans")) {
  stop("Python module 'fcmeans' (PyPI package 'fuzzy-c-means') is not ",
       "importable via reticulate. Install it in whatever Python ",
       "reticulate is pointed at (e.g. `pip install fuzzy-c-means`), or ",
       "point RETICULATE_PYTHON at an environment that already has it, ",
       "before running the real loop -- do not proceed on a guess.")
}
source_python(adaptive_fcm_py)   # exposes adaptive_fcm_ssr()
py_fcm_version <- tryCatch(
  as.character(import("importlib.metadata")$version("fuzzy-c-means")),
  error = function(e) "unknown (could not read package metadata)"
)
cat("  numpy: available\n")
cat("  fcmeans (fuzzy-c-means): available, version", py_fcm_version, "\n\n")

# ============================================================================
# 2. compute_spectral_species_richness_fcm() -- mirrors
#    compute_spectral_species_richness()'s masking/PCA/subsample-rep-loop
#    SHAPE exactly (AnnualSpectralDiversity.R, Section 5); the clustering
#    step itself is FCM + validity-function search instead of RF-proximity +
#    k-means + nearest-centroid. See header notes 3-4 for how richness and
#    Shannon's H' deliberately differ in what they measure, not just in
#    algorithm.
# ============================================================================
compute_spectral_species_richness_fcm <- function(r, veg_mask, c_min = 2, c_max = 20,
                                                    fuzziness = 2.0, n_subsample = 2500,
                                                    n_pc = 4, n_reps = 20, seed_base = 0) {
  # SAME masking + PCA recipe as compute_spectral_species_richness()
  # (AnnualSpectralDiversity.R, Section 5) -- identical PC-space inputs for
  # both the k-means sweep and this FCM sweep.
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  this_n_subsample <- n_subsample
  if (nrow(vals) < this_n_subsample) this_n_subsample <- nrow(vals)
  if (nrow(vals) <= c_max) {
    return(list(richness = NA_real_, shannon_h = NA_real_, shannon_effective = NA_real_,
                winning_c_per_rep = integer(0)))
  }

  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pcs_all <- pca$x[, 1:n_pc]

  richness_reps  <- numeric(n_reps)
  shannon_reps   <- numeric(n_reps)
  winning_c_reps <- integer(n_reps)
  for (rep_i in seq_len(n_reps)) {
    sub_idx <- sample(seq_len(nrow(pcs_all)), this_n_subsample)
    pcs_sub <- pcs_all[sub_idx, , drop = FALSE]

    res <- adaptive_fcm_ssr(pcs_sub, c_min = as.integer(c_min), c_max = as.integer(c_max),
                             fuzziness = fuzziness,
                             random_state = as.integer(seed_base + rep_i))

    winning_c_reps[rep_i] <- res$winning_c
    # per instruction: argmax_c L(c) IS this rep's richness value -- no
    # full-population reassignment step (see header note 3)
    richness_reps[rep_i]  <- res$winning_c

    # Shannon's H', SAME FORMULA as the k-means version, over the winning
    # c's HARD assignment of the SUBSAMPLE (see header note 4 for the
    # population-scope caveat). Python's hard_assignment is 0-indexed;
    # +1L just to look like R factor labels, doesn't affect table()/p_i.
    hard <- unlist(res$hard_assignment) + 1L
    tbl <- table(hard)
    p_i <- as.numeric(tbl[tbl > 0]) / length(hard)
    shannon_reps[rep_i] <- -sum(p_i * log(p_i))
  }

  list(richness = mean(richness_reps),
       shannon_h = mean(shannon_reps),
       shannon_effective = mean(exp(shannon_reps)),
       winning_c_per_rep = winning_c_reps)
}

# ============================================================================
# 3. Runtime estimate, BEFORE the real loop -- calibration probe against
#    synthetic data shaped like the real subsample (n_subsample_pixels x
#    n_pc_ssr), at the SAME c_min/c_max/fuzziness/n_reps this script will
#    actually use. This is a rough guide, not a guarantee: real reflectance
#    PC data may need more/fewer FCM iterations to converge than random
#    synthetic data, so treat this as a floor, not an exact prediction.
# ============================================================================
cat("==== Runtime estimate (calibration probe on synthetic-shaped data, before the real loop) ====\n")
set.seed(fcm_seed_base)
calib_pcs <- matrix(rnorm(n_subsample_pixels * n_pc_ssr), ncol = n_pc_ssr)
t_calib <- Sys.time()
invisible(adaptive_fcm_ssr(calib_pcs, c_min = as.integer(c_min), c_max = as.integer(c_max),
                            fuzziness = fcm_fuzziness_w, random_state = 0L))
sec_per_rep <- as.numeric(difftime(Sys.time(), t_calib, units = "secs"))

n_target_site_years <- sum(vapply(test_sites, function(s) length(s$years), integer(1)))
sec_per_site_year <- sec_per_rep * n_reps_ssr
sec_total <- sec_per_site_year * n_target_site_years

cat("  One rep (", c_max - c_min + 1, " candidate-c fits, n=", n_subsample_pixels,
    " px, ", n_pc_ssr, " PCs): ", round(sec_calib <- sec_per_rep, 1), " sec\n", sep = "")
cat("  Projected per site-year (", n_reps_ssr, " reps): ",
    round(sec_per_site_year / 60, 1), " min\n", sep = "")
cat("  Projected TOTAL across ", n_target_site_years, " target site-years: ",
    round(sec_total / 60, 1), " min (", round(sec_total / 3600, 2), " hr)\n", sep = "")
cat("  (Calibration only -- real reflectance PC data may converge slower/faster",
    " than this synthetic probe; not automatically reducing n_reps or the",
    " c-range to hit a time budget, per instruction.)\n\n", sep = "")

# ============================================================================
# 4. Main loop -- mirrors KMeans_test.R's site/year loop structure exactly
# ============================================================================
sweep_results_fcm <- list()

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

    cat("  -- adaptive FCM (c in [", c_min, ",", c_max, "], w = ", fcm_fuzziness_w,
        ", ", n_reps_ssr, " reps) --\n", sep = "")
    t1 <- Sys.time()
    job_seed_base <- fcm_seed_base + length(sweep_results_fcm) * 1000L
    fcm_result <- compute_spectral_species_richness_fcm(
      r, veg_mask, c_min = c_min, c_max = c_max, fuzziness = fcm_fuzziness_w,
      n_subsample = n_subsample_pixels, n_pc = n_pc_ssr, n_reps = n_reps_ssr,
      seed_base = job_seed_base
    )
    elapsed <- round(difftime(Sys.time(), t1, units = "secs"), 1)

    mean_winning_c <- mean(fcm_result$winning_c_per_rep)
    max_winning_c  <- max(fcm_result$winning_c_per_rep)

    cat("     richness:", round(fcm_result$richness, 2),
        " (", round(100 * fcm_result$richness / c_max, 1), "% of cap)",
        " | mean/max winning c:", round(mean_winning_c, 2), "/", max_winning_c,
        " | H':", round(fcm_result$shannon_h, 4),
        " | exp(H'):", round(fcm_result$shannon_effective, 2),
        " (", elapsed, "sec )\n", sep = "")
    if (max_winning_c == c_max) {
      cat("     !! at least one rep's winning c hit the c_max =", c_max,
          "cap -- possible sign the cap needs raising for this site-year.\n")
    }

    sweep_results_fcm[[length(sweep_results_fcm) + 1]] <- tibble(
      tower_id = job$tower_id, neon_site = neon_site, year = yr,
      veg_pixels = n_veg,
      richness = fcm_result$richness,
      richness_pct_of_cap = 100 * fcm_result$richness / c_max,
      mean_winning_c = mean_winning_c,
      max_winning_c = max_winning_c,
      shannon_h = fcm_result$shannon_h,
      shannon_effective = fcm_result$shannon_effective,
      elapsed_sec = as.numeric(elapsed)
    )
  }
}

sweep_df_fcm <- bind_rows(sweep_results_fcm)
cat("\n\n==== Full FCM sweep results ====\n")
print(sweep_df_fcm, n = Inf)

# ---- ABBY vs SRER gap, 2017 -- same shape as KMeans_test.R's final block ----
cat("\n==== ABBY vs SRER gap, adaptive-FCM SSR (2017) ====\n")
sweep_df_fcm %>%
  filter(year == "2017") %>%
  select(tower_id, richness, mean_winning_c, max_winning_c, shannon_h, shannon_effective) %>%
  pivot_wider(names_from = tower_id,
              values_from = c(richness, mean_winning_c, max_winning_c, shannon_h, shannon_effective)) %>%
  print()

cat("\nCompare this table against KMeans_test.R's own 2017 ABBY-vs-SRER print (same site-years,",
    " same PC-space) once both scripts have finished -- fixed-k SSR showed ABBY ~ SRER (no",
    " separation) in prior diagnostics; the question this script exists to answer is whether",
    " adaptive-FCM richness/mean_winning_c shows real separation instead.\n", sep = "")
