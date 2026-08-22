# ============================================================================
# Diagnostic: raw L(c) curve inspection for SRER, to check whether the
# adaptive-FCM richness=2 floor result is a genuine local maximum of the
# validity function, or an artifact of c_min=2 truncating a curve that would
# keep decreasing if allowed to go lower (i.e. true structure is c=1,
# "no meaningful subdivision", not c=2).
#
# Must be run in the same session as AnnualSpectralDiversity.R and after
# FCM_test.R has been sourced (needs adaptive_fcm_ssr() from
# adaptive_fcm_ssr.py, plus site_year_jobs / get_tower_reflectance /
# compute_ndvi_raster / ndvi_thresh / n_subsample_pixels / n_pc_ssr already
# in the global environment).
# ============================================================================

library(reticulate)
# Import fcmeans directly here rather than relying on adaptive_fcm_ssr()'s
# wrapper, since that wrapper only returns the argmax -- we need L(c) itself
# for every candidate c, not just the winner.
fcmeans_mod <- import("fcmeans")
np <- import("numpy")

# ---- validity function L(c), directly from Wu et al. Eq. 4 ----
# L(c) = [ sum_i sum_j u_ij^w * ||v_i - xbar||^2 / (c-1) ]
#      / [ sum_i sum_j u_ij^w * ||x_j - v_i||^2 / (m-c) ]
compute_validity_L <- function(pcs, c, fuzziness, random_state) {
  m <- nrow(pcs)
  fcm <- fcmeans_mod$FCM(n_clusters = as.integer(c), m = fuzziness,
                         random_state = as.integer(random_state))
  fcm$fit(pcs)
  
  u <- fcm$u                     # membership matrix, m x c
  v <- fcm$centers                # c x n_pc
  x_bar <- colMeans(pcs)          # weighted centroid; using unweighted mean
  # here since FCM.fit doesn't expose sample
  # weights -- matches how adaptive_fcm_ssr.py
  # computes it (confirm this assumption
  # against that file if in doubt)
  
  u_w <- u ^ fuzziness
  
  # between-cluster term: sum_i sum_j u_ij^w * ||v_i - xbar||^2
  between <- 0
  for (i in seq_len(c)) {
    d2 <- sum((v[i, ] - x_bar)^2)
    between <- between + sum(u_w[, i]) * d2
  }
  between <- between / (c - 1)
  
  # within-cluster term: sum_i sum_j u_ij^w * ||x_j - v_i||^2
  within <- 0
  for (i in seq_len(c)) {
    d2_j <- rowSums(sweep(pcs, 2, v[i, ], "-")^2)
    within <- within + sum(u_w[, i] * d2_j)
  }
  within <- within / (m - c)
  
  list(L = between / within, between = between, within = within)
}

# ---- pull one SRER site-year, reuse the exact PCA/subsample recipe ----
srer_job <- keep(site_year_jobs, ~ .x$tower_id == "US-xSR" & .x$year == "2017")[[1]]
neon_site <- towers_df$neon_site[towers_df$Site.ID == srer_job$tower_id]

cat("==== Reading", srer_job$tower_id, "(", neon_site, ") 2017 for L(c) curve check ====\n")
data <- get_tower_reflectance(srer_job$tower_id, srer_job$files, buffer_m)
r <- data$raster; wl <- data$wavelengths
ndvi <- compute_ndvi_raster(r, wl)
veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)

r_masked <- mask(r, veg_mask)
vals <- values(r_masked, na.rm = TRUE)
pca <- prcomp(vals, center = TRUE, scale. = FALSE)
pcs_all <- pca$x[, 1:n_pc_ssr]

# ---- inspect L(c) across c = 2..20 for a handful of independent reps ----
n_check_reps <- 5
c_range <- 2:20

curve_results <- list()
for (rep_i in seq_len(n_check_reps)) {
  set.seed(90000 + rep_i)
  sub_idx <- sample(seq_len(nrow(pcs_all)), n_subsample_pixels)
  pcs_sub <- pcs_all[sub_idx, , drop = FALSE]
  
  cat("\n-- rep", rep_i, "--\n")
  L_vals <- numeric(length(c_range))
  for (ci in seq_along(c_range)) {
    c_val <- c_range[ci]
    res <- compute_validity_L(pcs_sub, c_val, fuzziness = 2.0,
                              random_state = 90000 + rep_i * 100 + c_val)
    L_vals[ci] <- res$L
    cat("   c =", sprintf("%2d", c_val), " L(c) =", format(res$L, digits = 6), "\n")
  }
  winning_c <- c_range[which.max(L_vals)]
  cat("   >> winning c for this rep:", winning_c,
      if (winning_c == min(c_range)) "  !! AT FLOOR c_min -- check if L(c) is still decreasing here" else "", "\n")
  
  curve_results[[rep_i]] <- tibble(rep = rep_i, c = c_range, L = L_vals)
}

curve_df <- bind_rows(curve_results)

# Quick check: is L(c) monotonically decreasing across the WHOLE range for
# most reps? If so, the floor is very likely truncating a curve that wants
# to go below c=2 (i.e. true answer might be "no meaningful clustering" /
# c=1), not settling on a genuine local peak at c=2.
monotone_check <- curve_df %>%
  group_by(rep) %>%
  summarise(is_monotone_decreasing = all(diff(L) < 0), .groups = "drop")
print(monotone_check)

cat("\nIf is_monotone_decreasing is TRUE for most/all reps, richness=2 is likely a\n",
    "floor artifact, not real structure -- c_min may need to allow 1, or SRER's\n",
    "true signal may be 'no meaningful spectral subdivision' rather than '2 groups'.\n",
    "If FALSE (L(c) rises then falls, peaking genuinely at c=2), the result stands.\n", sep = "")
#============================
#=============================
#==============================
# ============================================================================
# Diagnostic, part 2: (a) run the same L(c)-curve inspection on ABBY, to see
# whether its curve has a genuine interior peak (rises then falls) or is
# just "declining from the boundary, but at a higher plateau" like SRER's
# turned out to be; (b) extend c_min down to 1 for both sites, to check
# whether a genuinely homogeneous window (per the SRER result) actually
# prefers c=1 once it's allowed to.
#
# NOTE on c=1: Wu et al.'s L(c) = between/(c-1) / [within/(m-c)] is undefined
# at c=1 (between-cluster separation has no meaning with only one cluster --
# there's nothing to separate FROM). Rather than silently returning Inf/NaN,
# c=1 is scored separately here: within-cluster dispersion of the WHOLE
# population around its single global centroid (no separation term at all),
# and c=1 "wins" only if that within-1 dispersion is LOWER than the best
# within-cluster dispersion achieved by any c>=2 candidate -- i.e. c=1 wins
# when adding ANY split fails to meaningfully tighten cluster cohesion. This
# is a reasonable operationalization of "no meaningful subdivision" but is a
# genuine methodological choice, not something Wu et al. specify (they never
# needed to, since their c floor was implicitly >=2 in practice) -- flag this
# choice explicitly if used in anything beyond this diagnostic.
#
# Must be run in the same session as the original floor-check chunk (reuses
# fcmeans_mod, np, compute_validity_L(), site_year_jobs, get_tower_reflectance,
# compute_ndvi_raster, ndvi_thresh, n_subsample_pixels, n_pc_ssr).
# ============================================================================

# ---- within-cluster dispersion for c=1 (whole-population variance around
#      the global centroid) -- comparable in SCALE to the "within" component
#      of compute_validity_L()'s c>=2 output (i.e. NOT divided by (m-c) the
#      same way, since c=1 has no such term -- divide by m instead, the
#      natural c=1 analog) ----
compute_within_1 <- function(pcs) {
  m <- nrow(pcs)
  x_bar <- colMeans(pcs)
  within_1 <- sum(rowSums(sweep(pcs, 2, x_bar, "-")^2)) / m
  within_1
}

# ---- run the full c=1..20 check (c=1 scored via within-dispersion-only
#      comparison, c=2..20 via the standard L(c)) for one site, printing the
#      curve and flagging whether c=1 "wins" ----
run_floor_check <- function(tower_id, year, label, n_check_reps = 5, seed_offset = 0) {
  job <- keep(site_year_jobs, ~ .x$tower_id == tower_id & .x$year == year)[[1]]
  neon_site <- towers_df$neon_site[towers_df$Site.ID == job$tower_id]
  
  cat("\n\n======== ", label, ": ", tower_id, "(", neon_site, ")", year, " ========\n", sep = "")
  data <- get_tower_reflectance(job$tower_id, job$files, buffer_m)
  r <- data$raster; wl <- data$wavelengths
  ndvi <- compute_ndvi_raster(r, wl)
  veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)
  
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pcs_all <- pca$x[, 1:n_pc_ssr]
  
  c_range <- 2:20
  results <- list()
  
  for (rep_i in seq_len(n_check_reps)) {
    set.seed(seed_offset + rep_i)
    sub_idx <- sample(seq_len(nrow(pcs_all)), n_subsample_pixels)
    pcs_sub <- pcs_all[sub_idx, , drop = FALSE]
    
    cat("\n-- rep", rep_i, "--\n")
    within_1 <- compute_within_1(pcs_sub)
    cat("   c =  1  within-dispersion =", format(within_1, digits = 6),
        " (no L(c) -- see header note)\n")
    
    L_vals <- numeric(length(c_range))
    within_vals <- numeric(length(c_range))
    for (ci in seq_along(c_range)) {
      c_val <- c_range[ci]
      res <- compute_validity_L(pcs_sub, c_val, fuzziness = 2.0,
                                random_state = seed_offset + rep_i * 100 + c_val)
      L_vals[ci] <- res$L
      within_vals[ci] <- res$within
      cat("   c =", sprintf("%2d", c_val), " L(c) =", format(res$L, digits = 6),
          " (within =", format(res$within, digits = 6), ")\n")
    }
    
    # c=1 "wins" only if no c>=2 candidate achieves meaningfully tighter
    # within-cluster dispersion than the whole population's own baseline --
    # i.e. splitting into ANY number of groups fails to improve cohesion.
    best_within_ge2 <- min(within_vals)
    c1_wins <- within_1 <= best_within_ge2
    
    if (c1_wins) {
      winning_c <- 1
    } else {
      winning_c <- c_range[which.max(L_vals)]
    }
    cat("   >> winning c for this rep:", winning_c,
        if (winning_c == 1) "  (c=1 wins -- no split improved cohesion; homogeneous)"
        else if (winning_c == min(c_range)) "  !! at c=2 floor -- check curve shape"
        else "", "\n")
    
    results[[rep_i]] <- tibble(rep = rep_i, label = label, tower_id = tower_id,
                               c = c(1, c_range), L = c(NA, L_vals),
                               within = c(within_1, within_vals),
                               winning_c = winning_c)
  }
  
  bind_rows(results)
}

# ---- run both sites ----
abby_curve_df <- run_floor_check("US-xAB", "2017", "ABBY", n_check_reps = 5, seed_offset = 91000)
srer_curve_df <- run_floor_check("US-xSR", "2017", "SRER", n_check_reps = 5, seed_offset = 92000)

curve_df_both <- bind_rows(abby_curve_df, srer_curve_df)

# ---- summary: winning c per rep, per site, with the c=1-allowed floor ----
cat("\n\n==== Summary: winning c (c_min=1 allowed) ====\n")
curve_df_both %>%
  distinct(label, rep, winning_c) %>%
  group_by(label) %>%
  summarise(winning_c_values = paste(winning_c, collapse = ", "),
            mean_winning_c = mean(winning_c), .groups = "drop") %>%
  print()

# ---- shape check: for c>=2 only, does ABBY show a genuine interior peak
#      (L(c) rises then falls) vs. SRER's already-confirmed monotone-decline-
#      from-c=2 shape? ----
cat("\n==== Shape check (c=2..20 only): does L(c) have an interior peak? ====\n")
curve_df_both %>%
  filter(c >= 2) %>%
  group_by(label, rep) %>%
  summarise(peak_c = c[which.max(L)],
            peak_is_at_boundary = peak_c %in% c(min(c), max(c)),
            .groups = "drop") %>%
  print(n = Inf)

cat("\nIf peak_is_at_boundary is TRUE for most SRER reps but FALSE for most ABBY\n",
    "reps, that confirms ABBY has genuine interior structure while SRER's c>=2\n",
    "result was purely a floor artifact -- consistent with c=1 now winning for SRER.\n", sep = "")