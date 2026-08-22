# CompareSSR_AdaptiveFCM_vs_KMeans.R
#
# Side-by-side comparison of two DIFFERENT definitions of Spectral Species
# Richness (SSR) on real ABBY/SRER data:
#   1. EXISTING: fixed-n_clusters=50 RF-proximity + k-means + nearest-centroid
#      SSR, exactly as implemented in AnnualSpectralDiversity.R. Confirmed
#      (SESSION_LOG.md, synthetic validation) to saturate near the n_clusters
#      ceiling regardless of true underlying diversity -- a structural
#      property of forcing exactly n_clusters non-empty partitions on a
#      large, dense, continuous pixel cloud, NOT a bug in that
#      implementation. This script does not touch or recompute that logic
#      except to call it, unchanged, as the baseline.
#   2. NEW: self-adaptive fuzzy c-means (FCM) SSR (Wu et al. 2026, Computers
#      and Electronics in Agriculture 252:112108) -- cluster count c is
#      chosen PER WINDOW by maximizing a validity function L(c) that trades
#      fuzzy between-cluster separation against fuzzy within-cluster
#      dispersion (implemented in adaptive_fcm_ssr.py, Section 2 below).
#      This is a genuinely different metric definition, not a drop-in fix
#      for (1).
#
# Goal: does the adaptive-c approach resolve the saturation problem (i.e.
# separate ABBY, dense forest, from SRER, sparse desert, the way CV/CHA/
# Rao's Q already do) without excessive runtime cost? Not decided by this
# script -- see NON-GOALS below.
#
# ============================================================================
# CORRECTIONS to the task prompt this script was written against -- verified
# against the actual repo contents before writing anything, not assumed:
#
# 1. "Reuse FieldDiversity.R functions via import_functions_from() ... for
#    the existing pixel extraction pipeline." FALSE as stated.
#    FieldDiversity.R has NO pixel-extraction logic at all -- it is
#    ground-plot floristic diversity from percent-cover CSVs, and its own
#    header says explicitly "this script never reads reflectance values,
#    only H5 metadata." Its ONLY export relevant to any hyperspectral
#    pipeline is get_flight_acquisition_date() (a metadata-text-parsing
#    helper, used elsewhere for phenology date-matching, not extraction).
#    The real PC-space pixel extraction pipeline -- read_neon_h5_tile,
#    get_tower_reflectance, compute_ndvi_raster, discover_site_years, and
#    the existing nearest-centroid compute_spectral_species_richness()
#    itself -- lives in AnnualSpectralDiversity.R. That is what Section 1
#    below imports, via the SAME import_functions_from() bootstrap
#    ComputeLUE.R / ComputeLUE_Annual.R already established (parse() the
#    origin script, eval() only the named top-level `<-` bindings into a
#    dedicated env -- safe against AnnualSpectralDiversity.R's monolithic
#    top level, which would otherwise run its whole site-year pipeline the
#    moment it was source()'d).
#
# 2. "Use Python's fcmeans package ... the existing repo SSR step is
#    already implemented in Python for this same reason [[that Wu et al.
#    used fcmeans]]; reuse import_functions_from()-equivalent plumbing to
#    bridge back to R ... matching the existing repo pattern for CV/SV/CHA
#    in R + SSR in Python." FALSE as stated. SSR is NOT implemented in
#    Python anywhere in this repo -- compute_spectral_species_richness()
#    in AnnualSpectralDiversity.R is pure R (randomForest + kmeans +
#    nearest-centroid). The only prior Python bridging in this repo is
#    Rao's Q via pyGNDiv/reticulate, and only in the OLDER
#    ComputeSpecBiodiv.R -- AnnualSpectralDiversity.R deliberately DROPPED
#    that bridge for Rao's Q (its own Section 6 header: paRao's global
#    distance matrix didn't scale to ~1m continuous reflectance data) in
#    favor of a pure-R local moving-window implementation. So there is no
#    existing "SSR in Python" pattern to match, and no pre-existing
#    R<->Python plumbing for SSR to reuse.
#    Still, `fcmeans` (the PyPI `fuzzy-c-means` package) IS what Wu et al.
#    themselves used, and reproducing their method means using the same
#    fitting code, not a different FCM implementation that could converge
#    to different partitions. So this script DOES bridge to Python for the
#    FCM step, via reticulate::source_python() -- the same bridging
#    MECHANISM ComputeSpecBiodiv.R used for pyGNDiv, applied fresh here.
#    See adaptive_fcm_ssr.py for the Python side; Section 2 below for the
#    R-side wrapper.
#
# 3. Output granularity ("one row per plot/window combination ... plot ID,
#    plot size / window ID"). The existing nearest-centroid SSR pipeline
#    (AnnualSpectralDiversity.R) computes exactly ONE richness value per
#    tower-YEAR, from the entire vegetation-masked buffer_m buffer as a
#    single population -- it has no sub-buffer plot/window concept to be
#    directly comparable against. Introducing a NEW finer-than-site-year
#    windowing scheme (some window size/stride) would be a methodological
#    decision with no established convention in this repo to inherit, and
#    would make the two SSR values NOT comparable to the existing CSV (a
#    different spatial support), contradicting the stated goal of direct
#    comparability. This script therefore treats the whole per-tower-year
#    vegetation-masked buffer as the single "window" for BOTH methods
#    (window_id is always "whole_buffer"), matching the existing output's
#    granularity exactly. If per-window (sub-buffer) SSR is actually wanted
#    -- which would be a closer match to how Wu et al. describe their own
#    plot-scale use case -- that needs a window size/stride decided
#    explicitly first; flagged here as an open decision, not guessed at.
#
# ============================================================================
# NON-GOALS (explicit, from the task this script was written against):
#   - Does NOT touch or modify the existing nearest-centroid SSR
#     implementation or its output files. Called read-only, as a baseline.
#   - Does NOT touch the PCA-vs-fixed-band-index question.
#   - Does NOT decide whether to adopt adaptive-FCM SSR as a new default.
#     That is a decision for a future session, once real ABBY/SRER results
#     are in.
#
# ASSUMPTIONS -- CONFIRM THESE ON THE SERVER BEFORE TRUSTING A FULL RUN
# (also restated in the session's final report):
#   - fcm_fuzziness_w = 2.0 (Section 0): the standard/default FCM fuzziness
#     exponent in the literature. NOT re-derived or tuned for this data.
#   - c_min/c_max = 2/20 (Section 0): per the task's explicit instruction,
#     deliberately capped well under the old n_clusters=50 ceiling.
#   - fcm_n_subsample = 1000, fcm_n_reps = 5 (Section 0): SMALLER than the
#     existing k-means step's 2500 subsample / 20 reps. Fitting FCM for
#     19 candidate c values (2..20) per rep is far more expensive per rep
#     than one k-means fit, so both were lowered to keep overnight runtime
#     plausible. This is a real precision/cost tradeoff, not a neutral
#     default -- revisit once the real per-window timing (Section 7 output)
#     is in.
#   - hyperspec_dir / existing_ssr_csv / out_csv (Section 0): copied from
#     AnnualSpectralDiversity.R's own path conventions (absolute D: paths --
#     matches every other script in this repo's DataAnalysis pipeline, even
#     though it doesn't match CLAUDE.md's unfilled "relative paths only"
#     template, which no script in this repo currently follows). Confirm
#     these still point at the right places on the server.
#   - Python / fcmeans availability on the server is UNCONFIRMED -- see the
#     dependency note at the bottom of this header.
#
# ============================================================================
# DEPENDENCIES TO CONFIRM ON THE SERVER:
#   - R: reticulate (already a repo dependency via ComputeSpecBiodiv.R's
#     pyGNDiv bridge), terra, rhdf5, randomForest, cluster, dplyr, stringr,
#     purrr, tibble, readr -- all already used elsewhere in this repo.
#   - Python, reachable by reticulate: numpy, and `fuzzy-c-means` (PyPI
#     package name; imported in Python as `fcmeans`). Tested in the
#     sandbox this script was written in against fuzzy-c-means==2.0.2 --
#     NOT the same Python environment the server uses; reticulate must be
#     pointed (RETICULATE_PYTHON env var, or reticulate::use_python())
#     at a Python that has both installed, with a python shared library
#     reticulate can dynamically load (a statically-linked/pyenv-shim
#     python3 without libpythonX.Y.so will fail at py_config() -- ran into
#     exactly this in the sandbox; a system Python
#     (e.g. /usr/bin/python3-style install) with numpy + fuzzy-c-means
#     pip-installed worked).
#
# Runtime: dominated by the FCM step (Section 2/7) -- 19 candidate-c fits
# per rep x fcm_n_reps reps x every site-year, vs. one k-means fit x 20
# reps for baseline. Section 7 times both and prints a running multiplier
# so the ~10-20x per-window cost estimate (from replacing one k=50 fit with
# 19 fits at c<=20) can be checked against real timing rather than assumed.
# Results are written incrementally after every site-year (same crash-safe
# pattern as AnnualSpectralDiversity.R), so an interrupted overnight run
# still leaves a usable partial CSV.

library(terra)
library(rhdf5)
library(randomForest)
library(cluster)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readr)
library(reticulate)

# ============================================================================
# 0. Config -- CONFIRM before trusting a full run (see header)
# ============================================================================
hyperspec_dir    <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
existing_ssr_csv <- "D:/projects/moore/SpectralBiodiversity/Data/spectral_diversity_by_year.csv"
out_csv          <- "D:/projects/moore/SpectralBiodiversity/Data/ssr_comparison_kmeans_vs_adaptiveFCM.csv"

spectral_diversity_script <- "./Code/AnnualSpectralDiversity.R"
adaptive_fcm_py           <- "./Code/DataAnalysis/adaptive_fcm_ssr.py"
neonsites_path            <- "./Data/NEONsites.csv"

target_site_ids <- c("US-xAB", "US-xSR")   # ABBY, SRER -- the motivating comparison

c_min           <- 2
c_max           <- 20     # deliberately << the old n_clusters=50 -- see header
fcm_fuzziness_w <- 2.0    # ASSUMPTION -- see header
fcm_n_subsample <- 1000   # ASSUMPTION, < k-means' 2500 -- see header
fcm_n_reps      <- 5      # ASSUMPTION, < k-means' 20 reps -- see header
fcm_seed_base   <- 20260101   # arbitrary fixed base, only for run-to-run reproducibility

# Only the two inputs the SYNTHETIC gate (Section 4 below) actually needs are
# checked here. neonsites_path/hyperspec_dir are checked later (Section 4.5/5),
# AFTER the synthetic gate has run -- the gate must run standalone on
# synthetic data with no dependency on real Data/ being present (per the
# task this script was written against), so it cannot be blocked by a
# top-of-script check for real-data inputs.
required_inputs <- c(spectral_diversity_script, adaptive_fcm_py)
missing_inputs  <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Required input(s) not found:\n  ", paste(missing_inputs, collapse = "\n  "))
}

# ============================================================================
# 1. Import (not duplicate) the real extraction + baseline SSR pipeline
#    from AnnualSpectralDiversity.R -- see CORRECTION 1 in the header above
#    for why this is NOT FieldDiversity.R.
# ============================================================================
import_functions_from <- function(script_path, names, env = new.env(parent = parent.frame())) {
  exprs <- parse(script_path)
  for (nm in names) {
    found <- FALSE
    for (e in exprs) {
      if (is.call(e) && length(e) >= 3 && identical(e[[1]], as.name("<-")) &&
          is.name(e[[2]]) && identical(as.character(e[[2]]), nm)) {
        eval(e, envir = env)
        found <- TRUE
        break
      }
    }
    if (!found) {
      stop("import_functions_from(): could not find a top-level '", nm,
           " <- ...' binding in ", script_path,
           " -- has it been renamed/removed upstream? Not guessing.")
    }
  }
  env
}

# Called as spec_fns$name(...) at every call site below, NOT dumped into the
# global environment -- matches the established call style in ComputeLUE.R /
# ComputeLUE_Annual.R (the only two prior uses of import_functions_from() in
# this repo), which keeps every imported name's origin visible at its call
# site rather than silently shadowing/being shadowed by a same-named local.
spec_fns <- import_functions_from(
  spectral_diversity_script,
  c("bad_band_ranges", "read_neon_h5_tile", "mask_bad_bands", "nearest_band",
    "get_tower_reflectance", "compute_ndvi_raster", "extract_year_from_path",
    "discover_site_years", "assign_nearest_centroid",
    "compute_spectral_species_richness",
    # plain config constants, imported (not hand-copied with a "matches X"
    # comment, unlike buffer_m/ndvi_thresh below) specifically so the
    # baseline recompute (Section 6) can never silently drift from the real
    # pipeline's actual current calibration
    "n_clusters", "n_subsample_pixels", "n_pc_ssr", "n_reps_ssr")
)

# buffer_m/ndvi_thresh: hand-copied, not imported -- matches the established
# idiom in ComputeLUE.R/ComputeLUE_Annual.R (both hand-copy these same two
# constants from AnnualSpectralDiversity.R with a "matches AnnualSpectral-
# Diversity.R" comment rather than importing them). Revisit if
# AnnualSpectralDiversity.R's own values ever change.
buffer_m    <- 500    # matches AnnualSpectralDiversity.R's single site-year buffer
ndvi_thresh <- 0.4    # veg-mask threshold, matches AnnualSpectralDiversity.R

cat("Imported", length(ls(spec_fns)), "bindings from", spectral_diversity_script, "\n")
cat("  (buffer_m =", buffer_m, ", ndvi_thresh =", ndvi_thresh,
    ", n_clusters =", spec_fns$n_clusters, ", n_subsample_pixels =", spec_fns$n_subsample_pixels,
    ", n_pc_ssr =", spec_fns$n_pc_ssr, ", n_reps_ssr =", spec_fns$n_reps_ssr, ")\n")

# ============================================================================
# 2. Adaptive-FCM SSR: R-side wrapper around adaptive_fcm_ssr.py
# ============================================================================
source_python(adaptive_fcm_py)   # exposes adaptive_fcm_ssr(), fcm_validity_L() -- see CORRECTION 2

# get_window_pcs(): the SAME "mask -> values(na.rm=TRUE) -> prcomp(center=TRUE,
# scale.=FALSE) -> first n_pc columns" recipe compute_spectral_species_richness()
# uses internally (AnnualSpectralDiversity.R, Section 5), re-derived here
# rather than imported because that script exposes it only INLINE inside
# compute_spectral_species_richness(), not as a separate top-level function --
# there is nothing to import_functions_from(). Re-deriving these 4 lines (vs.
# the real H5/raster assembly this script does NOT re-derive) keeps both SSR
# methods clustering in the IDENTICAL PC-space, which matters here: the point
# of this comparison is to isolate "how many clusters / which clustering
# algorithm", not confound it with two different PCA bases. A cleaner
# long-term fix would be for AnnualSpectralDiversity.R to expose this as its
# own function; out of scope for this script (non-goal: don't touch that
# file).
get_window_pcs <- function(r, veg_mask, n_pc) {
  r_masked <- mask(r, veg_mask)
  vals <- values(r_masked, na.rm = TRUE)
  pca <- prcomp(vals, center = TRUE, scale. = FALSE)
  pca$x[, 1:n_pc]
}

# adaptive_fcm_ssr_window(): fcm_n_reps independent subsample-and-search
# reps over one window's PC-space pixel cloud, each rep calling into
# adaptive_fcm_ssr.py's adaptive_fcm_ssr() (candidate sweep c in
# [c_min, c_max], argmax_c L(c)). Mirrors
# compute_spectral_species_richness()'s "PCA once, subsample per rep"
# structure. Returns richness (mean winning c across reps, ROUNDED -- c is
# discrete, an unrounded mean is not a valid cluster count), the full
# per-rep winning-c vector (for cap-bunching diagnostics), and elapsed time.
adaptive_fcm_ssr_window <- function(pcs_all, c_min, c_max, fuzziness,
                                     n_subsample, n_reps, seed_base) {
  n_total <- nrow(pcs_all)
  if (n_total <= c_max) {
    return(list(richness = NA_real_, winning_c_per_rep = integer(0),
                hit_cap_frac = NA_real_, elapsed_sec = 0,
                status = paste0("skipped: only ", n_total,
                                 " pixels, needs > c_max (", c_max, ")")))
  }
  this_n_subsample <- min(n_subsample, n_total)

  winning_c_per_rep <- integer(n_reps)
  hit_cap_per_rep    <- logical(n_reps)
  t0 <- Sys.time()
  for (rep_i in seq_len(n_reps)) {
    sub_idx <- sample(seq_len(n_total), this_n_subsample)
    pcs_sub <- pcs_all[sub_idx, , drop = FALSE]
    res <- adaptive_fcm_ssr(pcs_sub, c_min = as.integer(c_min), c_max = as.integer(c_max),
                             fuzziness = fuzziness,
                             random_state = as.integer(seed_base + rep_i))
    winning_c_per_rep[rep_i] <- res$winning_c
    hit_cap_per_rep[rep_i]   <- res$hit_cap
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(richness = round(mean(winning_c_per_rep)),
       winning_c_per_rep = winning_c_per_rep,
       hit_cap_frac = mean(hit_cap_per_rep),
       elapsed_sec = elapsed,
       status = "ok")
}

# ============================================================================
# 3. SECTION 4 (task numbering) -- synthetic validation, HARD GATE
# ============================================================================
# Re-runs the same kind of synthetic scenarios used to diagnose the original
# k-means saturation bug (SESSION_LOG.md: "5 well-separated Gaussian blobs;
# a near-homogeneous single cluster", 10000 pixels, 4 PCs) through THIS
# script's actual adaptive_fcm_ssr_window() wrapper, at the SAME
# c_min/c_max/fuzziness/n_subsample/n_reps configured for the real run above
# (Section 0) -- not separately-tuned parameters, so this gate is a genuine
# test of the settings that will actually run on the server.
cat("\n==== SECTION 4: synthetic validation (hard gate, must pass before any real-data run) ====\n")

make_synthetic_blobs <- function(n_pixels = 10000, n_pc = 4, n_true_groups = 5,
                                  blob_sd = 2.0, seed = 1) {
  set.seed(seed)
  centers <- matrix(runif(n_true_groups * n_pc, -150, 150), nrow = n_true_groups)
  per_blob <- ceiling(n_pixels / n_true_groups)
  X <- do.call(rbind, lapply(seq_len(n_true_groups), function(i) {
    matrix(rnorm(per_blob * n_pc, sd = blob_sd), ncol = n_pc) +
      matrix(centers[i, ], nrow = per_blob, ncol = n_pc, byrow = TRUE)
  }))
  X[seq_len(n_pixels), , drop = FALSE]
}

make_synthetic_homogeneous <- function(n_pixels = 10000, n_pc = 4, sd = 2.0, seed = 2) {
  set.seed(seed)
  matrix(rnorm(n_pixels * n_pc, sd = sd), ncol = n_pc)
}

synthetic_n_pc <- spec_fns$n_pc_ssr   # reuse the SAME PC dimensionality as the real run

blob_pcs  <- make_synthetic_blobs(n_pixels = 10000, n_pc = synthetic_n_pc,
                                   n_true_groups = 5, seed = 1)
homog_pcs <- make_synthetic_homogeneous(n_pixels = 10000, n_pc = synthetic_n_pc, seed = 2)

cat("  Running adaptive-FCM SSR on 5-blob synthetic data (", nrow(blob_pcs), " px, ",
    synthetic_n_pc, " PCs)...\n", sep = "")
blob_result <- adaptive_fcm_ssr_window(blob_pcs, c_min, c_max, fcm_fuzziness_w,
                                        fcm_n_subsample, fcm_n_reps, fcm_seed_base)
cat("    winning c per rep:", paste(blob_result$winning_c_per_rep, collapse = ", "),
    "-> richness =", blob_result$richness,
    "(", round(blob_result$elapsed_sec, 1), "sec )\n")

cat("  Running adaptive-FCM SSR on homogeneous synthetic data (", nrow(homog_pcs), " px, ",
    synthetic_n_pc, " PCs)...\n", sep = "")
homog_result <- adaptive_fcm_ssr_window(homog_pcs, c_min, c_max, fcm_fuzziness_w,
                                         fcm_n_subsample, fcm_n_reps, fcm_seed_base + 1000L)
cat("    winning c per rep:", paste(homog_result$winning_c_per_rep, collapse = ", "),
    "-> richness =", homog_result$richness,
    "(", round(homog_result$elapsed_sec, 1), "sec )\n")

cat("\n  ---- Gate checks ----\n")
gate_separation_ok <- isTRUE(homog_result$richness < blob_result$richness)
gate_homog_low_ok  <- isTRUE(homog_result$richness <= 5)   # generous margin above the "ideally 1-3" bar
gate_blob_reasonable_ok <- isTRUE(blob_result$richness >= 3 && blob_result$richness <= 8)

cat("  1. homogeneous richness (", homog_result$richness,
    ") < blobs richness (", blob_result$richness, "): ",
    if (gate_separation_ok) "PASS" else "FAIL", "\n", sep = "")
cat("  2. homogeneous richness <= 5 (not near the c_max=", c_max, " cap): ",
    if (gate_homog_low_ok) "PASS" else "FAIL", "\n", sep = "")
cat("  3. (informational only, not gating) blobs richness in a plausible range of the true group",
    " count (5): ", blob_result$richness, " -- ",
    if (gate_blob_reasonable_ok) "within [3,8]" else "OUTSIDE [3,8], investigate before trusting real-data results",
    "\n", sep = "")

if (!gate_separation_ok || !gate_homog_low_ok) {
  stop("SYNTHETIC VALIDATION GATE FAILED -- the homogeneous case did not come in ",
       "meaningfully below the separated-blobs case (or landed near the c_max cap). ",
       "This means something is wrong with the validity-function implementation ",
       "(adaptive_fcm_ssr.py) or this wrapper, NOT the real data. STOPPING before ",
       "any real-data run, per the task's hard-gate instruction. Investigate ",
       "adaptive_fcm_ssr.py's fcm_validity_L() and the fcmeans fit before proceeding.")
}
cat("\n  SYNTHETIC VALIDATION GATE: PASSED. Proceeding to real-data comparison.\n\n")

# ============================================================================
# 4.5. Site crosswalk (reused pattern from every other script in this repo)
#      -- real-data inputs are checked HERE, not at the top of the script,
#      so the synthetic gate above can never be blocked by their absence.
# ============================================================================
if (!file.exists(neonsites_path)) {
  stop("Required input not found: ", neonsites_path,
       " -- this is expected in a sandbox with no server data access; the ",
       "synthetic gate above already ran and passed. Re-run this script on ",
       "the server, where neonsites_path points at the real file.")
}
towers_df <- read.csv(neonsites_path, fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site), Site.ID %in% target_site_ids)

if (nrow(towers_df) == 0) {
  stop("None of target_site_ids (", paste(target_site_ids, collapse = ", "),
       ") found in ", neonsites_path, " -- check Site.ID values match.")
}
cat("Target sites resolved:\n")
print(as.data.frame(towers_df[, c("Site.ID", "neon_site")]))

# ============================================================================
# 5. Discover site-year jobs for target sites only (mirrors
#    AnnualSpectralDiversity.R's own site_year_jobs construction, Section 7,
#    restricted here to target_site_ids)
# ============================================================================
if (!dir.exists(hyperspec_dir)) {
  stop("Hyperspectral tile directory not found: ", hyperspec_dir,
       " -- this is expected in a sandbox with no server data access; the ",
       "synthetic gate above already ran and passed. Re-run this script on ",
       "the server, where hyperspec_dir points at the real archive.")
}

site_year_jobs <- list()
for (tid in towers_df$Site.ID) {
  df <- spec_fns$discover_site_years(tid)
  if (nrow(df) == 0) next
  for (yr in sort(unique(df$year))) {
    site_year_jobs[[length(site_year_jobs) + 1]] <- list(
      tower_id = tid, year = yr, files = df$file[df$year == yr]
    )
  }
}
cat("Total site-year jobs for target sites:", length(site_year_jobs), "\n")
if (length(site_year_jobs) == 0) {
  stop("No site-years with hyperspectral data found under ", hyperspec_dir,
       " for target_site_ids -- nothing to process.")
}

# existing baseline SSR, if cached (see CORRECTION-adjacent header note:
# "pull existing values rather than recompute, if cached")
existing_ssr <- if (file.exists(existing_ssr_csv)) {
  read_csv(existing_ssr_csv, show_col_types = FALSE)
} else {
  cat("No cached SSR csv found at", existing_ssr_csv, "-- baseline will be recomputed for every site-year.\n")
  NULL
}

# ============================================================================
# 6. Main loop: one row per tower-YEAR (== one "window", see CORRECTION 3)
# ============================================================================
results <- tibble(
  tower_id = character(), neon_site = character(), year = character(),
  window_id = character(), n_veg_pixels = integer(),
  ssr_kmeans_fixed50 = double(), ssr_kmeans_baseline_source = character(),
  ssr_kmeans_runtime_sec = double(),
  ssr_adaptive_fcm = double(), fcm_winning_c_per_rep = character(),
  fcm_hit_cap_frac = double(), ssr_adaptive_fcm_runtime_sec = double(),
  status = character()
)

total_kmeans_runtime_sec <- 0
total_fcm_runtime_sec    <- 0

for (j in seq_along(site_year_jobs)) {
  job <- site_year_jobs[[j]]
  tower_id <- job$tower_id; yr <- job$year
  neon_site <- towers_df$neon_site[towers_df$Site.ID == tower_id][1]
  cat("\n==== [", j, "/", length(site_year_jobs), "] ", tower_id, " (", neon_site,
      ") - ", yr, " ====\n", sep = "")

  row <- tryCatch({
    data <- spec_fns$get_tower_reflectance(tower_id, job$files, buffer_m)
    if (is.null(data)) {
      tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
             window_id = "whole_buffer", n_veg_pixels = NA_integer_,
             ssr_kmeans_fixed50 = NA_real_, ssr_kmeans_baseline_source = NA_character_,
             ssr_kmeans_runtime_sec = NA_real_,
             ssr_adaptive_fcm = NA_real_, fcm_winning_c_per_rep = NA_character_,
             fcm_hit_cap_frac = NA_real_, ssr_adaptive_fcm_runtime_sec = NA_real_,
             status = "no reflectance data")
    } else {
      r <- data$raster; wl <- data$wavelengths
      ndvi <- spec_fns$compute_ndvi_raster(r, wl)
      veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)
      n_veg <- sum(!is.na(values(veg_mask)))

      # ---- baseline: pull cached value if present, else recompute (unchanged) ----
      cached_row <- NULL
      if (!is.null(existing_ssr)) {
        cached_row <- existing_ssr %>%
          filter(tower_id == !!tower_id, as.character(year) == !!yr)
      }
      if (!is.null(cached_row) && nrow(cached_row) >= 1 &&
          !is.na(cached_row$spectral_species_richness[1])) {
        ssr_kmeans <- cached_row$spectral_species_richness[1]
        kmeans_source <- "cached"
        kmeans_runtime <- NA_real_
        cat("  baseline SSR (fixed k=50): ", round(ssr_kmeans, 2), " [cached from ",
            existing_ssr_csv, "]\n", sep = "")
      } else {
        cat("  recomputing baseline SSR (fixed k=50, unchanged implementation)...\n")
        t0 <- Sys.time()
        km_result <- spec_fns$compute_spectral_species_richness(
          r, veg_mask, spec_fns$n_clusters, spec_fns$n_subsample_pixels,
          spec_fns$n_pc_ssr, spec_fns$n_reps_ssr)
        kmeans_runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        ssr_kmeans <- km_result$richness
        kmeans_source <- "recomputed"
        cat("    SSR:", round(ssr_kmeans, 2), "(", round(kmeans_runtime, 1), "sec )\n")
      }
      if (!is.na(kmeans_runtime)) total_kmeans_runtime_sec <<- total_kmeans_runtime_sec + kmeans_runtime

      # ---- adaptive-FCM SSR ----
      cat("  computing adaptive-FCM SSR (c in [", c_min, ",", c_max, "], w = ",
          fcm_fuzziness_w, ", ", fcm_n_reps, " reps)...\n", sep = "")
      pcs_all <- get_window_pcs(r, veg_mask, spec_fns$n_pc_ssr)
      fcm_result <- adaptive_fcm_ssr_window(pcs_all, c_min, c_max, fcm_fuzziness_w,
                                             fcm_n_subsample, fcm_n_reps, fcm_seed_base)
      total_fcm_runtime_sec <<- total_fcm_runtime_sec + fcm_result$elapsed_sec
      cat("    winning c per rep:", paste(fcm_result$winning_c_per_rep, collapse = ", "),
          "-> richness =", fcm_result$richness,
          " (hit_cap_frac =", round(fcm_result$hit_cap_frac, 2), ",",
          round(fcm_result$elapsed_sec, 1), "sec )\n")
      if (isTRUE(fcm_result$hit_cap_frac > 0)) {
        cat("    !! NOTE: winning c hit the c_max=", c_max, " cap in ",
            round(100 * fcm_result$hit_cap_frac), "% of reps for this window -- ",
            "possible sign the cap needs raising here.\n", sep = "")
      }

      tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
             window_id = "whole_buffer", n_veg_pixels = as.integer(n_veg),
             ssr_kmeans_fixed50 = ssr_kmeans, ssr_kmeans_baseline_source = kmeans_source,
             ssr_kmeans_runtime_sec = kmeans_runtime,
             ssr_adaptive_fcm = fcm_result$richness,
             fcm_winning_c_per_rep = paste(fcm_result$winning_c_per_rep, collapse = ";"),
             fcm_hit_cap_frac = fcm_result$hit_cap_frac,
             ssr_adaptive_fcm_runtime_sec = fcm_result$elapsed_sec,
             status = "success")
    }
  }, error = function(e) {
    cat("  unexpected error:", conditionMessage(e), "\n")
    tibble(tower_id = tower_id, neon_site = neon_site, year = yr,
           window_id = "whole_buffer", n_veg_pixels = NA_integer_,
           ssr_kmeans_fixed50 = NA_real_, ssr_kmeans_baseline_source = NA_character_,
           ssr_kmeans_runtime_sec = NA_real_,
           ssr_adaptive_fcm = NA_real_, fcm_winning_c_per_rep = NA_character_,
           fcm_hit_cap_frac = NA_real_, ssr_adaptive_fcm_runtime_sec = NA_real_,
           status = paste("error:", conditionMessage(e)))
  })

  results <- bind_rows(results, row)

  # checkpoint after EVERY site-year -- crash-safe, matches
  # AnnualSpectralDiversity.R's incremental-write convention
  write.csv(results, out_csv, row.names = FALSE)
}

cat("\nDone. Wrote", nrow(results), "rows to", out_csv, "\n")

# ============================================================================
# 7. Summary: ABBY vs. SRER separation, cap-bunching, runtime tradeoff
# ============================================================================
cat("\n==== SUMMARY ====\n")

ok_results <- results %>% filter(status == "success")

cat("\n---- ABBY vs. SRER separation (the original motivating comparison) ----\n")
if (nrow(ok_results) > 0) {
  site_summary <- ok_results %>%
    group_by(neon_site) %>%
    summarise(n_site_years = n(),
              mean_ssr_kmeans = mean(ssr_kmeans_fixed50, na.rm = TRUE),
              mean_ssr_adaptive_fcm = mean(ssr_adaptive_fcm, na.rm = TRUE),
              .groups = "drop")
  print(as.data.frame(site_summary))
  cat("\nFixed-k=50 SSR showed ABBY ~ SRER ~ 50 (no separation) in prior synthetic/",
      "real diagnostics -- compare mean_ssr_kmeans above to confirm that persists ",
      "on this real data. Adaptive-FCM SSR should show real separation between ",
      "ABBY (dense forest) and SRER (sparse desert) if it resolves the ",
      "saturation problem -- compare mean_ssr_adaptive_fcm above; a 2-6x gap ",
      "(matching CV/CHA/Rao's Q's prior separation) would be the strongest ",
      "positive signal.\n", sep = "")
} else {
  cat("No successful site-years to summarize.\n")
}

cat("\n---- c_max cap bunching ----\n")
cap_hits <- ok_results %>% filter(fcm_hit_cap_frac > 0)
if (nrow(cap_hits) > 0) {
  cat("!! ", nrow(cap_hits), " of ", nrow(ok_results),
      " site-years had at least one rep's winning c hit the cap (c_max = ", c_max,
      ") -- this is a sign c_max may need to be raised, NOT that 20 is confirmed ",
      "sufficient. See fcm_hit_cap_frac per row in ", out_csv, ".\n", sep = "")
  print(as.data.frame(cap_hits[, c("tower_id", "year", "fcm_hit_cap_frac")]))
} else {
  cat("No site-year hit the c_max =", c_max, "cap in any rep.\n")
}

cat("\n---- Runtime tradeoff ----\n")
cat("Total k-means baseline runtime (RECOMPUTED site-years only; cached ones",
    "contribute 0 here since they weren't timed this run):",
    round(total_kmeans_runtime_sec, 1), "sec\n")
cat("Total adaptive-FCM runtime (all site-years):", round(total_fcm_runtime_sec, 1), "sec\n")
if (total_kmeans_runtime_sec > 0) {
  cat("Adaptive-FCM / k-means runtime multiplier (recomputed site-years only):",
      round(total_fcm_runtime_sec / total_kmeans_runtime_sec, 1), "x\n")
} else {
  cat("(No site-years had their k-means baseline recomputed this run -- all were ",
      "cached -- so no direct same-run runtime multiplier is available. Compare ",
      "total_fcm_runtime_sec above against AnnualSpectralDiversity.R's own logged ",
      "per-site-year SSR timing instead.)\n", sep = "")
}
cat("\nSee", out_csv, "for the full per-site-year comparison, including per-window",
    "winning-c logs (fcm_winning_c_per_rep) for later inspection.\n")
