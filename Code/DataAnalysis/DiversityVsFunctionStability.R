#!/usr/bin/env Rscript
# Spectral diversity (Rao's Q) vs. temporal stability of annual ecosystem
# function (NEE, RECO, GPP, ET), one point per tower.
# Usage: Rscript Code/DataAnalysis/DiversityVsFunctionStability.R [spectral_csv] [lue_csv] [data_dir] [fig_dir]
#   spectral_csv default "./Data/spectral_diversity_by_year.csv"
#   lue_csv      default "./Data/lue_by_tower_year.csv" (produced by ComputeAnnualLUE.R
#                -- only its tower_id/year/nee_annual/reco_annual/gpp_annual/et_annual
#                columns are used; see Section 2 for why this source is still correct
#                even though the file name says "lue")
#   data_dir     default "./Data"    (diversity_stability_by_tower.csv,
#                                     diversity_stability_tests.csv,
#                                     diversity_stability_attrition.csv)
#   fig_dir      default "./Figures" (diversity_vs_function_stability.png)
# Also runnable line-by-line in RStudio with the project root as the working
# directory -- every default path below is relative to that root.
#
# ==============================================================================
# LUE REMOVED (this revision) -- a prior version plotted n ~ 8 towers in every
# panel (expected ~45). Root cause: that version filtered the ENTIRE
# lue_by_tower_year.csv to `status == "ok"` before splitting into the 5
# per-function stability calculations. But status in that file reflects
# ONLY whether LUE itself could be computed (it requires BOTH a flux match
# AND a QC-passed MODIS fPAR composite -- see
# Code/DataAnalysis/ComputeAnnualLUE.R's `status <- case_when(is.na(fpar_annual)
# ~ "missing fpar data", is.na(gpp_annual) | is.na(par_annual) ~ "missing flux
# data", apar_annual == 0 ~ "divide by zero", TRUE ~ "ok")`, ~line 362). A
# tower-year with complete, valid flux data but no usable MODIS fPAR composite
# got status "missing fpar data" and was dropped entirely by that filter --
# which makes sense for LUE (it genuinely needs fPAR) but is wrong for
# NEE/RECO/GPP/ET, which don't depend on fPAR/MODIS at all.
#
# Investigated whether nee_annual/reco_annual/gpp_annual/et_annual themselves
# are ever blanked or dropped for fPAR-missing tower-years: they are NOT.
# ComputeAnnualLUE.R builds `combined <- full_join(fpar_annual_tbl,
# flux_annual_tbl, by = c("tower_id", "year"))` (~line 353) -- a full_join, so
# every tower-year present in flux_annual_tbl (which comes from
# load_annual_flux_data(), reading the real annual AmeriFlux files directly --
# entirely independent of MODIS) is kept, with its flux columns carried
# through UNCHANGED regardless of whether a matching MODIS row exists. `status`
# and `lue_annual` are the only columns affected by fPAR availability (~line
# 359-369); gpp_annual/reco_annual/et_annual/nee_annual are never masked by it.
# So lue_by_tower_year.csv IS a correct source for these four flux variables,
# as long as nothing filters on `status` before using them -- confirmed by
# reading the script, not run against real data (this sandbox still has no
# Data/ directory; see SESSION_LOG.md). If the real file's flux columns turn
# out to be sparser than this reasoning predicts, the attrition table below
# (Section 2) will show it directly, and the fix at that point would be to
# read the annual AmeriFlux files directly instead of trusting this file --
# not implemented here since the join-level reasoning above indicates it
# shouldn't be necessary.
#
# LUE, and the MODIS/fPAR/PAR dependence that came with it, is removed
# entirely from this script: no lue_annual CV, no implausibility check, no
# fifth plot column. The figure is now 2 rows (Rao's Q NIRv, Rao's Q
# all-bands) x 4 columns (NEE, RECO, GPP, ET) = 8 panels, 8 tests.
# ==============================================================================

required_pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "stringr", "purrr")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(stringr); library(purrr)
})

args         <- commandArgs(trailingOnly = TRUE)
spectral_csv <- if (length(args) >= 1) args[1] else "./Data/spectral_diversity_by_year.csv"
lue_csv      <- if (length(args) >= 2) args[2] else "./Data/lue_by_tower_year.csv"
data_dir     <- if (length(args) >= 3) args[3] else "./Data"
fig_dir      <- if (length(args) >= 4) args[4] else "./Figures"

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

cv_mean_eps  <- 1e-8  # |mean| below this is treated as "zero" for CV% (guards div-by-~0)
CV_MIN_YEARS <- 3     # minimum non-NA flux years required to compute a tower's CV for a function

diversity_map <- tibble::tibble(
  col   = c("raoq_nirv", "raoq_allbands"),
  label = c("Rao's Q (NIRv)", "Rao's Q (all bands)")
)
function_map <- tibble::tibble(
  col     = c("nee", "reco", "gpp", "et"),
  label   = c("NEE", "RECO", "GPP", "ET"),
  src_col = c("nee_annual", "reco_annual", "gpp_annual", "et_annual")
)

attrition <- tibble::tibble(stage = character(), description = character(), n_towers = integer())
add_attrition <- function(stage, description, n_towers) {
  attrition <<- bind_rows(attrition, tibble::tibble(stage = stage, description = description, n_towers = n_towers))
}

# ============================================================================
# 1. Spectral diversity input -- columns confirmed against
#    Code/AnnualSpectralDiversity.R's `results` tibble (~line 966) and its
#    status assignment (~line 1038: `status <- if (tower_id %in%
#    single_year_tower_ids) "success (single-year site)" else "success"`) --
#    so BOTH values count as a successful row, not just "success".
# ============================================================================
if (!file.exists(spectral_csv)) stop("Input csv not found: ", spectral_csv)

spectral_required_cols <- c("tower_id", "neon_site", "year", "raoq_nirv", "raoq_allbands", "status")
spectral_raw <- read_csv(spectral_csv, show_col_types = FALSE)
missing_cols <- setdiff(spectral_required_cols, names(spectral_raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in ", spectral_csv, ": ", paste(missing_cols, collapse = ", "),
       "\nColumns present in file: ", paste(names(spectral_raw), collapse = ", "))
}
spectral_raw <- spectral_raw %>% mutate(year = as.integer(year))

cat("\n=== spectral_diversity_by_year.csv status values ===\n")
print(table(spectral_raw$status, useNA = "ifany"))

spectral_success_vals <- c("success", "success (single-year site)")
spectral <- spectral_raw %>% filter(status %in% spectral_success_vals)
cat("Kept ", nrow(spectral), " successful spectral row(s) (status %in% c(",
    paste(sprintf('"%s"', spectral_success_vals), collapse = ", "), ")).\n", sep = "")

dup_spec <- spectral %>% count(tower_id, year) %>% filter(n > 1)
if (nrow(dup_spec) > 0) {
  offending <- spectral %>% semi_join(dup_spec, by = c("tower_id", "year")) %>% arrange(tower_id, year)
  print(offending)
  stop("Duplicate tower_id x year rows found in ", spectral_csv,
       " (per-tower diversity mean assumes one value per tower-year). See offending rows above.")
}

spectral_valid <- spectral %>% filter(!is.na(raoq_nirv))
n_spectral_towers <- n_distinct(spectral_valid$tower_id)
add_attrition("spectral", "Towers in spectral_diversity_by_year.csv with status success/single-year AND raoq_nirv non-NA", n_spectral_towers)

# ============================================================================
# 2. Flux input -- columns confirmed by reading
#    Code/DataAnalysis/ComputeAnnualLUE.R's own output `transmute()` call
#    (~line 372): tower_id, year, fpar_annual, n_fpar_obs_used,
#    n_fpar_obs_total, gpp_annual, reco_annual, et_annual, nee_annual,
#    par_annual, lue_annual, status. Only tower_id/year/nee_annual/
#    reco_annual/gpp_annual/et_annual are used below -- status is NOT
#    filtered on (see the file header for why) and fpar/par/lue columns are
#    not required. NEE_VUT_REF (AmeriFlux) is negative for net uptake -- that
#    sign convention is NOT flipped anywhere in ComputeAnnualLUE.R or here.
# ============================================================================
if (!file.exists(lue_csv)) {
  stop(lue_csv, " not found. Run Code/DataAnalysis/ComputeAnnualLUE.R first to produce it.")
}

flux_required_cols <- c("tower_id", "year", "nee_annual", "reco_annual", "gpp_annual", "et_annual")
flux_raw <- read_csv(lue_csv, show_col_types = FALSE)
missing_cols <- setdiff(flux_required_cols, names(flux_raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in ", lue_csv, ": ", paste(missing_cols, collapse = ", "),
       "\nColumns present in file: ", paste(names(flux_raw), collapse = ", "))
}
flux_raw <- flux_raw %>% mutate(year = as.integer(year))

if ("status" %in% names(flux_raw)) {
  cat("\n=== lue_by_tower_year.csv status values (NOT used to filter nee/reco/gpp/et -- see file header) ===\n")
  print(table(flux_raw$status, useNA = "ifany"))
} else {
  cat("\nNo status column in ", lue_csv, " -- not required, not used.\n", sep = "")
}
cat("NEE sign convention: negative = net uptake (AmeriFlux convention), not flipped.\n")

dup_flux <- flux_raw %>% count(tower_id, year) %>% filter(n > 1)
if (nrow(dup_flux) > 0) {
  offending <- flux_raw %>% semi_join(dup_flux, by = c("tower_id", "year")) %>% arrange(tower_id, year)
  print(offending)
  stop("Duplicate tower_id x year rows found in ", lue_csv,
       " (per-tower stability CV assumes one value per tower-year). See offending rows above.")
}

n_flux_towers <- n_distinct(flux_raw$tower_id)
add_attrition("flux", "Towers in lue_by_tower_year.csv (any status -- status not used to filter)", n_flux_towers)

# ---- intersection + ID-mismatch check ----
spectral_towers_set <- unique(spectral_valid$tower_id)
flux_towers_set     <- unique(flux_raw$tower_id)
n_intersection <- length(intersect(spectral_towers_set, flux_towers_set))
add_attrition("intersection", "Towers present in both (exact tower_id match)", n_intersection)

only_in_spectral <- setdiff(spectral_towers_set, flux_towers_set)
only_in_flux     <- setdiff(flux_towers_set, spectral_towers_set)
cat("\n=== Towers in only one file ===\n")
cat("Only in spectral (", length(only_in_spectral), "): ", paste(sort(only_in_spectral), collapse = ", "), "\n", sep = "")
cat("Only in flux (", length(only_in_flux), "): ", paste(sort(only_in_flux), collapse = ", "), "\n", sep = "")

normalize_id <- function(x) sub("^US-", "", trimws(toupper(x)))
mismatch_pairs <- tibble::tibble(spectral_id = character(), flux_id = character())
for (a in only_in_spectral) {
  hit <- only_in_flux[normalize_id(only_in_flux) == normalize_id(a)]
  if (length(hit) > 0) {
    mismatch_pairs <- bind_rows(mismatch_pairs, tibble::tibble(spectral_id = a, flux_id = hit))
  }
}
if (nrow(mismatch_pairs) > 0) {
  cat("\n!! Likely ID MISMATCHES (same tower, different case/whitespace/\"US-\" prefix spelling):\n")
  print(as.data.frame(mismatch_pairs))
} else if (length(only_in_spectral) > 0 || length(only_in_flux) > 0) {
  cat("\nNo case/whitespace/\"US-\"-prefix mismatches detected among the IDs found in only one file ",
      "(they appear to be genuinely different towers, not spelling variants).\n", sep = "")
}

# ============================================================================
# 3. Ecosystem mapping -- COPIED VERBATIM (bin map, levels, palette, file
#    fallback, join key, unmapped-code stop()) from
#    Code/DataAnalysis/RaoQYearVariability.R Section 2, per instruction to
#    reuse that exact logic so the two scripts' ecosystem figures stay
#    consistent. If that script's version ever diverges from this block,
#    RaoQYearVariability.R is the source of truth -- not this copy.
# ============================================================================
bin_map    <- c(DBF = "Forest", EBF = "Forest", ENF = "Forest", MF = "Forest",
                 GRA = "Grassland",
                 OSH = "Shrub/Savanna", SAV = "Shrub/Savanna",
                 CRO = "Cropland/Mosaic", CVM = "Cropland/Mosaic",
                 WET = "Wetland")
bin_levels <- c("Forest", "Grassland", "Shrub/Savanna", "Cropland/Mosaic", "Wetland")
bin_colors <- c(Forest = "#2a9d8f", Grassland = "#e9c46a",
                 `Shrub/Savanna` = "#f4a261", `Cropland/Mosaic` = "#8ab17d",
                 Wetland = "#277da1")

sites_candidates <- c("./Data/NEONsites.csv", "./NEONsites_Footprints.csv")
sites_path <- NA_character_
veg_col    <- NA_character_
for (p in sites_candidates) {
  if (!file.exists(p)) next
  hdr <- names(read.csv(p, fileEncoding = "UTF-8-BOM", nrows = 0, check.names = TRUE))
  vcol <- hdr[grepl("veg.?type", hdr, ignore.case = TRUE)]
  if (length(vcol) == 1) {
    sites_path <- p
    veg_col <- vcol
    break
  }
}
if (is.na(sites_path)) {
  stop("Could not find a Veg Type column in either candidate ecosystem file: ",
       paste(sites_candidates, collapse = ", "),
       ". Checked each for existence and for exactly one column name matching ",
       "/veg.?type/i. Not guessing which file or column to use.")
}
cat("\nUsing ecosystem/Veg Type file:", sites_path, "(column '", veg_col, "')\n", sep = "")

sites_raw <- read.csv(sites_path, fileEncoding = "UTF-8-BOM")

observed_codes <- unique(sites_raw[[veg_col]])
unmapped_codes <- setdiff(observed_codes, names(bin_map))
if (length(unmapped_codes) > 0) {
  stop("'", veg_col, "' in ", sites_path, " contains IGBP code(s) not covered by the ",
       "Forest/Grassland/Shrub-Savanna/Cropland-Mosaic/Wetland bin mapping: ",
       paste(unmapped_codes, collapse = ", "),
       ". Extend the mapping deliberately rather than dropping these sites silently.")
}

sites_lookup <- sites_raw %>%
  transmute(tower_id = Site.ID,
            ecosystem_bin = bin_map[.data[[veg_col]]]) %>%
  distinct(tower_id, .keep_all = TRUE)

# ============================================================================
# 4. Per-tower diversity (mean across all available spectral years) and
#    per-tower stability (temporal CV% across all available flux years, no
#    status filter -- just per-function non-NA). These two periods are NOT
#    aligned -- a tower's diversity mean can come from entirely different
#    years than its stability CV; see the figure caption and console below.
# ============================================================================
diversity_tbl <- spectral_valid %>%
  group_by(tower_id) %>%
  summarise(n_spectral_years = n(),
            neon_site = dplyr::first(neon_site),
            mean_raoq_nirv = mean(raoq_nirv, na.rm = TRUE),
            mean_raoq_allbands = mean(raoq_allbands, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(mean_raoq_nirv, mean_raoq_allbands), ~ ifelse(is.nan(.x), NA_real_, .x)))

stability_long <- flux_raw %>%
  select(tower_id, year, nee = nee_annual, reco = reco_annual, gpp = gpp_annual, et = et_annual) %>%
  pivot_longer(cols = c(nee, reco, gpp, et), names_to = "function_name", values_to = "value")

stability_tbl <- stability_long %>%
  filter(!is.na(value)) %>%
  group_by(tower_id, function_name) %>%
  filter(n() >= CV_MIN_YEARS) %>%
  summarise(n_years_used = n(), mean = mean(value), sd = sd(value),
            any_pos = any(value > 0), any_neg = any(value < 0), .groups = "drop") %>%
  mutate(
    sign_change = any_pos & any_neg,
    cv_flagged = (abs(mean) < sd) | sign_change,
    cv_raw = 100 * sd / abs(mean),
    cv_guarded_zero = !is.finite(cv_raw) | abs(mean) < cv_mean_eps,
    cv_pct = ifelse(cv_guarded_zero, NA_real_, cv_raw)
  )

cat("\n=== Per-function stability summary (>= ", CV_MIN_YEARS, " non-NA years, no status filter) ===\n", sep = "")
for (fn in function_map$col) {
  sub <- stability_tbl %>% filter(function_name == fn)
  n_guarded <- sum(sub$cv_guarded_zero & !sub$sign_change)
  flagged <- sub %>% filter(cv_flagged)
  cat(sprintf("%s: %d tower(s) with >= %d years; %d mean~0-guarded; %d flagged unstable",
              toupper(fn), nrow(sub), CV_MIN_YEARS, n_guarded, nrow(flagged)))
  if (nrow(flagged) > 0) {
    cat(" (", paste(flagged$tower_id, collapse = ", "), ")", sep = "")
  }
  cat("\n")
  add_attrition(paste0("flux_years_", fn),
                paste0(toupper(fn), ": towers with >= ", CV_MIN_YEARS, " non-NA years (no status filter)"),
                nrow(sub))
}

cat("\n=== Attrition table ===\n")
print(attrition)
write_csv(attrition, file.path(data_dir, "diversity_stability_attrition.csv"))

stability_wide <- stability_tbl %>%
  select(tower_id, function_name, n_years_used, mean, sd, cv_pct, cv_flagged) %>%
  pivot_wider(id_cols = tower_id, names_from = function_name,
              values_from = c(n_years_used, mean, sd, cv_pct, cv_flagged),
              names_glue = "{function_name}_{.value}")

all_towers <- union(diversity_tbl$tower_id, flux_towers_set)
master <- tibble(tower_id = all_towers) %>%
  left_join(diversity_tbl, by = "tower_id") %>%
  left_join(stability_wide, by = "tower_id") %>%
  left_join(sites_lookup, by = "tower_id")

unmatched_towers <- master %>% filter(is.na(ecosystem_bin)) %>% pull(tower_id)
if (length(unmatched_towers) > 0) {
  cat("\n!! ", length(unmatched_towers), " tower(s) with no ecosystem match in ", sites_path,
      " (kept with ecosystem_bin = NA, excluded from the figure):\n  ",
      paste(unmatched_towers, collapse = ", "), "\n", sep = "")
} else {
  cat("\nAll towers matched an ecosystem bin.\n")
}

out_cols <- c("tower_id", "neon_site", "ecosystem_bin", "n_spectral_years",
              "mean_raoq_nirv", "mean_raoq_allbands")
for (fn in function_map$col) {
  out_cols <- c(out_cols, paste0(fn, c("_n_years_used", "_mean", "_sd", "_cv_pct", "_cv_flagged")))
}
write_csv(master %>% select(all_of(out_cols)), file.path(data_dir, "diversity_stability_by_tower.csv"))

# ============================================================================
# 5. Rank-order tests (Spearman rho, Kendall tau) for each of the 8
#    diversity-metric x function panels, BH FDR within each diversity metric
#    across its 4 functions (figure shows raw rho/p only -- FDR columns are
#    CSV-only, per instruction).
# ============================================================================
compute_rank_tests <- function(x, y) {
  ok <- complete.cases(x, y)
  n_ok <- sum(ok)
  if (n_ok < 5) {
    return(tibble(n = n_ok, rho = NA_real_, rho_p = NA_real_, tau = NA_real_, tau_p = NA_real_))
  }
  sp <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  kt <- cor.test(x[ok], y[ok], method = "kendall", exact = FALSE)
  tibble(n = n_ok, rho = unname(sp$estimate), rho_p = sp$p.value,
         tau = unname(kt$estimate), tau_p = kt$p.value)
}

panel_grid <- tibble()
for (dm in diversity_map$col) for (fn in function_map$col) {
  panel_grid <- bind_rows(panel_grid, tibble(diversity_metric = dm, `function` = fn))
}

test_results <- panel_grid %>%
  rowwise() %>%
  mutate(compute_rank_tests(master[[paste0("mean_", diversity_metric)]], master[[paste0(`function`, "_cv_pct")]])) %>%
  ungroup() %>%
  group_by(diversity_metric) %>%
  mutate(rho_p_fdr = p.adjust(rho_p, method = "BH"), tau_p_fdr = p.adjust(tau_p, method = "BH")) %>%
  ungroup()

cat("\n=== Rank-order tests (raw p; BH FDR within each diversity metric; 8 tests total) ===\n")
print(test_results)

write_csv(
  test_results %>% select(diversity_metric, `function`, n, rho, rho_p, tau, tau_p, rho_p_fdr, tau_p_fdr),
  file.path(data_dir, "diversity_stability_tests.csv")
)

# ============================================================================
# 6. Figure -- 2 (diversity metrics) x 4 (functions) = 8 panels, built as
#    individual ggplot objects and assembled with patchwork (not
#    facet_grid/facet_wrap): each panel needs its own log10-vs-linear y-axis
#    decision (a per-function property shared down both rows of that
#    column), while x is shared within a row (same diversity metric) -- a
#    combination facet_grid's row/column-scale semantics can't express
#    directly, since its free_x/free_y are tied to columns/rows respectively,
#    the opposite pairing needed here.
# ============================================================================
shape_flagged <- c(`FALSE` = 16, `TRUE` = 1)  # filled vs. hollow

# y-axis scale/log decision is a per-FUNCTION (column) property: computed
# from ALL towers with a valid cv_pct for that function, independent of which
# diversity metric it's paired with in a given panel, so both rows of a
# column share the same y treatment.
function_y_info <- function_map %>%
  rowwise() %>%
  mutate({
    cv_col <- paste0(`col`, "_cv_pct")
    cv_pos <- master[[cv_col]][!is.na(master[[cv_col]]) & master[[cv_col]] > 0]
    use_log <- length(cv_pos) > 0 && (max(cv_pos) / min(cv_pos) > 50)
    tibble(use_log = use_log)
  }) %>%
  ungroup()

cat("\n=== Figure y-axis scale per function ===\n")
for (i in seq_len(nrow(function_y_info))) {
  cat(function_y_info$label[i], ": ", if (function_y_info$use_log[i]) "log10 (CV spans > 50x)" else "linear", "\n", sep = "")
}

# x range is a per-DIVERSITY-METRIC (row) property, shared across that row's
# 4 panels regardless of each panel's own NA pattern.
diversity_x_range <- diversity_map %>%
  rowwise() %>%
  mutate({
    v <- master[[paste0("mean_", col)]]
    rng <- range(v, na.rm = TRUE)
    pad <- 0.05 * diff(rng)
    tibble(xmin = rng[1] - pad, xmax = rng[2] + pad)
  }) %>%
  ungroup()

fmt_p_raw <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))

make_panel <- function(diversity_col, function_col) {
  d_label <- diversity_map$label[diversity_map$col == diversity_col]
  f_label <- function_map$label[function_map$col == function_col]
  x_rng   <- diversity_x_range %>% filter(col == diversity_col)
  use_log <- function_y_info$use_log[function_y_info$col == function_col]

  d <- master %>%
    transmute(x = .data[[paste0("mean_", diversity_col)]],
              y = .data[[paste0(function_col, "_cv_pct")]],
              flagged = factor(.data[[paste0(function_col, "_cv_flagged")]], levels = c(FALSE, TRUE)),
              ecosystem_bin = factor(ecosystem_bin, levels = bin_levels)) %>%
    filter(!is.na(x), !is.na(y))

  r <- test_results %>% filter(diversity_metric == diversity_col, `function` == function_col)
  subtitle_txt <- if (is.na(r$n) || r$n < 5) {
    "insufficient data"
  } else {
    sprintf("Spearman rho = %.2f (%s, n = %d)", r$rho, fmt_p_raw(r$rho_p), r$n)
  }

  p <- ggplot(d, aes(x = x, y = y)) +
    labs(title = f_label, subtitle = subtitle_txt,
         x = d_label, y = if (use_log) paste0(f_label, " CV (%) [log10]") else paste0(f_label, " CV (%)")) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7.5, colour = "grey30")) +
    coord_cartesian(xlim = c(x_rng$xmin, x_rng$xmax))

  if (nrow(d) > 0) {
    p <- p +
      geom_smooth(method = "lm", se = TRUE, formula = y ~ x, colour = "#c0392b", linewidth = 0.6, na.rm = TRUE) +
      geom_point(aes(colour = ecosystem_bin, shape = flagged), size = 2.4, alpha = 0.8, na.rm = TRUE) +
      scale_colour_manual(values = bin_colors, drop = FALSE) +
      scale_shape_manual(values = shape_flagged, drop = FALSE)
  }
  if (use_log) p <- p + scale_y_log10()
  p
}

panels <- list()
for (dm in diversity_map$col) for (fn in function_map$col) {
  panels[[length(panels) + 1]] <- make_panel(dm, fn) +
    guides(colour = guide_legend(title = "Ecosystem", override.aes = list(size = 3)),
           shape = guide_legend(title = "Flagged\nunstable CV", override.aes = list(size = 3)))
}

# plot_layout(guides = "collect") merges the 8 panels' identical
# colour/shape scales into one shared legend (patchwork's standard
# de-duplication) instead of rendering a legend inside every panel.
fig <- wrap_plots(panels, nrow = 2, ncol = 4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Spectral diversity vs. temporal stability of annual ecosystem function",
    caption = str_wrap(paste(
      "Stability = temporal CV (%) across years (lower = more stable); diversity = per-tower mean",
      "Rao's Q over its available spectral years. Years are NOT aligned between the two datasets --",
      "diversity and stability can reflect different year ranges for the same tower. Hollow points =",
      "CV unreliable (|mean| < sd, i.e. CV > 100%, or the annual values changed sign across years).",
      "NEE is negative for net uptake (AmeriFlux sign convention), not flipped. Raw p-values shown;",
      "8 tests, no multiplicity correction in this figure (see diversity_stability_tests.csv for BH FDR)."
    ), width = 150),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0))
  ) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "diversity_vs_function_stability.png"), fig,
       width = 15.5, height = 6.5, dpi = 200, bg = "white", limitsize = FALSE)

cat("\nSaved CSVs to", data_dir, "and figure to", fig_dir, "\n")
