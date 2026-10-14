#!/usr/bin/env Rscript
# Spectral diversity (Rao's Q) vs. temporal stability of annual ecosystem
# function (NEE, RECO, GPP, ET, LUE), one point per tower.
# Usage: Rscript Code/DataAnalysis/DiversityVsFunctionStability.R [spectral_csv] [lue_csv] [data_dir] [fig_dir]
#   spectral_csv default "./Data/spectral_diversity_by_year.csv"
#   lue_csv      default "./Data/lue_by_tower_year.csv" (produced by ComputeAnnualLUE.R)
#   data_dir     default "./Data"    (diversity_stability_by_tower.csv, diversity_stability_tests.csv)
#   fig_dir      default "./Figures" (diversity_vs_function_stability.png)
# Also runnable line-by-line in RStudio with the project root as the working
# directory -- every default path below is relative to that root.

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

cv_mean_eps <- 1e-8  # |mean| below this is treated as "zero" for CV% (guards div-by-~0)

diversity_map <- tibble::tibble(
  col   = c("raoq_nirv", "raoq_allbands"),
  label = c("Rao's Q (NIRv)", "Rao's Q (all bands)")
)
function_map <- tibble::tibble(
  col     = c("nee", "reco", "gpp", "et", "lue"),
  label   = c("NEE", "RECO", "GPP", "ET", "LUE"),
  src_col = c("nee_annual", "reco_annual", "gpp_annual", "et_annual", "lue_annual")
)

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

# ============================================================================
# 2. LUE/flux input -- columns confirmed by reading
#    Code/DataAnalysis/ComputeAnnualLUE.R's own output `transmute()` call
#    (~line 372): tower_id, year, fpar_annual, n_fpar_obs_used,
#    n_fpar_obs_total, gpp_annual, reco_annual, et_annual, nee_annual,
#    par_annual, lue_annual, status -- status in {"ok", "missing fpar data",
#    "missing flux data", "divide by zero"} (same script, ~line 362-367).
#    NEE_VUT_REF (AmeriFlux) is negative for net uptake -- that sign
#    convention is NOT flipped anywhere in ComputeAnnualLUE.R or here.
# ============================================================================
if (!file.exists(lue_csv)) {
  stop(lue_csv, " not found. Run Code/DataAnalysis/ComputeAnnualLUE.R first to produce it.")
}

lue_required_cols <- c("tower_id", "year", "fpar_annual", "n_fpar_obs_used", "n_fpar_obs_total",
                        "gpp_annual", "reco_annual", "et_annual", "nee_annual", "par_annual",
                        "lue_annual", "status")
lue_raw <- read_csv(lue_csv, show_col_types = FALSE)
missing_cols <- setdiff(lue_required_cols, names(lue_raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in ", lue_csv, ": ", paste(missing_cols, collapse = ", "),
       "\nColumns present in file: ", paste(names(lue_raw), collapse = ", "))
}
lue_raw <- lue_raw %>% mutate(year = as.integer(year))

cat("\n=== lue_by_tower_year.csv status values ===\n")
print(table(lue_raw$status, useNA = "ifany"))

lue <- lue_raw %>% filter(status == "ok")
cat("Kept ", nrow(lue), " tower-year(s) with status == \"ok\" (NEE sign convention: negative = net ",
    "uptake, AmeriFlux convention, not flipped).\n", sep = "")

dup_lue <- lue %>% count(tower_id, year) %>% filter(n > 1)
if (nrow(dup_lue) > 0) {
  offending <- lue %>% semi_join(dup_lue, by = c("tower_id", "year")) %>% arrange(tower_id, year)
  print(offending)
  stop("Duplicate tower_id x year rows found in ", lue_csv,
       " (per-tower stability CV assumes one value per tower-year). See offending rows above.")
}

lue_med <- median(lue$lue_annual, na.rm = TRUE)
cat("\nlue_annual (status == \"ok\"): min = ", round(min(lue$lue_annual, na.rm = TRUE), 4),
    ", median = ", round(lue_med, 4), ", max = ", round(max(lue$lue_annual, na.rm = TRUE), 4), "\n", sep = "")
implausible_lue <- lue %>% filter(!is.na(lue_annual) & (lue_annual < 0 | lue_annual > 10 * lue_med))
if (nrow(implausible_lue) > 0) {
  cat("!! ", nrow(implausible_lue), " physically implausible lue_annual value(s) (negative or > 10x ",
      "the median) -- kept in the analysis, not dropped:\n", sep = "")
  print(as.data.frame(implausible_lue %>% select(tower_id, year, lue_annual)))
} else {
  cat("No implausible lue_annual values (negative or > 10x median) found.\n")
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
#    per-tower stability (temporal CV% across all available flux years).
#    These two periods are NOT aligned -- a tower's diversity mean can come
#    from entirely different years than its stability CV; see the figure
#    caption and the console summary below.
# ============================================================================
diversity_tbl <- spectral %>%
  group_by(tower_id) %>%
  summarise(n_spectral_years = n(),
            neon_site = dplyr::first(neon_site),
            mean_raoq_nirv = mean(raoq_nirv, na.rm = TRUE),
            mean_raoq_allbands = mean(raoq_allbands, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(mean_raoq_nirv, mean_raoq_allbands), ~ ifelse(is.nan(.x), NA_real_, .x)))

n_spectral_towers <- n_distinct(spectral$tower_id)
n_flux_towers     <- n_distinct(lue$tower_id)
n_intersection    <- length(intersect(spectral$tower_id, lue$tower_id))
cat("\n=== Tower counts ===\n")
cat("n spectral towers:", n_spectral_towers, "\n")
cat("n flux (LUE, status == \"ok\") towers:", n_flux_towers, "\n")
cat("n intersection (both diversity and flux data available):", n_intersection, "\n")

stability_long <- lue %>%
  select(tower_id, year, nee = nee_annual, reco = reco_annual, gpp = gpp_annual,
         et = et_annual, lue = lue_annual) %>%
  pivot_longer(cols = c(nee, reco, gpp, et, lue), names_to = "function_name", values_to = "value")

stability_tbl <- stability_long %>%
  filter(!is.na(value)) %>%
  group_by(tower_id, function_name) %>%
  filter(n() >= 3) %>%
  summarise(n_years_used = n(), mean = mean(value), sd = sd(value),
            any_pos = any(value > 0), any_neg = any(value < 0), .groups = "drop") %>%
  mutate(
    sign_change = any_pos & any_neg,
    cv_flagged = (abs(mean) < sd) | sign_change,
    cv_raw = 100 * sd / abs(mean),
    cv_guarded_zero = !is.finite(cv_raw) | abs(mean) < cv_mean_eps,
    cv_pct = ifelse(cv_guarded_zero, NA_real_, cv_raw)
  )

cat("\n=== Per-function stability summary ===\n")
for (fn in function_map$col) {
  sub <- stability_tbl %>% filter(function_name == fn)
  n_guarded <- sum(sub$cv_guarded_zero & !sub$sign_change)
  flagged <- sub %>% filter(cv_flagged)
  cat(sprintf("%s: %d tower(s) with >= 3 years; %d mean~0-guarded; %d flagged unstable",
              toupper(fn), nrow(sub), n_guarded, nrow(flagged)))
  if (nrow(flagged) > 0) {
    cat(" (", paste(flagged$tower_id, collapse = ", "), ")", sep = "")
  }
  cat("\n")
}

stability_wide <- stability_tbl %>%
  select(tower_id, function_name, n_years_used, mean, sd, cv_pct, cv_flagged) %>%
  pivot_wider(id_cols = tower_id, names_from = function_name,
              values_from = c(n_years_used, mean, sd, cv_pct, cv_flagged),
              names_glue = "{function_name}_{.value}")

all_towers <- union(diversity_tbl$tower_id, unique(lue$tower_id))
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
# 5. Rank-order tests (Spearman rho, Kendall tau) for each of the 10
#    diversity-metric x function panels, BH FDR within each diversity metric
#    across its 5 functions (figure shows raw rho/p only -- FDR columns are
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

cat("\n=== Rank-order tests (raw p; BH FDR within each diversity metric) ===\n")
print(test_results)

write_csv(
  test_results %>% select(diversity_metric, `function`, n, rho, rho_p, tau, tau_p, rho_p_fdr, tau_p_fdr),
  file.path(data_dir, "diversity_stability_tests.csv")
)

# ============================================================================
# 6. Figure -- 2 (diversity metrics) x 5 (functions) = 10 panels, built as
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
# 5 panels regardless of each panel's own NA pattern.
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

# plot_layout(guides = "collect") merges the 10 panels' identical
# colour/shape scales into one shared legend (patchwork's standard
# de-duplication) instead of rendering a legend inside every panel.
fig <- wrap_plots(panels, nrow = 2, ncol = 5) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Spectral diversity vs. temporal stability of annual ecosystem function",
    caption = str_wrap(paste(
      "Stability = temporal CV (%) across years (lower = more stable); diversity = per-tower mean",
      "Rao's Q over its available spectral years. Years are NOT aligned between the two datasets --",
      "diversity and stability can reflect different year ranges for the same tower. Hollow points =",
      "CV unreliable (|mean| < sd, i.e. CV > 100%, or the annual values changed sign across years).",
      "NEE is negative for net uptake (AmeriFlux sign convention), not flipped. Raw p-values shown;",
      "10 tests, no multiplicity correction in this figure (see diversity_stability_tests.csv for BH FDR)."
    ), width = 160),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0))
  ) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "diversity_vs_function_stability.png"), fig,
       width = 19, height = 6.5, dpi = 200, bg = "white", limitsize = FALSE)

cat("\nSaved CSVs to", data_dir, "and figure to", fig_dir, "\n")
