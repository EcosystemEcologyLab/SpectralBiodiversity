#!/usr/bin/env Rscript
# Per-site temporal variability (CV%) across 7 spectral diversity metrics,
# by ecosystem (IGBP) type, across NEON sites.
# Usage: Rscript Code/DataAnalysis/RaoQYearVariability.R [csv] [data_dir] [fig_dir]
#   csv      default "./Data/spectral_diversity_by_year.csv"
#   data_dir default "./Data"    (spectral_site_temporal_cv.csv, raoq_year_anova_by_ecosystem.csv)
#   fig_dir  default "./Figures" (spectral_cv_by_method.png, raoq_cv_by_ecosystem.png, raoq_cv_by_site.png)
# Also runnable line-by-line in RStudio with the project root as the working
# directory -- every default path below is relative to that root.

required_pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

args     <- commandArgs(trailingOnly = TRUE)
csv      <- if (length(args) >= 1) args[1] else "./Data/spectral_diversity_by_year.csv"
data_dir <- if (length(args) >= 2) args[2] else "./Data"
fig_dir  <- if (length(args) >= 3) args[3] else "./Figures"

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

if (!file.exists(csv)) stop("Input csv not found: ", csv)

# ---- metric mapping (col -> display label), fixed order used by Figure 1 ----
# Evidence for these exact column names: Code/AnnualSpectralDiversity.R's
# `results` tibble definition (section 8, ~line 966) and the row it writes to
# spectral_diversity_by_year.csv (~line 1040): tower_id, neon_site, year, cv,
# chv, chv_standardized, cha, spectral_species_richness, shannon_h,
# shannon_effective, ssr_mean_winning_c, ssr_max_winning_c,
# n_reps_gap_check_triggered, raoq_ndvi, raoq_nirv, raoq_allbands, status.
# `spectral_species_richness` is the SSR metric itself; ssr_mean_winning_c /
# ssr_max_winning_c / n_reps_gap_check_triggered are per-rep diagnostics, NOT
# used here. CHV uses the raw `chv` column, not `chv_standardized` (that
# z-score is computed across ALL site-years for cross-metric comparison in
# CompareSpectralVsFieldDiversity.R and is not the natural per-site-year CHV
# value for a within-site year-to-year variability question).
metric_map <- tibble::tribble(
  ~col,                         ~label,
  "raoq_nirv",                  "Rao's Q (NIRv)",
  "raoq_allbands",              "Rao's Q (all bands)",
  "shannon_h",                  "Shannon H'",
  "cha",                        "CHA",
  "chv",                        "CHV",
  "cv",                         "CV",
  "spectral_species_richness",  "SSR"
)

id_cols       <- c("tower_id", "neon_site", "year")
required_cols <- c(id_cols, metric_map$col)

# ---- ecosystem bins (IGBP -> bin), per spec. Any IGBP code present in the
# site file that isn't in this mapping stops the run rather than silently
# dropping that site (see Section 3 below). ----
bin_map    <- c(DBF = "Forest", EBF = "Forest", ENF = "Forest", MF = "Forest",
                 GRA = "Grassland",
                 OSH = "Shrub/Savanna", SAV = "Shrub/Savanna",
                 CRO = "Cropland/Mosaic", CVM = "Cropland/Mosaic",
                 WET = "Wetland")
bin_levels <- c("Forest", "Grassland", "Shrub/Savanna", "Cropland/Mosaic", "Wetland")
bin_colors <- c(Forest = "#2a9d8f", Grassland = "#e9c46a",
                 `Shrub/Savanna` = "#f4a261", `Cropland/Mosaic` = "#8ab17d",
                 Wetland = "#277da1")

cv_mean_eps <- 1e-8  # |mean| below this is treated as "zero" for CV% (guards div-by-~0)

fmt_p <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
stars_from_p <- function(p) dplyr::case_when(
  is.na(p)   ~ "n/a",
  p < 0.001  ~ "***",
  p < 0.01   ~ "**",
  p < 0.05   ~ "*",
  TRUE       ~ "ns"
)

# ============================================================================
# 1. Load + validate the spectral diversity CSV
# ============================================================================
dat <- read_csv(csv, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "),
       "\nColumns present in file: ", paste(names(dat), collapse = ", "))
}

dat <- dat %>% mutate(year = as.integer(year))

status_col <- intersect(c("status", "ok", "error"), names(dat))
if (length(status_col) > 0) {
  status_col <- status_col[1]
  cat("\n=== Status column found:", status_col, "===\n")
  print(table(dat[[status_col]], useNA = "ifany"))
  if (status_col %in% c("status", "ok")) {
    success_vals <- c("success", "ok", "TRUE", TRUE, 1)
    dat <- dat %>% filter(.data[[status_col]] %in% success_vals)
  } else {
    dat <- dat %>% filter(is.na(.data[[status_col]]) | .data[[status_col]] == "")
  }
  cat("Filtered to successful rows using column '", status_col, "'\n", sep = "")
} else {
  cat("\nNo status/success column found (checked: status, ok, error). Proceeding without filtering.\n")
}

dups <- dat %>% count(neon_site, year) %>% filter(n > 1)
if (nrow(dups) > 0) {
  offending <- dat %>% semi_join(dups, by = c("neon_site", "year")) %>%
    arrange(neon_site, year)
  print(offending)
  stop("Duplicate neon_site x year rows found (per-site CV and the ANOVA assume ",
       "one value per site-year). See offending rows printed above.")
}

cat("\n=== Post-filter summary ===\n")
cat("n site-years:", nrow(dat), "\n")
cat("n sites:", n_distinct(dat$neon_site), "\n")
cat("year range:", paste(range(dat$year, na.rm = TRUE), collapse = " - "), "\n\n")

# ============================================================================
# 2. Ecosystem mapping -- Veg Type from NEONsites.csv or NEONsites_Footprints.csv
# ============================================================================
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
cat("Using ecosystem/Veg Type file:", sites_path, "(column '", veg_col, "')\n", sep = "")

sites_raw <- read.csv(sites_path, fileEncoding = "UTF-8-BOM")

observed_codes <- unique(sites_raw[[veg_col]])
unmapped_codes <- setdiff(observed_codes, names(bin_map))
if (length(unmapped_codes) > 0) {
  stop("'", veg_col, "' in ", sites_path, " contains IGBP code(s) not covered by the ",
       "Forest/Grassland/Shrub-Savanna/Cropland-Mosaic/Wetland bin mapping: ",
       paste(unmapped_codes, collapse = ", "),
       ". Extend the mapping deliberately rather than dropping these sites silently.")
}

# tower_id IS Site.ID (AmeriFlux-style code), same join key
# Code/CanopyDeltaByEcosystem.R uses -- NOT the 4-letter neon_site code.
sites_lookup <- sites_raw %>%
  transmute(tower_id = Site.ID,
            ecosystem_bin = bin_map[.data[[veg_col]]]) %>%
  distinct(tower_id, .keep_all = TRUE)

dat <- dat %>% left_join(sites_lookup, by = "tower_id")

unmatched <- dat %>% filter(is.na(ecosystem_bin)) %>% distinct(tower_id, neon_site)
if (nrow(unmatched) > 0) {
  cat("\n!! ", nrow(unmatched), " site(s) in the spectral data have no ecosystem match in ",
      sites_path, " (excluded from the ecosystem-binned figures, kept in Figure 1):\n", sep = "")
  print(as.data.frame(unmatched))
} else {
  cat("\nAll sites in the spectral data matched an ecosystem bin.\n")
}

# ============================================================================
# 3. Per-site temporal CV for all 7 metrics (sites with >= 3 years, non-NA)
# ============================================================================
dat_long <- dat %>%
  select(tower_id, neon_site, ecosystem_bin, year, all_of(metric_map$col)) %>%
  pivot_longer(cols = all_of(metric_map$col), names_to = "metric", values_to = "value")

cv_tbl <- dat_long %>%
  filter(!is.na(value)) %>%
  group_by(metric, site = neon_site, ecosystem_bin) %>%
  filter(n() >= 3) %>%
  summarise(n_years = n(), mean_val = mean(value), sd_val = sd(value), .groups = "drop") %>%
  mutate(
    cv_raw = 100 * sd_val / mean_val,
    cv_guarded = !is.finite(cv_raw) | abs(mean_val) < cv_mean_eps,
    cv_pct = ifelse(cv_guarded, NA_real_, cv_raw),
    metric_label = factor(metric_map$label[match(metric, metric_map$col)], levels = metric_map$label)
  )

guard_counts <- cv_tbl %>% filter(cv_guarded) %>% count(metric, name = "n_guarded")
if (nrow(guard_counts) > 0) {
  cat("\nCV guard (mean ~ 0) triggered for:\n")
  for (i in seq_len(nrow(guard_counts))) {
    cat("  ", guard_counts$metric[i], ": ", guard_counts$n_guarded[i], " site(s) set to NA\n", sep = "")
  }
} else {
  cat("\nCV guard (mean ~ 0) did not trigger for any metric/site.\n")
}

write_csv(
  cv_tbl %>% transmute(metric, site, ecosystem_bin, n_years, mean = mean_val, sd = sd_val, cv_pct),
  file.path(data_dir, "spectral_site_temporal_cv.csv")
)

# ============================================================================
# 4. FIGURE 1 -- per-site temporal CV by method, all 7 metrics
# ============================================================================
cat("\n=== Per-site temporal CV by method ===\n")
method_stats <- cv_tbl %>%
  group_by(metric_label) %>%
  summarise(n_sites = n_distinct(site), min = min(cv_pct, na.rm = TRUE),
            median = median(cv_pct, na.rm = TRUE), max = max(cv_pct, na.rm = TRUE), .groups = "drop")
print(method_stats)

cv_pos1 <- cv_tbl$cv_pct[!is.na(cv_tbl$cv_pct) & cv_tbl$cv_pct > 0]
use_log1 <- length(cv_pos1) > 0 && (max(cv_pos1) / min(cv_pos1)) > 50
cat("\nFigure 1 y-axis: ", if (use_log1) "log10 (max/min positive CV > 50x)" else "linear", "\n", sep = "")

cv_tbl <- cv_tbl %>%
  group_by(metric_label) %>%
  mutate(q1 = quantile(cv_pct, 0.25, na.rm = TRUE), q3 = quantile(cv_pct, 0.75, na.rm = TRUE),
         iqr = q3 - q1, lower = q1 - 1.5 * iqr, upper = q3 + 1.5 * iqr,
         is_outlier = !is.na(cv_pct) & (cv_pct < lower | cv_pct > upper)) %>%
  ungroup()
outliers1 <- cv_tbl %>% filter(is_outlier)

p1 <- ggplot(cv_tbl, aes(x = metric_label, y = cv_pct)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey30", na.rm = TRUE) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2, colour = "#2c6e91", na.rm = TRUE) +
  geom_text_repel(data = outliers1, aes(label = site), size = 2.5, max.overlaps = Inf,
                   min.segment.length = 0, segment.size = 0.3) +
  { if (use_log1) scale_y_log10() } +
  labs(x = NULL,
       y = if (use_log1) "Per-site temporal CV (%) [log10 scale]" else "Per-site temporal CV (%)",
       title = "Per-site temporal variability by method") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(fig_dir, "spectral_cv_by_method.png"), p1, width = 9, height = 6, dpi = 200)

# ============================================================================
# 5. Year-effect ANOVA per Rao's Q metric x ecosystem bin (site-year level)
# ============================================================================
raoq_cols <- c("raoq_nirv", "raoq_allbands")

run_ecosystem_anova <- function(metric_col, bin) {
  d <- dat %>%
    filter(ecosystem_bin == bin) %>%
    select(site = neon_site, year, y = all_of(metric_col)) %>%
    filter(!is.na(y)) %>%
    group_by(site) %>% filter(n() >= 2) %>% ungroup()

  n_sites <- n_distinct(d$site); n_years <- n_distinct(d$year); n_obs <- nrow(d)
  if (n_sites < 3 || n_years < 2) {
    return(tibble(n_sites = n_sites, n_years = n_years, n_obs = n_obs,
                   F = NA_real_, df1 = NA_real_, df2 = NA_real_, p_raw = NA_real_, status = "n/a"))
  }

  fit_res <- tryCatch({
    an <- anova(lm(y ~ site + factor(year), data = d %>% mutate(site = factor(site))))
    tibble(n_sites = n_sites, n_years = n_years, n_obs = n_obs,
           F = an["factor(year)", "F value"], df1 = an["factor(year)", "Df"],
           df2 = an["Residuals", "Df"], p_raw = an["factor(year)", "Pr(>F)"], status = "ok")
  }, error = function(e) {
    message(sprintf("[%s / %s] ANOVA failed: %s", metric_col, bin, conditionMessage(e)))
    tibble(n_sites = n_sites, n_years = n_years, n_obs = n_obs,
           F = NA_real_, df1 = NA_real_, df2 = NA_real_, p_raw = NA_real_, status = "n/a")
  })
  fit_res
}

anova_res <- tidyr::expand_grid(metric = raoq_cols, bin = bin_levels) %>%
  rowwise() %>%
  mutate(run_ecosystem_anova(metric, bin)) %>%
  ungroup() %>%
  group_by(metric) %>%
  mutate(p_fdr = {
    out <- rep(NA_real_, n())
    ok_idx <- which(status == "ok")
    if (length(ok_idx) > 0) out[ok_idx] <- p.adjust(p_raw[ok_idx], method = "BH")
    out
  }) %>%
  ungroup() %>%
  mutate(stars = stars_from_p(p_fdr),
         metric_label = factor(metric_map$label[match(metric, metric_map$col)], levels = metric_map$label))

cat("\n=== Year-effect ANOVA by ecosystem (site-adjusted, BH FDR within metric) ===\n")
print(anova_res %>% select(metric, bin, n_sites, n_years, n_obs, F, df1, df2, p_raw, p_fdr, stars, status))

write_csv(
  anova_res %>% transmute(metric, bin, n_sites, n_years, n_obs, F, df1, df2, p_raw, p_fdr, stars, status),
  file.path(data_dir, "raoq_year_anova_by_ecosystem.csv")
)

# ============================================================================
# 6. FIGURE 2 -- ecosystem boxplots (Rao's Q NIRv / all bands), per-site CV
# ============================================================================
cv_raoq <- cv_tbl %>%
  filter(metric %in% raoq_cols, !is.na(ecosystem_bin)) %>%
  mutate(ecosystem_bin = factor(ecosystem_bin, levels = bin_levels))

box_n <- cv_raoq %>% group_by(metric_label, ecosystem_bin) %>% summarise(n = n(), .groups = "drop")
star_lab <- anova_res %>% select(metric_label, bin, stars) %>%
  mutate(ecosystem_bin = factor(bin, levels = bin_levels))

p2 <- ggplot(cv_raoq, aes(x = ecosystem_bin, y = cv_pct, fill = ecosystem_bin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.1, na.rm = TRUE) +
  geom_text(data = box_n, aes(x = ecosystem_bin, y = -Inf, label = paste0("n = ", n), fill = NULL),
             vjust = -0.4, size = 2.8, colour = "grey30") +
  geom_text(data = star_lab, aes(x = ecosystem_bin, y = Inf, label = stars, fill = NULL),
             vjust = 1.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = bin_colors, guide = "none") +
  facet_wrap(~ metric_label, scales = "free_y") +
  labs(x = NULL, y = "Per-site temporal CV (%)",
       title = "Per-site temporal variability of Rao's Q, by ecosystem",
       caption = "Stars = site-adjusted year effect (two-way ANOVA, BH FDR within metric); boxes = per-site temporal CV.") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.caption = element_text(size = 7, colour = "grey40"))

ggsave(file.path(fig_dir, "raoq_cv_by_ecosystem.png"), p2, width = 11, height = 6, dpi = 200)

# ============================================================================
# 7. FIGURE 3 -- per-site CV bars, two stacked panels (one per Rao's Q metric),
#    sites ordered by ecosystem bin then descending CV within bin (ordering
#    computed independently per panel/metric, so each panel gets its own
#    site order -- not achievable with a single shared facet axis, hence
#    patchwork rather than facet_wrap here).
# ============================================================================
bar_plot_for_metric <- function(mcol) {
  d <- cv_raoq %>% filter(metric == mcol) %>%
    arrange(ecosystem_bin, desc(cv_pct)) %>%
    mutate(site = factor(site, levels = unique(site)))
  mlabel <- metric_map$label[metric_map$col == mcol]

  ggplot(d, aes(x = site, y = cv_pct, fill = ecosystem_bin)) +
    geom_col(na.rm = TRUE) +
    scale_fill_manual(values = bin_colors, drop = FALSE) +
    labs(x = NULL, y = "Per-site temporal CV (%)", fill = "Ecosystem", title = mlabel) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
          plot.title = element_text(face = "bold"))
}

p3 <- (bar_plot_for_metric("raoq_nirv") / bar_plot_for_metric("raoq_allbands")) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Per-site temporal CV of Rao's Q, ordered by ecosystem bin then descending CV")

ggsave(file.path(fig_dir, "raoq_cv_by_site.png"), p3, width = max(10, 0.22 * n_distinct(cv_raoq$site)),
       height = 9, dpi = 200, limitsize = FALSE)

cat("\nSaved CSVs to", data_dir, "and figures to", fig_dir, "\n")
