# ============================================================================
# CompareSpectralVsFieldDiversity.R
#
# Compares spectral diversity metrics against field-measured floristic
# diversity, using PER-TOWER AVERAGED values (not per-site-year) to avoid
# pseudoreplication -- repeated years from the same tower are not
# independent observations, and would otherwise inflate the appearance of
# rank-order agreement. Reports BOTH Spearman's rho and Kendall's tau for
# each comparison (see chat for why both are worth checking together --
# short version: tau is more robust to the tied values common in these
# metrics, e.g. repeated small-integer SSR values).
#
# DESIGNED TO BE RE-RUN AS-IS as AnnualSpectralDiversity.R's output grows --
# always reads both CSVs fresh, re-averages by tower, no changes needed
# between runs.
# ============================================================================

library(tidyverse)
library(patchwork)

# ---- file paths -- EDIT if your paths differ ----
field_csv    <- "./Data/NEON_FieldData/field_diversity_long.csv"
spectral_csv <- "spectral_diversity_by_year.csv"
out_dir      <- "."

# ============================================================================
# 1. Load data fresh
# ============================================================================
field <- read_csv(field_csv, show_col_types = FALSE)
spectral <- read_csv(spectral_csv, show_col_types = FALSE) %>%
  filter(status == "success")

cat("Loaded", nrow(spectral), "completed spectral site-years\n")

# ============================================================================
# 2. Resolve field data to one row per tower-year (plot_scope == "all",
#    prefer peak_flight else average per_bout -- unchanged from before)
# ============================================================================
field_resolved <- field %>%
  filter(plot_scope == "all", status == "ok") %>%
  group_by(tower_id, neon_site, year) %>%
  group_modify(~ {
    peak <- filter(.x, temporal_scope == "peak_flight")
    src <- if (nrow(peak) > 0) peak else .x
    tibble(
      floristic_richness             = mean(src$floristic_richness, na.rm = TRUE),
      floristic_shannon_mean         = mean(src$floristic_shannon_mean, na.rm = TRUE),
      floristic_shannon_gamma        = mean(src$floristic_shannon_gamma, na.rm = TRUE),
      floristic_richness_canopy      = mean(src$floristic_richness_canopy, na.rm = TRUE),
      floristic_shannon_mean_canopy  = mean(src$floristic_shannon_mean_canopy, na.rm = TRUE),
      floristic_shannon_gamma_canopy = mean(src$floristic_shannon_gamma_canopy, na.rm = TRUE)
    )
  }) %>%
  ungroup()

# ============================================================================
# 3. Join spectral + field on tower-year
# ============================================================================
joined_by_year <- spectral %>%
  inner_join(field_resolved, by = c("tower_id", "year", "neon_site"))

cat("Joined:", nrow(joined_by_year), "tower-years with both spectral and field data.\n")

# ============================================================================
# 4. Average per tower across all its available years -- ONE point per site.
#    This is what the correlation tests below actually run on, to avoid
#    treating repeated years from the same tower as independent observations.
# ============================================================================
metric_cols <- c(
  "spectral_species_richness", "cv", "cha", "chv_standardized",
  "shannon_h", "shannon_effective",
  "raoq_ndvi", "raoq_nirv", "raoq_allbands",
  "floristic_richness", "floristic_shannon_mean", "floristic_shannon_gamma",
  "floristic_richness_canopy", "floristic_shannon_mean_canopy", "floristic_shannon_gamma_canopy"
)

# NOTE on chv_standardized: it's a z-score computed across ALL completed
# spectral rows in spectral_diversity_by_year.csv (site-year level, before
# this script averages by tower) -- so it's already standardized upstream,
# and averaging those z-scores per tower here is just aggregating an
# already-comparable quantity, not re-standardizing anything itself.

joined <- joined_by_year %>%
  group_by(tower_id, neon_site) %>%
  summarise(
    n_years_averaged = n(),
    years_included   = paste(sort(year), collapse = ", "),
    across(all_of(metric_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

cat("Averaged to", nrow(joined), "unique towers (from", nrow(joined_by_year), "tower-years).\n")
if (nrow(joined) < 5) {
  cat("\n!! Only", nrow(joined), "towers so far -- both rho and tau will be UNSTABLE\n",
      "   with this few points. Treat as preliminary.\n\n", sep = "")
}

# ============================================================================
# 5. Metric pairs -- SSR only for richness; Shannon + Rao's Q for diversity.
#    Each pair is doubled into an "Unfiltered" and a "Canopy-filtered"
#    variant (field_var + "_canopy", label + " (canopy)") so the two can be
#    compared directly in the improvement summary below (section 7).
# ============================================================================
richness_pairs <- tribble(
  ~field_var,            ~spectral_var,               ~label,
  "floristic_richness",  "spectral_species_richness", "SSR vs. Floristic Richness",
  "floristic_richness",  "cv",                         "CV vs. Floristic Richness",
  "floristic_richness",  "cha",                        "CHA vs. Floristic Richness",
  "floristic_richness",  "chv_standardized",           "CHV (standardized) vs. Floristic Richness"
) %>%
  mutate(field_type = "Unfiltered")

diversity_pairs <- tribble(
  ~field_var,                 ~spectral_var,       ~label,
  "floristic_shannon_mean",   "shannon_h",         "Spectral Shannon H' vs. Floristic Shannon (mean)",
  "floristic_shannon_mean",   "shannon_effective", "Spectral Effective Diversity vs. Floristic Shannon (mean)",
  "floristic_shannon_mean",   "raoq_ndvi",         "Rao's Q (NDVI) vs. Floristic Shannon (mean)",
  "floristic_shannon_mean",   "raoq_nirv",         "Rao's Q (NIRv) vs. Floristic Shannon (mean)",
  "floristic_shannon_mean",   "raoq_allbands",     "Rao's Q (all-bands) vs. Floristic Shannon (mean)",
  "floristic_shannon_gamma",  "shannon_h",         "Spectral Shannon H' vs. Floristic Shannon (gamma)",
  "floristic_shannon_gamma",  "shannon_effective", "Spectral Effective Diversity vs. Floristic Shannon (gamma)",
  "floristic_shannon_gamma",  "raoq_ndvi",         "Rao's Q (NDVI) vs. Floristic Shannon (gamma)",
  "floristic_shannon_gamma",  "raoq_nirv",         "Rao's Q (NIRv) vs. Floristic Shannon (gamma)",
  "floristic_shannon_gamma",  "raoq_allbands",     "Rao's Q (all-bands) vs. Floristic Shannon (gamma)"
) %>%
  mutate(field_type = "Unfiltered")

make_canopy_variant <- function(pairs) {
  pairs %>%
    mutate(field_var  = paste0(field_var, "_canopy"),
           label      = paste0(label, " (canopy)"),
           field_type = "Canopy-filtered")
}

richness_pairs  <- bind_rows(richness_pairs,  make_canopy_variant(richness_pairs))
diversity_pairs <- bind_rows(diversity_pairs, make_canopy_variant(diversity_pairs))

all_pairs <- bind_rows(
  richness_pairs %>% mutate(group = "Richness"),
  diversity_pairs %>% mutate(group = "Diversity")
)

# ============================================================================
# 6. Rank-order tests: Spearman's rho AND Kendall's tau, side by side
# ============================================================================
compute_rank_tests <- function(df, field_var, spectral_var) {
  x <- df[[field_var]]
  y <- df[[spectral_var]]
  ok <- complete.cases(x, y)
  n_ok <- sum(ok)
  if (n_ok < 4) {
    return(tibble(n = n_ok, rho = NA_real_, rho_p = NA_real_,
                  tau = NA_real_, tau_p = NA_real_))
  }
  sp <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  kt <- suppressWarnings(cor.test(x[ok], y[ok], method = "kendall"))
  tibble(n = n_ok,
         rho = unname(sp$estimate), rho_p = sp$p.value,
         tau = unname(kt$estimate), tau_p = kt$p.value)
}

results <- all_pairs %>%
  rowwise() %>%
  mutate(compute_rank_tests(joined, field_var, spectral_var)) %>%
  ungroup()

# ---- FDR (Benjamini-Hochberg) correction, applied SEPARATELY WITHIN each
# group (Richness, Diversity -- each now covering both the Unfiltered and
# Canopy-filtered variant of every pair) -- NOT pooled across both.
# Richness and Diversity are distinct research questions here (different
# field metric, different spectral metric sets, reported/interpreted
# separately throughout this analysis), so they shouldn't share one
# correction family: pooling would let one group's results influence the
# other's adjusted significance just by sharing rank position in a combined
# list (e.g. your strongest Richness result and strongest Diversity result
# could land at adjacent ranks and get tied to the identical adjusted p,
# even though they're unrelated questions) -- an arbitrary coupling that
# per-group correction avoids by keeping each family's threshold
# self-contained. Applied to rho_p and tau_p independently within each
# group, same rationale as before (related but distinct tests, not pooled
# with each other either). ----
results <- results %>%
  group_by(group) %>%
  mutate(
    rho_p_fdr = p.adjust(rho_p, method = "BH"),
    tau_p_fdr = p.adjust(tau_p, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(group, desc(abs(rho)))

cat("\n==== Rank-order tests (per-tower averaged, n =", nrow(joined), "towers) ====\n")
cat("     (rho_p_fdr / tau_p_fdr are Benjamini-Hochberg FDR-corrected WITHIN each\n",
    "      group separately -- Richness's and Diversity's tests (each now including both\n",
    "      Unfiltered and Canopy-filtered variants) are NOT pooled together, since they're\n",
    "      distinct research questions. Use these, not the raw p-values, to judge significance.)\n\n", sep = "")
print(results %>% select(group, label, n, rho, rho_p, rho_p_fdr, tau, tau_p, tau_p_fdr), n = Inf)

write_csv(results, file.path(out_dir, "rank_order_summary_by_tower.csv"))
write_csv(joined, file.path(out_dir, "joined_data_averaged_by_tower.csv"))

# ============================================================================
# 7. Canopy-filtering improvement summary -- does canopy-filtering the field
#    metric improve rank-order agreement with each spectral metric, relative
#    to the unfiltered version of the same comparison? pair_key strips the
#    "_canopy" suffix from field_var so "floristic_richness" and
#    "floristic_richness_canopy" (and likewise the Shannon mean/gamma
#    variants) are matched up as the same underlying comparison.
# ============================================================================
improvement_summary <- results %>%
  mutate(pair_key = str_remove(field_var, "_canopy$"),
         type_key = if_else(field_type == "Unfiltered", "unfiltered", "canopy")) %>%
  select(group, pair_key, spectral_var, type_key, rho, tau) %>%
  pivot_wider(names_from = type_key, values_from = c(rho, tau),
              names_glue = "{type_key}_{.value}") %>%
  mutate(rho_delta = canopy_rho - unfiltered_rho,
         tau_delta = canopy_tau - unfiltered_tau) %>%
  arrange(group, desc(rho_delta))

cat("\n==== Canopy-filtering improvement summary (per group, sorted by rho_delta desc) ====\n")
cat("     (rho_delta / tau_delta = canopy-filtered stat minus unfiltered stat for the SAME\n",
    "      field metric vs. spectral metric pair. Positive means canopy-filtering improved\n",
    "      rank-order agreement with that spectral metric; negative means it hurt.)\n\n", sep = "")
print(improvement_summary, n = Inf)

write_csv(improvement_summary, file.path(out_dir, "canopy_filtering_improvement_summary.csv"))

# ============================================================================
# 8. Plots -- one clean, larger panel per pair, stacked
# ============================================================================
make_panel <- function(field_var, spectral_var, label) {
  d <- joined %>%
    select(tower_id, neon_site, x = all_of(field_var), y = all_of(spectral_var)) %>%
    filter(!is.na(x), !is.na(y))
  
  r <- results %>% filter(field_var == !!field_var, spectral_var == !!spectral_var, label == !!label)
  stat_lab <- sprintf("Spearman rho = %.2f (p = %.3f, FDR p = %.3f)   |   Kendall tau = %.2f (p = %.3f, FDR p = %.3f)   |   n = %d",
                      r$rho, r$rho_p, r$rho_p_fdr, r$tau, r$tau_p, r$tau_p_fdr, r$n)
  
  ggplot(d, aes(x = x, y = y)) +
    geom_point(size = 3, alpha = 0.75, color = "#2c6e91") +
    geom_smooth(method = "lm", se = TRUE, color = "#c0392b", linewidth = 0.6,
                formula = y ~ x, na.rm = TRUE) +
    labs(title = label, subtitle = stat_lab, x = field_var, y = spectral_var) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, color = "grey30"))
}

richness_plots <- richness_pairs %>%
  select(field_var, spectral_var, label) %>%
  pmap(make_panel)

diversity_plots <- diversity_pairs %>%
  select(field_var, spectral_var, label) %>%
  pmap(make_panel)

richness_stack <- wrap_plots(richness_plots, ncol = 1) +
  plot_annotation(title = "Field Richness vs. Spectral Richness-Type Metrics (per-tower average)")
ggsave(file.path(out_dir, "richness_by_tower.png"), richness_stack,
       width = 7, height = 5 * length(richness_plots), dpi = 150, limitsize = FALSE)

diversity_stack <- wrap_plots(diversity_plots, ncol = 1) +
  plot_annotation(title = "Field Shannon Diversity vs. Spectral Diversity Metrics (per-tower average)")
ggsave(file.path(out_dir, "diversity_by_tower.png"), diversity_stack,
       width = 7, height = 5 * length(diversity_plots), dpi = 150, limitsize = FALSE)

cat("\nSaved:\n",
    "  - rank_order_summary_by_tower.csv\n",
    "  - joined_data_averaged_by_tower.csv\n",
    "  - richness_by_tower.png\n",
    "  - diversity_by_tower.png (tall, one panel per pair)\n", sep = "")

# ============================================================================
# 9. Table figure: ALL tests, columns grouped by statistic (all rho columns
#    together, then all tau columns together, raw-before-FDR within each),
#    with FDR-significant cells highlighted independently for rho and tau
#    (NOT requiring both to agree -- each gets its own flag). A small gap
#    separates the rho block from the tau block, matching the existing gap
#    between the Richness/Diversity row groups. Built in ggplot2 (no new
#    dependency) rather than a package like gt, since gt's image export
#    needs webshot2/Chromium -- an extra system dependency this script's
#    environment shouldn't need to rely on.
# ============================================================================

fdr_alpha <- 0.05   # significance threshold applied to the FDR-corrected
# p-values only -- change here if you want a stricter/
# looser cutoff; does not affect the raw p columns,
# which are shown for reference only and never flagged.

# ---- build the display table: one row per test, columns grouped by
# statistic (all rho columns together, then all tau columns together),
# raw-before-FDR within each group. n is dropped as a column here since
# it's constant across all tests (same per-tower joined dataset) -- shown
# once in the subtitle instead. ----
n_values <- unique(results$n)
n_note <- if (length(n_values) == 1) {
  sprintf("n = %d towers", n_values)
} else {
  sprintf("n ranges %d-%d towers across tests (varies by metric availability)",
          min(n_values), max(n_values))
}

table_df <- results %>%
  arrange(group, desc(abs(rho))) %>%
  mutate(
    row_label = label,
    rho_fmt      = sprintf("%.2f", rho),
    rho_p_fmt    = sprintf("%.3f", rho_p),
    rho_p_fdr_fmt = sprintf("%.3f", rho_p_fdr),
    tau_fmt      = sprintf("%.2f", tau),
    tau_p_fmt    = sprintf("%.3f", tau_p),
    tau_p_fdr_fmt = sprintf("%.3f", tau_p_fdr),
    rho_fdr_sig  = rho_p_fdr < fdr_alpha,
    tau_fdr_sig  = tau_p_fdr < fdr_alpha
  )

# row order top-to-bottom (ggplot y-axis plots bottom-to-top by default, so
# reverse for a natural reading order with the first row at the top)
row_order <- table_df$row_label
table_df <- table_df %>% mutate(row_label = factor(row_label, levels = rev(row_order)))

# rho columns grouped together (raw, raw p, FDR p), then tau columns
# grouped together (raw, raw p, FDR p) -- n column removed
col_order <- c("rho", "rho_p", "rho_p_fdr", "tau", "tau_p", "tau_p_fdr")
col_labels <- c("rho", "rho p", "rho p (FDR)", "tau", "tau p", "tau p (FDR)")

# long format: one row per (test x column) cell, with the text to display,
# which stat family (rho/tau) it belongs to (used only to create the visual
# gap below), and whether that specific cell should be flagged as
# FDR-significant
cell_long <- table_df %>%
  transmute(
    row_label, group,
    rho = rho_fmt, rho_p = rho_p_fmt, rho_p_fdr = rho_p_fdr_fmt,
    tau = tau_fmt, tau_p = tau_p_fmt, tau_p_fdr = tau_p_fdr_fmt,
    rho_fdr_sig, tau_fdr_sig
  ) %>%
  pivot_longer(cols = all_of(col_order), names_to = "column", values_to = "value") %>%
  mutate(
    stat_family = if_else(str_starts(column, "rho"), "rho", "tau"),  # which side of the gap this column belongs to
    column = factor(column, levels = col_order, labels = col_labels),
    is_fdr_col = column %in% c("rho p (FDR)", "tau p (FDR)"),
    is_sig = case_when(
      column == "rho p (FDR)" ~ rho_fdr_sig,
      column == "tau p (FDR)" ~ tau_fdr_sig,
      TRUE ~ FALSE
    )
  )

# ---- render as a ggplot table: tiles for background (light highlight only
# on FDR-significant cells), text for values, group labels down the left
# (rotated vertical, like a standard y-axis title). Columns are also
# faceted by stat_family (rho vs. tau) purely to get a visual gap between
# the two blocks via panel.spacing.x, matching the existing gap between the
# Richness/Diversity row groups -- the column-facet strip itself is hidden
# since the column labels underneath already say "rho"/"tau". ----
table_fig <- ggplot(cell_long, aes(x = column, y = row_label)) +
  geom_tile(aes(fill = is_sig), color = "grey85", linewidth = 0.4) +
  geom_text(aes(label = value, fontface = ifelse(is_sig, "bold", "plain")), size = 5) +
  scale_fill_manual(values = c(`TRUE` = "#d9f2d0", `FALSE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  facet_grid(rows = vars(group), cols = vars(stat_family),
             scales = "free", space = "free", switch = "y") +
  labs(x = NULL, y = NULL,
       title = "Rank-Order Field vs. Spectral Diversity",
       subtitle = sprintf("Highlighted cells: FDR-corrected p < %.2f   |   %s",
                          fdr_alpha, n_note)) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.8, "cm"),   # gap between the rho and tau column blocks
    panel.spacing.y = unit(0.3, "cm"),   # gap between the Richness/Diversity row blocks
    axis.text.y = element_text(hjust = 1, size = 13),
    axis.text.x = element_text(face = "bold", size = 13),
    strip.placement = "outside",
    strip.text.x = element_blank(),       # hide the rho/tau column-facet strip -- redundant with column labels
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 14),
    plot.title = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )

ggsave(file.path(out_dir, "rank_order_table_figure.png"), table_fig,
       width = 9.5, height = 0.5 * nrow(table_df) + 2, dpi = 150)

cat("  - rank_order_table_figure.png (all tests, FDR-significant cells highlighted)\n")