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
# TWO INDEPENDENT PASSES: this script runs the SAME 14 tests (4 Richness +
# 10 Diversity) twice over the SAME tower population -- once against the
# unfiltered field metrics ("all observations" pass) and once against the
# canopy-filtered field metrics ("filtered observations" pass) -- each with
# its own independent Benjamini-Hochberg FDR correction within its own
# Richness/Diversity groups. The two passes are NEVER pooled into one FDR
# family; see section 6 for why. This replaces an earlier single-pass
# design that pooled both variants into one 8-test Richness family and one
# 20-test Diversity family before correcting.
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
data_dir     <- "./Data"      # CSV outputs
fig_dir      <- "./Figures"   # PNG outputs
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

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
#    Both passes below run over this SAME per-tower table (same population,
#    same n) -- only which field columns get tested differs between passes.
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
#    These are the BASE (unfiltered) pairs -- 4 Richness + 10 Diversity, 14
#    total, unchanged from the original single-pass design.
#    make_canopy_variant() re-points field_var at the "_canopy" column for
#    the exact same spectral_var/label pairing -- this is how the filtered
#    pass below reuses these definitions instead of hand-writing a second
#    set of tribbles.
# ============================================================================
richness_pairs <- tribble(
  ~field_var,            ~spectral_var,               ~label,
  "floristic_richness",  "spectral_species_richness", "SSR vs. Floristic Richness",
  "floristic_richness",  "cv",                         "CV vs. Floristic Richness",
  "floristic_richness",  "cha",                        "CHA vs. Floristic Richness",
  "floristic_richness",  "chv_standardized",           "CHV (standardized) vs. Floristic Richness"
)

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
)

make_canopy_variant <- function(pairs) {
  pairs %>%
    mutate(field_var = paste0(field_var, "_canopy"),
           label     = paste0(label, " (canopy)"))
}

# ============================================================================
# 6. Rank-order tests: Spearman's rho AND Kendall's tau, side by side.
#    compute_rank_tests() is UNCHANGED from the original single-pass script
#    -- same per-test Spearman/Kendall logic, same n<4 guard.
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

# ---- run_rank_order_pass(): one full pass (14 tests: 4 Richness + 10
# Diversity) against whichever richness/diversity pairs it's handed, with
# FDR correction applied SEPARATELY within Richness and Diversity (same
# rationale as before -- they're distinct research questions with different
# field/spectral metric sets, so pooling would let one group's results
# influence the other's adjusted significance just by sharing rank position
# in a combined list). Called TWICE below: once on the unfiltered pairs,
# once on the canopy-filtered pairs -- each call's FDR correction sees ONLY
# its own 4-test Richness family and its own 10-test Diversity family, never
# the other pass's tests. This is what makes the two passes independent: a
# test's adjusted p-value here depends only on the other 3 (Richness) or 9
# (Diversity) tests in the SAME pass, not on all 7/19 sibling tests across
# both passes the way the old pooled 8-/20-test families did. ----
run_rank_order_pass <- function(richness, diversity, joined) {
  bind_rows(
    richness %>% mutate(group = "Richness"),
    diversity %>% mutate(group = "Diversity")
  ) %>%
    rowwise() %>%
    mutate(compute_rank_tests(joined, field_var, spectral_var)) %>%
    ungroup() %>%
    group_by(group) %>%
    mutate(
      rho_p_fdr = p.adjust(rho_p, method = "BH"),
      tau_p_fdr = p.adjust(tau_p, method = "BH")
    ) %>%
    ungroup() %>%
    arrange(group, desc(abs(rho)))
}

results_all      <- run_rank_order_pass(richness_pairs, diversity_pairs, joined)
results_filtered <- run_rank_order_pass(make_canopy_variant(richness_pairs),
                                         make_canopy_variant(diversity_pairs), joined)

cat("\n==== PASS 1: All observations (unfiltered field metrics, n =", nrow(joined), "towers) ====\n")
cat("     (rho_p_fdr / tau_p_fdr are BH-corrected within this pass's own Richness (n=4)\n",
    "      and Diversity (n=10) families only -- NOT pooled with the filtered pass below.)\n\n", sep = "")
print(results_all %>% select(group, label, n, rho, rho_p, rho_p_fdr, tau, tau_p, tau_p_fdr), n = Inf)

cat("\n==== PASS 2: Filtered observations (canopy-filtered field metrics, n =", nrow(joined), "towers) ====\n")
cat("     (rho_p_fdr / tau_p_fdr are BH-corrected within THIS pass's own Richness (n=4)\n",
    "      and Diversity (n=10) families only -- independent of pass 1 above.)\n\n", sep = "")
print(results_filtered %>% select(group, label, n, rho, rho_p, rho_p_fdr, tau, tau_p, tau_p_fdr), n = Inf)

write_csv(results_all,      file.path(data_dir, "rank_order_summary_all_observations.csv"))
write_csv(results_filtered, file.path(data_dir, "rank_order_summary_filtered_observations.csv"))
write_csv(joined,           file.path(data_dir, "joined_data_averaged_by_tower.csv"))

# ============================================================================
# 7. Canopy-filtering improvement summary -- does canopy-filtering the field
#    metric improve rank-order agreement with each spectral metric, relative
#    to the unfiltered version of the same comparison? UNCHANGED pivot/delta
#    logic from the original single-pass script -- only its INPUT changed,
#    from one pooled 28-row results table to bind_rows(results_all,
#    results_filtered) (28 rows total: 14 + 14, same as before). rho/tau
#    themselves are computed identically either way (compute_rank_tests() is
#    untouched), so rho_delta/tau_delta here are numerically IDENTICAL to
#    what the old pooled approach would have produced -- only the FDR
#    columns (not used in this delta) differ between the two designs.
#    pair_key strips the "_canopy" suffix from field_var so
#    "floristic_richness" and "floristic_richness_canopy" (and likewise the
#    Shannon mean/gamma variants) are matched up as the same underlying
#    comparison.
# ============================================================================
improvement_summary <- bind_rows(
    results_all      %>% mutate(field_type = "Unfiltered"),
    results_filtered %>% mutate(field_type = "Canopy-filtered")
  ) %>%
  mutate(pair_key = str_remove(field_var, "_canopy$"),
         type_key = if_else(field_type == "Unfiltered", "unfiltered", "canopy")) %>%
  select(group, pair_key, spectral_var, type_key, rho, tau) %>%
  pivot_wider(names_from = type_key, values_from = c(rho, tau),
              names_glue = "{type_key}_{.value}") %>%
  mutate(rho_delta = canopy_rho - unfiltered_rho,
         tau_delta = canopy_tau - unfiltered_tau) %>%
  arrange(group, desc(rho_delta))

cat("\n==== Canopy-filtering improvement summary (per group, sorted by rho_delta desc) ====\n")
cat("     (rho_delta / tau_delta = filtered-pass stat minus all-observations-pass stat for the\n",
    "      SAME field metric vs. spectral metric pair. Positive means canopy-filtering improved\n",
    "      rank-order agreement with that spectral metric; negative means it hurt.)\n\n", sep = "")
print(improvement_summary, n = Inf)

write_csv(improvement_summary, file.path(data_dir, "canopy_filtering_improvement_summary.csv"))

# ============================================================================
# 8. Scatter panels -- one clean, larger panel per pair, stacked. Unchanged
#    from the original script except: (a) the per-pair stat label now pulls
#    from results_combined (results_all + results_filtered stacked back
#    together purely as a lookup table for this plot's subtitle text, not a
#    shared FDR family -- each row still carries the FDR value it got from
#    its own pass), and (b) output paths now point at fig_dir.
# ============================================================================
results_combined <- bind_rows(results_all, results_filtered)

make_panel <- function(field_var, spectral_var, label) {
  d <- joined %>%
    select(tower_id, neon_site, x = all_of(field_var), y = all_of(spectral_var)) %>%
    filter(!is.na(x), !is.na(y))

  r <- results_combined %>% filter(field_var == !!field_var, spectral_var == !!spectral_var, label == !!label)
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

richness_pairs_both  <- bind_rows(richness_pairs,  make_canopy_variant(richness_pairs))
diversity_pairs_both <- bind_rows(diversity_pairs, make_canopy_variant(diversity_pairs))

richness_plots <- richness_pairs_both %>%
  select(field_var, spectral_var, label) %>%
  pmap(make_panel)

diversity_plots <- diversity_pairs_both %>%
  select(field_var, spectral_var, label) %>%
  pmap(make_panel)

richness_stack <- wrap_plots(richness_plots, ncol = 1) +
  plot_annotation(title = "Field Richness vs. Spectral Richness-Type Metrics (per-tower average)")
ggsave(file.path(fig_dir, "richness_by_tower.png"), richness_stack,
       width = 7, height = 5 * length(richness_plots), dpi = 150, limitsize = FALSE, bg = "white")

diversity_stack <- wrap_plots(diversity_plots, ncol = 1) +
  plot_annotation(title = "Field Shannon Diversity vs. Spectral Diversity Metrics (per-tower average)")
ggsave(file.path(fig_dir, "diversity_by_tower.png"), diversity_stack,
       width = 7, height = 5 * length(diversity_plots), dpi = 150, limitsize = FALSE, bg = "white")

cat("\nSaved:\n",
    "  - ", data_dir, "/rank_order_summary_all_observations.csv\n",
    "  - ", data_dir, "/rank_order_summary_filtered_observations.csv\n",
    "  - ", data_dir, "/joined_data_averaged_by_tower.csv\n",
    "  - ", fig_dir, "/richness_by_tower.png\n",
    "  - ", fig_dir, "/diversity_by_tower.png (tall, one panel per pair)\n", sep = "")

# ============================================================================
# 9. Table figure -- ALL tests in one pass, columns grouped by statistic
#    (all rho columns together, then all tau columns together, raw-before-
#    FDR within each), FDR-significant cells highlighted independently for
#    rho and tau. Refactored into make_table_figure(results_df, title_txt)
#    so the SAME rendering code produces both passes' table figures rather
#    than duplicating the ggplot-building logic twice. Built in ggplot2 (no
#    new dependency) rather than a package like gt, since gt's image export
#    needs webshot2/Chromium -- an extra system dependency this script's
#    environment shouldn't need to rely on.
# ============================================================================
make_table_figure <- function(results_df, title_txt, fdr_alpha = 0.05) {
  # n is dropped as a column here since it's constant (or near-constant)
  # across tests within a pass -- shown once in the subtitle instead.
  n_values <- unique(results_df$n)
  n_note <- if (length(n_values) == 1) {
    sprintf("n = %d towers", n_values)
  } else {
    sprintf("n ranges %d-%d towers across tests (varies by metric availability)",
            min(n_values), max(n_values))
  }

  table_df <- results_df %>%
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

  # row order top-to-bottom (ggplot y-axis plots bottom-to-top by default,
  # so reverse for a natural reading order with the first row at the top)
  row_order <- table_df$row_label
  table_df <- table_df %>% mutate(row_label = factor(row_label, levels = rev(row_order)))

  # rho columns grouped together (raw, raw p, FDR p), then tau columns
  # grouped together (raw, raw p, FDR p)
  col_order <- c("rho", "rho_p", "rho_p_fdr", "tau", "tau_p", "tau_p_fdr")
  col_labels <- c("rho", "rho p", "rho p (FDR)", "tau", "tau p", "tau p (FDR)")

  cell_long <- table_df %>%
    transmute(
      row_label, group,
      rho = rho_fmt, rho_p = rho_p_fmt, rho_p_fdr = rho_p_fdr_fmt,
      tau = tau_fmt, tau_p = tau_p_fmt, tau_p_fdr = tau_p_fdr_fmt,
      rho_fdr_sig, tau_fdr_sig
    ) %>%
    pivot_longer(cols = all_of(col_order), names_to = "column", values_to = "value") %>%
    mutate(
      stat_family = if_else(str_starts(column, "rho"), "rho", "tau"),
      column = factor(column, levels = col_order, labels = col_labels),
      is_fdr_col = column %in% c("rho p (FDR)", "tau p (FDR)"),
      is_sig = case_when(
        column == "rho p (FDR)" ~ rho_fdr_sig,
        column == "tau p (FDR)" ~ tau_fdr_sig,
        TRUE ~ FALSE
      )
    )

  ggplot(cell_long, aes(x = column, y = row_label)) +
    geom_tile(aes(fill = is_sig), color = "grey85", linewidth = 0.4) +
    geom_text(aes(label = value, fontface = ifelse(is_sig, "bold", "plain")), size = 5) +
    scale_fill_manual(values = c(`TRUE` = "#d9f2d0", `FALSE` = "white"), guide = "none") +
    scale_x_discrete(position = "top") +
    facet_grid(rows = vars(group), cols = vars(stat_family),
               scales = "free", space = "free", switch = "y") +
    # title is wrapped: with facet_grid + a wide row-label column, ggplot
    # allots the title only the panel's own width (not the full device
    # width including that column) -- an unwrapped long title clips off
    # the right edge of the PNG rather than spanning over the row labels.
    labs(x = NULL, y = NULL,
         title = str_wrap(title_txt, width = 40),
         subtitle = sprintf("Highlighted cells: FDR-corrected p < %.2f   |   %s",
                            fdr_alpha, n_note)) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.8, "cm"),
      panel.spacing.y = unit(0.3, "cm"),
      axis.text.y = element_text(hjust = 1, size = 13),
      axis.text.x = element_text(face = "bold", size = 13),
      strip.placement = "outside",
      strip.text.x = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 14),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
}

table_fig_all <- make_table_figure(
  results_all, "Rank-Order Field vs. Spectral Diversity — All Observations")
table_fig_filtered <- make_table_figure(
  results_filtered, "Rank-Order Field vs. Spectral Diversity — Canopy-Filtered Observations")

ggsave(file.path(fig_dir, "rank_order_table_all_observations.png"), table_fig_all,
       width = 13, height = 0.5 * nrow(results_all) + 2, dpi = 150, bg = "white")
# filtered-pass row labels carry a " (canopy)" suffix (longer than the
# all-observations labels), which leaves less panel width at a fixed
# device width and crowds the rho/tau column headers together -- a bit
# more width keeps them legible.
ggsave(file.path(fig_dir, "rank_order_table_filtered_observations.png"), table_fig_filtered,
       width = 15, height = 0.5 * nrow(results_filtered) + 2, dpi = 150, bg = "white")

cat("  - ", fig_dir, "/rank_order_table_all_observations.png (pass 1, FDR-significant cells highlighted)\n",
    "  - ", fig_dir, "/rank_order_table_filtered_observations.png (pass 2, FDR-significant cells highlighted)\n", sep = "")

# ============================================================================
# 10. Summary figure -- how did canopy-filtering change rank-order agreement,
#     per spectral metric within each group? Visualizes rho_delta from
#     improvement_summary (section 7) directly -- no delta is recomputed
#     here, this is purely a plot of the existing column. Diversity has two
#     field-metric variants (Shannon mean and gamma) per spectral metric, so
#     bars are labeled with the spectral metric name plus that variant where
#     it applies, to avoid two different deltas sharing one label.
# ============================================================================
delta_plot_df <- improvement_summary %>%
  mutate(
    field_variant = case_when(
      pair_key == "floristic_shannon_mean"  ~ " (mean)",
      pair_key == "floristic_shannon_gamma" ~ " (gamma)",
      TRUE ~ ""
    ),
    bar_label = paste0(spectral_var, field_variant),
    bar_label = fct_reorder(bar_label, rho_delta)
  )

delta_summary_plot <- ggplot(delta_plot_df, aes(x = bar_label, y = rho_delta, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  coord_flip() +
  # subtitle wrapped -- at this plot's width, the unwrapped sentence runs
  # off the right edge of the PNG rather than staying inside the canvas
  # (same fix as make_table_figure()'s title above).
  labs(
    x = NULL, y = "rho delta (filtered pass − all-observations pass)",
    fill = "Group",
    title = "Effect of Canopy-Filtering on Rank-Order Agreement",
    subtitle = str_wrap("Positive = canopy-filtering IMPROVED agreement with the spectral metric; negative = it HURT agreement.", width = 70)
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9.5, color = "grey30"))

ggsave(file.path(fig_dir, "canopy_filtering_rho_delta_summary.png"), delta_summary_plot,
       width = 8, height = max(4, 0.4 * nrow(delta_plot_df) + 1.5), dpi = 150, bg = "white")

cat("  - ", fig_dir, "/canopy_filtering_rho_delta_summary.png (rho_delta per spectral metric, from canopy_filtering_improvement_summary.csv)\n", sep = "")
