# ============================================================================
# CanopyDeltaByEcosystem.R
#
# Does canopy-filtering's effect on field-measured floristic diversity
# concentrate in certain ecosystem types (working hypothesis, from an
# earlier full-dataset scan: dense-canopy forest sites like GRSM, UNDE,
# HARV)? For each tower, computes the delta between the unfiltered and
# canopy-filtered floristic richness/Shannon-diversity metrics, attaches
# each tower's IGBP ecosystem type from NEONsites.csv, and plots both
# deltas as one bar per tower colored by ecosystem.
#
# DIVERSITY METRIC CHOICE: uses floristic_shannon_MEAN (mean across plots),
# not floristic_shannon_gamma (pooled-plot diversity), as this figure's
# default -- flagged here as a choice, not a settled question. gamma could
# be swapped in trivially (see canopy_metric_cols below) if wanted later.
# ============================================================================

library(tidyverse)

# ============================================================================
# 0. Paths -- fail loudly if inputs are missing, same convention as
#    FieldDiversity.R / CompareSpectralVsFieldDiversity.R. Run from the
#    project root (Rscript Code/CanopyDeltaByEcosystem.R).
# ============================================================================
compare_script <- "./Code/CompareSpectralVsFieldDiversity.R"
neonsites_path <- "./Data/NEONsites.csv"
out_dir        <- "."

required_inputs <- c(compare_script, neonsites_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("CanopyDeltaByEcosystem.R: missing required input(s): ",
       paste(missing_inputs, collapse = ", "),
       ". This script needs CompareSpectralVsFieldDiversity.R (its ",
       "field-resolution logic is imported below, not reimplemented) and ",
       "NEONsites.csv (for the Veg Type / IGBP ecosystem lookup).")
}

# ============================================================================
# 1. Reuse (not reimplement) the field-resolution logic from
#    CompareSpectralVsFieldDiversity.R -- same import_functions_from()
#    convention already established in ComputeLUE.R / ComputeLUE_Annual.R /
#    CompareSSR_AdaptiveFCM_vs_KMeans.R for pulling logic out of a
#    monolithic script that is not a function library (a plain source()
#    would re-run its own downstream rank-order tests, plots, and ggsave()
#    calls as a side effect -- this only evals the specific top-level
#    bindings named below). field_csv -> field -> field_resolved are
#    imported IN DEPENDENCY ORDER so each later eval can see the earlier
#    binding in the same env.
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

field_fns <- import_functions_from(
  compare_script,
  c("field_csv", "field", "field_resolved")
)
field_resolved <- field_fns$field_resolved

cat("Imported field_resolved from", compare_script, "--", nrow(field_resolved),
    "tower-years (peak_flight preferred over per_bout, identical rule to the",
    "spectral comparison).\n")

# ============================================================================
# 2. Per-tower averaging -- the SAME operation CompareSpectralVsFieldDiversity.R
#    applies (group_by(tower_id, neon_site) %>% summarise(across(..., mean,
#    na.rm = TRUE))), run directly on field_resolved instead of on its
#    joined_by_year. That script's version runs on joined_by_year because it
#    needs the spectral join for its own comparisons; this script never
#    touches spectral_diversity_by_year.csv at all, so there is nothing to
#    join before averaging -- applying the identical aggregation call
#    straight to field_resolved is the direct adaptation, not a separate
#    reimplementation of the averaging logic itself.
# ============================================================================
canopy_metric_cols <- c("floristic_richness", "floristic_richness_canopy",
                         "floristic_shannon_mean", "floristic_shannon_mean_canopy")
# swap "floristic_shannon_mean"/"_canopy" for "floristic_shannon_gamma"/"_canopy"
# above (and in the plot call in section 6) to use pooled-plot diversity instead.

per_tower <- field_resolved %>%
  group_by(tower_id, neon_site) %>%
  summarise(n_years_averaged = n(),
            across(all_of(canopy_metric_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

cat("Averaged to", nrow(per_tower), "unique towers.\n")

# ============================================================================
# 3. Investigate NEONsites.csv's Veg Type column BEFORE mapping anything.
#    Print the actual column name, unique values, counts, and sample rows --
#    do NOT assume numeric IGBP codes (1-17), 3-letter abbreviations
#    (ENF/DBF/...), or full names until confirmed here.
# ============================================================================
sites_raw <- read.csv(neonsites_path, fileEncoding = "UTF-8-BOM")

veg_col <- names(sites_raw)[str_detect(names(sites_raw), regex("veg.?type", ignore_case = TRUE))]
if (length(veg_col) != 1) {
  stop("Expected exactly one Veg Type column in ", neonsites_path,
       ", found ", length(veg_col), ": ", paste(veg_col, collapse = ", "),
       ". Not guessing which one is intended.")
}

cat("\n==== NEONsites.csv Veg Type investigation ====\n")
cat("Column name (as read into R):", veg_col, "\n")
cat("nrow:", nrow(sites_raw), "\n")
cat("\nUnique values (", length(unique(sites_raw[[veg_col]])), "):\n", sep = "")
print(sort(unique(sites_raw[[veg_col]])))
cat("\nValue counts:\n")
print(table(sites_raw[[veg_col]], useNA = "always"))
cat("\nSample rows:\n")
print(head(sites_raw, 5))

# ---- ASSUMPTION TO CONFIRM AGAINST THE PRINTOUT ABOVE: the values are
# standard 3-letter IGBP land-cover abbreviations (matching AmeriFlux's own
# site-page vegetation-abbreviation convention -- consistent with this
# file's "Hub" column reading "AmeriFlux" where that column exists), NOT
# numeric 1-17 codes and NOT full names. If the printed values above don't
# look like this, this mapping is wrong and must be revisited before
# trusting any ecosystem grouping in the figures below -- the stop() right
# after the legend catches an unrecognized value rather than silently
# leaving it unmapped. ----
igbp_legend <- tribble(
  ~code, ~ecosystem,
  "ENF", "Evergreen Needleleaf Forest",
  "EBF", "Evergreen Broadleaf Forest",
  "DNF", "Deciduous Needleleaf Forest",
  "DBF", "Deciduous Broadleaf Forest",
  "MF",  "Mixed Forest",
  "CSH", "Closed Shrublands",
  "OSH", "Open Shrublands",
  "WSA", "Woody Savannas",
  "SAV", "Savannas",
  "GRA", "Grasslands",
  "WET", "Permanent Wetlands",
  "CRO", "Croplands",
  "URB", "Urban and Built-up",
  "CVM", "Cropland/Natural Vegetation Mosaic",
  "SNO", "Snow and Ice",
  "BSV", "Barren or Sparsely Vegetated",
  "WAT", "Water Bodies"
)

observed_codes <- unique(sites_raw[[veg_col]])
unmapped_codes <- setdiff(observed_codes, igbp_legend$code)
if (length(unmapped_codes) > 0) {
  stop("NEONsites.csv's '", veg_col, "' column contains value(s) not in the ",
       "standard IGBP 17-class abbreviation legend assumed above: ",
       paste(unmapped_codes, collapse = ", "),
       ". Confirm the real encoding (numeric? full names? a different ",
       "vocabulary entirely?) against the printout above before extending ",
       "this mapping -- not guessing.")
}

sites_lookup <- sites_raw %>%
  transmute(tower_id  = Site.ID,
            site_name = Site.Name,
            veg_code  = .data[[veg_col]],
            ecosystem = igbp_legend$ecosystem[match(.data[[veg_col]], igbp_legend$code)])

cat("\nMapped", nrow(sites_lookup), "sites to", n_distinct(sites_lookup$ecosystem),
    "IGBP ecosystem labels (of", nrow(igbp_legend), "in the full legend).\n")

# ============================================================================
# 4. Join per-tower metrics to ecosystem type on tower_id -- tower_id IS
#    Site.ID (AmeriFlux-style code, e.g. "US-xAB"), matching the SAME
#    site_xwalk convention FieldDiversity.R already establishes
#    (site_xwalk <- towers_df %>% transmute(tower_id = Site.ID, neon_site)),
#    not the 4-letter NEON code embedded in Site.Name.
# ============================================================================
per_tower_eco <- per_tower %>%
  left_join(sites_lookup, by = "tower_id")

unmatched <- per_tower_eco %>% filter(is.na(ecosystem))
if (nrow(unmatched) > 0) {
  cat("\n!!", nrow(unmatched), "tower(s) have field data but no ecosystem match in",
      neonsites_path, "-- listed below, plotted under an explicit \"Unknown\"",
      "category rather than silently dropped:\n")
  print(as.data.frame(unmatched %>% select(tower_id, neon_site)))
  per_tower_eco <- per_tower_eco %>%
    mutate(ecosystem = if_else(is.na(ecosystem), "Unknown", ecosystem))
} else {
  cat("\nAll", nrow(per_tower_eco), "towers matched an ecosystem type in", neonsites_path, ".\n")
}

# ============================================================================
# 5. Deltas. Positive richness_delta/shannon_delta = canopy-filtering
#    EXCLUDED species/lowered diversity relative to the unfiltered metric.
#    pct_filtered guards against a floristic_richness of 0 (division by
#    zero) by reporting NA rather than Inf/NaN for that tower.
# ============================================================================
deltas <- per_tower_eco %>%
  mutate(
    richness_delta = floristic_richness - floristic_richness_canopy,
    shannon_delta  = floristic_shannon_mean - floristic_shannon_mean_canopy,
    pct_filtered   = if_else(floristic_richness > 0,
                              100 * (floristic_richness - floristic_richness_canopy) / floristic_richness,
                              NA_real_),
    tower_label    = coalesce(neon_site, tower_id)
  )

write_csv(deltas, file.path(out_dir, "canopy_delta_by_tower.csv"))

cat("\nSummary of richness_delta by ecosystem (mean, n towers):\n")
print(deltas %>% group_by(ecosystem) %>%
        summarise(mean_richness_delta = mean(richness_delta, na.rm = TRUE),
                  mean_pct_filtered   = mean(pct_filtered, na.rm = TRUE),
                  n_towers = n(), .groups = "drop") %>%
        arrange(desc(mean_richness_delta)))

# ============================================================================
# 6. Figures -- one bar per tower, ALL towers included (zero-delta bars are
#    part of what shows whether the effect concentrates by ecosystem type),
#    sorted by delta descending, fill = IGBP ecosystem, tower/site labels
#    rotated for the 44+-bar width.
# ============================================================================
make_delta_plot <- function(df, delta_col, y_lab, title_txt, annotate_pct = FALSE) {
  d <- df %>%
    mutate(tower_label = fct_reorder(tower_label, .data[[delta_col]], .desc = TRUE))

  p <- ggplot(d, aes(x = tower_label, y = .data[[delta_col]], fill = ecosystem)) +
    geom_col() +
    labs(x = NULL, y = y_lab, fill = "Ecosystem (IGBP)", title = title_txt) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          plot.title = element_text(face = "bold"))

  if (annotate_pct) {
    d_lab <- d %>%
      mutate(pct_lab   = if_else(!is.na(pct_filtered) & abs(pct_filtered) > 0.05,
                                  sprintf("%.1f%%", pct_filtered), ""),
             lab_vjust = if_else(.data[[delta_col]] >= 0, -0.3, 1.1))
    p <- p + geom_text(data = d_lab, aes(label = pct_lab, vjust = lab_vjust),
                        size = 2.2, color = "grey20")
  }
  p
}

richness_plot <- make_delta_plot(
  deltas, "richness_delta",
  "Richness delta (unfiltered − canopy-filtered)",
  "Canopy-Filtering Effect on Field Richness, by Tower and Ecosystem",
  annotate_pct = TRUE
)

diversity_plot <- make_delta_plot(
  deltas, "shannon_delta",
  "Shannon diversity delta, mean (unfiltered − canopy-filtered)",
  "Canopy-Filtering Effect on Field Shannon Diversity, by Tower and Ecosystem",
  annotate_pct = FALSE
)

fig_width <- max(10, 0.22 * nrow(deltas))
ggsave(file.path(out_dir, "canopy_delta_richness_by_ecosystem.png"), richness_plot,
       width = fig_width, height = 6, dpi = 150, limitsize = FALSE)
ggsave(file.path(out_dir, "canopy_delta_diversity_by_ecosystem.png"), diversity_plot,
       width = fig_width, height = 6, dpi = 150, limitsize = FALSE)

cat("\nSaved:\n",
    "  - canopy_delta_by_tower.csv\n",
    "  - canopy_delta_richness_by_ecosystem.png\n",
    "  - canopy_delta_diversity_by_ecosystem.png\n", sep = "")
