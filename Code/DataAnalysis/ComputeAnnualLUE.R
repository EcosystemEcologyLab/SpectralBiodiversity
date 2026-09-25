# ============================================================================
# ComputeAnnualLUE.R
#
# Annual Light Use Efficiency (LUE) from REAL downloaded MODIS fPAR
# (MOD15A2H.061 8-day composites, AppEEARS point-sample output) joined
# against REAL annual-resolution AmeriFlux (YY) flux data -- a THIRD,
# independent LUE estimate alongside ComputeLUE.R's half-hourly regression
# (NDVI-derived fPAR) and ComputeLUE_Annual.R's annual ratio (also
# NDVI-derived fPAR). This script is the first to use MODIS fPAR directly as
# the canopy-absorption input, rather than an NDVI-derived proxy. It does
# NOT modify, call, or depend on ExtractMODIS.R, ComputeLUE.R, or
# ComputeLUE_Annual.R -- fully standalone except for reusing two small,
# already-correct helper functions from ComputeLUE_Annual.R (see SOURCING
# below), and reads the raw AppEEARS results CSV directly rather than going
# through ExtractMODIS.R's (untested, no-internet-in-this-sandbox) download
# pipeline.
#
# ==============================================================================
# REAL MODIS DATA -- CORRECTS UNTESTED ASSUMPTIONS IN ExtractMODIS.R
#
# The real downloaded file (Data/MODIS/.../NEON-MODIS-2013-2021-MOD15A2H-061-
# results.csv -- confirmed by direct inspection, not this sandbox: no Data/
# directory exists here, see VALIDATION below) corrects two assumptions
# ExtractMODIS.R had to guess at without internet access:
#   1. NO x0.01 scale factor: real Fpar_500m values are already fractional
#      (0.51, 0.6, 0.75, ...), unlike ExtractMODIS.R's documented-but-
#      unconfirmed fpar_scale_factor <- 0.01. NOT applied here.
#   2. FILL-VALUE TRAP: Fpar_500m contains literal 255.00 sentinels (302
#      occurrences in the real file) -- MODIS's "not processed" fill value,
#      left unconverted by AppEEARS. Any Fpar_500m > 1 is physically
#      impossible and is excluded as invalid, never averaged in as a real
#      255% fPAR observation (fpar_valid_range below).
#   3. QC is PRE-DECODED by AppEEARS into human-readable *_Description
#      columns (MOD15A2H_061_FparLai_QC_MODLAND_Description,
#      ..._CloudState_Description) -- this script uses those directly and
#      does NOT reimplement ExtractMODIS.R's bit-level decode_fpar_qc().
# ==============================================================================
#
# ==============================================================================
# CLOUD-STATE JUDGMENT CALL -- see include_undefined_cloud_state below
#
# The real QC CloudState_Description column has (at least) these values:
#   "Significant clouds NOT present (clear)"        -- unambiguously good
#   "Significant clouds WERE present"                -- always excluded
#   "Mixed cloud present in pixel"                   -- always excluded
#   "Cloud state not defined, assumed clear"         -- AMBIGUOUS: the
#     algorithm couldn't determine cloud state but defaulted to treating the
#     pixel as clear anyway. Whether to trust that default is a judgment
#     call, not a fact -- DEFAULTS TO EXCLUDED here (conservative: an
#     "assumed clear" pixel is a weaker guarantee than a confirmed-clear
#     one, and this script would rather under-use than silently accept
#     lower-confidence composites). Flip include_undefined_cloud_state to
#     TRUE below to include it instead -- one named constant, not a buried
#     conditional.
# ==============================================================================
#
# SOURCING, NOT COPY-PASTING: resolve_par_column() and load_annual_flux_data()
# are imported (not reimplemented) from ComputeLUE_Annual.R via
# import_functions_from() -- the same bootstrap/reuse convention this repo's
# other cross-script imports use (ComputeLUE_Annual.R importing from
# ComputeLUE.R/AnnualSpectralDiversity.R/FieldDiversity.R, CanopyDeltaByEcosystem.R
# importing from CompareSpectralVsFieldDiversity.R). Both are clean,
# standalone top-level functions in ComputeLUE_Annual.R, not entangled with
# its hyperspectral-specific state -- confirmed by reading that script in
# full before importing. Their free variables (annual_data_dir,
# sw_to_ppfd_factor) are resolved via normal lexical scoping against THIS
# script's own Section 0 globals, set to the same values, exactly the
# pattern ComputeLUE_Annual.R itself uses when importing from ComputeLUE.R.
# import_functions_from() itself is necessarily redefined here (it cannot
# import itself without circularity), same as every script in this repo.
#
# ET RESOLUTION: no function for this existed anywhere in the repo to
# import -- Code/NEON_FluxVariability.R does it inline (ET = LE_F_MDS, then
# `ET * 0.0864 / 2.45 * 365.25` to convert W/m2 -> mm/yr, lambda = 2.45
# MJ/kg). resolve_et_column() below is new, but uses that EXACT SAME
# conversion factor (0.0864 / 2.45 * 365.25), not a re-derived one, and
# checks for a literal ET-named column first in case one is actually present
# (investigate-before-assuming, per task instructions) before falling back
# to the confirmed LE_F_MDS conversion.
#
# ANNUAL AMERIFLUX SCHEMA -- confirmed by repo evidence, not this session's
# own file access (this sandbox has no Data/ directory at all -- see
# VALIDATION below). Code/NEON_FluxVariability.R contains real, WORKING
# column-selection code reading this exact directory
# ("./Data/NEON_Ameriflux/AnnualData") today:
#     year = TIMESTAMP, NEE = NEE_VUT_REF, GPP = GPP_NT_VUT_REF,
#     RECO = RECO_NT_VUT_REF, ET = LE_F_MDS (W/m2, converted as above)
# This is direct evidence of real column names in this repo's own real data,
# not an assumption -- but it is still NOT the same as this script itself
# having read a real file in this run. Section 1 below therefore still does
# its OWN runtime investigation (prints names()/head() of the first real
# annual file found) and stop()s with the actual column list if the expected
# names aren't there, rather than trusting the cited evidence blindly.
#
# VALIDATION: this sandbox has neither Data/MODIS/.../NEON-MODIS-...csv nor
# Data/NEON_Ameriflux/AnnualData/*.csv (confirmed via `find /` -- consistent
# with every prior session touching this repo's Data/ directory). Validated
# here against synthetic fixtures built to mirror BOTH the real MODIS
# structure described in the task (including literal 255.00 fill rows and
# all four QC description strings) and the real AmeriFlux annual schema
# confirmed above (including a deliberately wrong/missing-column fixture to
# confirm the Section 1 investigation guard actually stop()s rather than
# silently proceeding). See SESSION_LOG.md for what was checked.
# ============================================================================

library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(readr)

# ============================================================================
# 0. Config
# ============================================================================
modis_fpar_csv  <- "./Data/MODIS/4f0a5de4-0053-433c-87ab-33fc5fd3aadc/NEON-MODIS-2013-2021-MOD15A2H-061-results.csv"
annual_data_dir <- "./Data/NEON_Ameriflux/AnnualData"   # CONFIRMED real path -- see header
out_csv         <- "./Data/lue_by_tower_year.csv"        # matches this project's recent Data/-output convention

annual_flux_script <- "./Code/DataAnalysis/ComputeLUE_Annual.R"

sw_to_ppfd_factor <- 2.02   # only used if resolve_par_column() falls back to SW_IN_F/SW_IN --
                              # free variable inside the IMPORTED resolve_par_column(), see SOURCING above

le_to_et_factor <- 0.0864 / 2.45 * 365.25   # W/m2 -> mm/yr, EXACT SAME formula as
                                              # Code/NEON_FluxVariability.R (lambda = 2.45 MJ/kg)

fpar_valid_range <- c(0, 1)   # physically valid fPAR -- excludes the confirmed 255.00 fill
                                # sentinel (302 occurrences in the real file) and any other
                                # out-of-range value; NO rescaling applied (see header)

good_quality_desc    <- "Good quality (main algorithm with or without saturation)"
clear_cloud_desc     <- "Significant clouds NOT present (clear)"
undefined_cloud_desc <- "Cloud state not defined, assumed clear"

# JUDGMENT CALL -- see header. Change to TRUE to also accept
# undefined_cloud_desc composites; FALSE (default) excludes them.
include_undefined_cloud_state <- FALSE

acceptable_cloud_states <- if (include_undefined_cloud_state) {
  c(clear_cloud_desc, undefined_cloud_desc)
} else {
  c(clear_cloud_desc)
}

if (!file.exists(modis_fpar_csv)) {
  stop("Required input not found: ", modis_fpar_csv,
       " -- this is the real downloaded AppEEARS MOD15A2H point-sample CSV;",
       " if it has moved, update modis_fpar_csv above.")
}

cat("Cloud-state filter: '", clear_cloud_desc, "'",
    if (include_undefined_cloud_state) paste0(" + '", undefined_cloud_desc, "'") else "",
    " (include_undefined_cloud_state = ", include_undefined_cloud_state, ")\n", sep = "")

# ============================================================================
# 1. Import (not reimplement) two helpers from ComputeLUE_Annual.R -- see
#    SOURCING above for why this is safe and why import_functions_from()
#    can't import itself.
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

annual_fns <- import_functions_from(
  annual_flux_script,
  c("load_annual_flux_data", "resolve_par_column")
)

# ============================================================================
# 2. Real-file investigation (task step 3) -- print actual column names and
#    sample rows from a real annual AmeriFlux file BEFORE any logic below
#    depends on a specific column name. Reads the FIRST file directly
#    (separate from the full load_annual_flux_data() call below) so this
#    investigation always happens, and always happens first.
# ============================================================================
if (!dir.exists(annual_data_dir)) {
  stop("Required input directory not found: ", annual_data_dir,
       " -- confirmed real path (see header); if it has moved, update",
       " annual_data_dir above.")
}
annual_files <- list.files(annual_data_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(annual_files) == 0) {
  stop("No CSV files found under ", annual_data_dir, " -- nothing to process.")
}

sample_df <- read.csv(annual_files[1], na.strings = "-9999", stringsAsFactors = FALSE)
cat("\n==== Annual AmeriFlux file investigation (", basename(annual_files[1]), ") ====\n", sep = "")
cat("Columns (", ncol(sample_df), "):\n", sep = "")
print(names(sample_df))
cat("\nSample rows:\n")
print(head(sample_df, 3))

# Expected per repo evidence (Code/NEON_FluxVariability.R's own working
# column selection -- see header). Checked against what was ACTUALLY just
# printed above, not assumed: a real schema mismatch stops here with the
# real column list already visible, rather than silently guessing or
# failing deep inside the main loop.
flux_expected_cols <- c("TIMESTAMP", "GPP_NT_VUT_REF")
missing_flux_cols <- setdiff(flux_expected_cols, names(sample_df))
if (length(missing_flux_cols) > 0) {
  stop("Expected column(s) [", paste(missing_flux_cols, collapse = ", "),
       "] not found in the real annual AmeriFlux file (see columns printed",
       " above). Update this script's column names to match rather than",
       " guessing -- TIMESTAMP/GPP_NT_VUT_REF were expected per",
       " Code/NEON_FluxVariability.R's own working code against this same",
       " directory.")
}

reco_present <- "RECO_NT_VUT_REF" %in% names(sample_df)
nee_present  <- "NEE_VUT_REF" %in% names(sample_df)
cat("\nRECO_NT_VUT_REF present:", reco_present, " | NEE_VUT_REF present:", nee_present, "\n")
if (!reco_present) cat("  -- reco_annual will be NA for all rows (column not found).\n")
if (!nee_present)  cat("  -- nee_annual will be NA for all rows (column not found).\n")

# PAR-equivalent: PPFD_IN tried first, SW_IN_F/SW_IN fallback with
# sw_to_ppfd_factor -- see resolve_par_column() (imported) and header.
par_info <- annual_fns$resolve_par_column(sample_df)
if (is.null(par_info)) {
  stop("No PAR-equivalent column found (tried PPFD_IN, SW_IN_F, SW_IN) in the",
       " real annual AmeriFlux file -- inspect the columns printed above and",
       " extend resolve_par_column() in ", annual_flux_script,
       " rather than guessing.")
}
cat("PAR-equivalent column resolved to '", par_info$col, "' (par_source = '",
    par_info$source, "')",
    if (par_info$source == "sw_in_derived") paste0(" -- converting via sw_to_ppfd_factor = ", sw_to_ppfd_factor) else "",
    "\n", sep = "")

# ET: literal ET column tried first (in case one is actually present --
# investigate, don't assume), LE_F_MDS -> ET conversion as the confirmed
# real-data fallback (see header). NOT imported -- no such function existed
# anywhere in the repo to reuse (see header's ET RESOLUTION note).
resolve_et_column <- function(df) {
  et_direct_candidates <- c("ET", "ET_F_MDS")
  hit <- intersect(et_direct_candidates, names(df))
  if (length(hit) > 0) {
    return(list(col = hit[1], source = "et_direct", conversion = 1))
  }
  if ("LE_F_MDS" %in% names(df)) {
    return(list(col = "LE_F_MDS", source = "le_converted_to_et", conversion = le_to_et_factor))
  }
  NULL
}

et_info <- resolve_et_column(sample_df)
if (is.null(et_info)) {
  cat("No ET or LE column found (tried ET, ET_F_MDS, LE_F_MDS) -- et_annual",
      " will be NA for all rows.\n", sep = "")
} else {
  cat("ET-equivalent column resolved to '", et_info$col, "' (et_source = '",
      et_info$source, "')",
      if (et_info$source == "le_converted_to_et") paste0(" -- converting via le_to_et_factor = ", round(le_to_et_factor, 4)) else "",
      "\n", sep = "")
}

# ============================================================================
# 3. Load the FULL annual AmeriFlux dataset (all sites) and extract, per
#    tower-year: GPP, RECO, NEE, ET, PAR. Reuses load_annual_flux_data()
#    (imported) rather than re-reading the directory a second way.
# ============================================================================
annual_flux <- annual_fns$load_annual_flux_data()
annual_flux$year <- suppressWarnings(as.integer(annual_flux$TIMESTAMP))

flux_annual_tbl <- tibble(
  tower_id    = annual_flux$tower_id,
  year        = annual_flux$year,
  gpp_annual  = annual_flux$GPP_NT_VUT_REF,
  reco_annual = if (reco_present) annual_flux$RECO_NT_VUT_REF else NA_real_,
  nee_annual  = if (nee_present)  annual_flux$NEE_VUT_REF else NA_real_,
  par_annual  = annual_flux[[par_info$col]] * par_info$conversion,
  et_annual   = if (!is.null(et_info)) annual_flux[[et_info$col]] * et_info$conversion else NA_real_
)

cat("\nLoaded annual flux data:", nrow(flux_annual_tbl), "tower-years across",
    n_distinct(flux_annual_tbl$tower_id), "towers.\n")

# ============================================================================
# 4. Real MODIS fPAR: load, QC-filter, aggregate to annual per tower-year.
#    QC uses the PRE-DECODED description columns directly (no bit decoding --
#    see header). fpar_valid_range excludes the confirmed 255.00 fill trap.
# ============================================================================
modis_raw <- read_csv(modis_fpar_csv, show_col_types = FALSE)

modis_required_cols <- c("ID", "Date", "MOD15A2H_061_Fpar_500m",
                          "MOD15A2H_061_FparLai_QC_MODLAND_Description",
                          "MOD15A2H_061_FparLai_QC_CloudState_Description")
missing_modis_cols <- setdiff(modis_required_cols, names(modis_raw))
if (length(missing_modis_cols) > 0) {
  stop("Expected column(s) [", paste(missing_modis_cols, collapse = ", "),
       "] not found in ", modis_fpar_csv, " -- actual columns: ",
       paste(names(modis_raw), collapse = ", "),
       ". Update this script's column names to match rather than guessing.")
}

modis_raw$date <- as.Date(modis_raw$Date)
n_bad_dates <- sum(is.na(modis_raw$date) & !is.na(modis_raw$Date) & modis_raw$Date != "")
if (n_bad_dates > 0) {
  warning(n_bad_dates, " row(s) in ", modis_fpar_csv, " have a Date value that",
          " failed to parse -- inspect the real Date format before trusting",
          " year-level aggregation.")
}
modis_raw$year <- as.integer(format(modis_raw$date, "%Y"))

modis_qc <- modis_raw %>%
  mutate(
    fpar_val = as.numeric(MOD15A2H_061_Fpar_500m),
    qc_passed = MOD15A2H_061_FparLai_QC_MODLAND_Description == good_quality_desc &
      MOD15A2H_061_FparLai_QC_CloudState_Description %in% acceptable_cloud_states &
      !is.na(fpar_val) & fpar_val >= fpar_valid_range[1] & fpar_val <= fpar_valid_range[2]
  )

n_fill_excluded <- sum(!is.na(modis_qc$fpar_val) & modis_qc$fpar_val > fpar_valid_range[2], na.rm = TRUE)
cat("\nMODIS fPAR: ", nrow(modis_qc), " total 8-day composite rows, ",
    sum(modis_qc$qc_passed), " passed QC, ", n_fill_excluded,
    " excluded for being outside the valid fPAR range (", fpar_valid_range[1], "-",
    fpar_valid_range[2], ", catches the confirmed 255.00 fill sentinel).\n", sep = "")

fpar_annual_tbl <- modis_qc %>%
  group_by(tower_id = ID, year) %>%
  summarise(
    n_fpar_obs_total = n(),
    n_fpar_obs_used  = sum(qc_passed, na.rm = TRUE),
    fpar_annual = if (sum(qc_passed, na.rm = TRUE) > 0) mean(fpar_val[qc_passed], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

cat("Aggregated to", nrow(fpar_annual_tbl), "tower-years with at least one MODIS composite (",
    sum(!is.na(fpar_annual_tbl$fpar_annual)), " with >= 1 QC-passed composite).\n")

# ============================================================================
# 5. Join fPAR + flux per tower-year, compute LUE = GPP / (PAR * fPAR).
#    full_join keeps every tower-year present on EITHER side, so a
#    tower-year with MODIS data but no flux data (or vice versa) is still
#    reported with an explicit status rather than silently dropped.
# ============================================================================
combined <- full_join(fpar_annual_tbl, flux_annual_tbl, by = c("tower_id", "year"))

# status precedence (documented, not incidental): fpar missing is checked
# before flux missing when BOTH are missing -- an arbitrary but fixed
# tie-break, not a meaningful priority judgment. RECO/ET/NEE are reported
# for reference only and never affect status -- they aren't LUE inputs.
combined <- combined %>%
  mutate(
    apar_annual = par_annual * fpar_annual,
    status = case_when(
      is.na(fpar_annual) ~ "missing fpar data",
      is.na(gpp_annual) | is.na(par_annual) ~ "missing flux data",
      apar_annual == 0 ~ "divide by zero",
      TRUE ~ "ok"
    ),
    lue_annual = if_else(status == "ok", gpp_annual / apar_annual, NA_real_)
  )

result <- combined %>%
  transmute(
    tower_id, year,
    fpar_annual, n_fpar_obs_used, n_fpar_obs_total,
    gpp_annual, reco_annual, et_annual, nee_annual, par_annual,
    lue_annual, status
  ) %>%
  arrange(tower_id, year)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write_csv(result, out_csv)

cat("\n==== Status summary (", nrow(result), " tower-years) ====\n", sep = "")
print(result %>% count(status, sort = TRUE))

cat("\nSaved:", out_csv, "\n")
