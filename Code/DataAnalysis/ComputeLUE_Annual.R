# ComputeLUE_Annual.R
#
# Light Use Efficiency (LUE, epsilon) from ANNUAL-resolution AmeriFlux (YY)
# data, for the 45 NEON/AmeriFlux sites in this project. A SEPARATE, PARALLEL
# script to Code/DataAnalysis/ComputeLUE.R (the half-hourly, regression-slope
# version) -- ComputeLUE.R is NOT modified by this script and stays exactly
# as it was, in case the half-hourly approach is revisited later. Different
# output file (lue_by_year_annual.csv vs. lue_by_year.csv), so the two
# methodologically-distinct results can never overwrite or be confused with
# each other.
#
# ==============================================================================
# METHODOLOGICAL DISTINCTION FROM ComputeLUE.R -- READ BEFORE COMPARING OUTPUTS
#
# ComputeLUE.R computes LUE as the SLOPE of lm(GPP ~ APAR) across many
# midday half-hourly observations within a site-year, specifically to avoid
# ratio-of-means bias (see its header) and to get an R^2 QA diagnostic for
# free. Annual-resolution AmeriFlux data gives exactly ONE GPP value and ONE
# PAR-equivalent value per site-year -- there is no within-period variation
# left to fit a slope to, so the regression approach is categorically
# UNAVAILABLE here, not merely skipped. This script instead computes LUE as
# the SIMPLE RATIO:
#
#     LUE = GPP_annual / (PAR_annual * fPAR)
#
# This is a DELIBERATE, DIFFERENT methodological choice forced by the data
# resolution, not an oversight or a simplification of ComputeLUE.R's method.
# It does NOT have the same protection against ratio-of-means bias that the
# regression-slope approach was specifically built to avoid -- a ratio of two
# annual means/means-of-rates is exactly the "ratio of period means" pattern
# ComputeLUE.R's header describes as biased whenever the numerator and
# denominator co-vary within the averaging period (which annual GPP and PAR
# certainly do, tracking the seasonal cycle together). CONSEQUENCE: this
# script's `lue` values are NOT expected to be numerically identical to what
# a hypothetical regression-based annual approach would give, and the two
# scripts' `lue` columns must be treated as methodologically distinct
# outputs for comparison/sensitivity purposes, never as interchangeable or
# as a single merged "the" LUE value.
#
# A SECOND, DISTINCT caveat (not the ratio-vs-regression one above): the
# fPAR/NDVI value is still a SINGLE snapshot from one hyperspectral flight
# day, exactly as in ComputeLUE.R. ComputeLUE.R handles this by restricting
# its GPP/APAR data to the SAME calendar window the snapshot represents
# (growing_season_source / gs$months), so the flux data being fit is
# temporally consistent with the canopy state the snapshot captures. At
# annual resolution there is no such restriction available -- gpp_annual and
# par_annual are, by construction, whole-year aggregates covering dormant-
# season months right alongside peak-growing-season months, and cannot be
# subset the way half-hourly rows can. Dividing a whole-year GPP/PAR ratio by
# a single-day snapshot therefore conflates seasons in a way the half-hourly
# version specifically avoids -- this script's LUE is a coarser, more
# seasonally-blended quantity, not just a resolution-downgraded version of
# ComputeLUE.R's. `growing_season_source` is still resolved and reported
# (real flight month vs. fallback, for provenance/consistency with the other
# two scripts and because it identifies which month the fPAR snapshot
# represents), but it is NOT used as an active data-window filter here --
# there is nothing left at annual resolution to filter.
# ==============================================================================
#
# ----------------------------------------------------------------------------
# BEFORE-WRITING INVESTIGATION (per task instructions -- see SESSION_LOG.md):
#
# ANNUAL FILE PATH/STRUCTURE -- CONFIRMED, not guessed (unlike ComputeLUE.R's
# HH file path, which remains an unconfirmed best-guess): Code/NEON_FluxVariability.R
# already reads this exact directory successfully in this repo:
#   annual_data_dir <- "./Data/NEON_Ameriflux/AnnualData"   -- one YY CSV per site
# That script's own column selection (`select(site, year = TIMESTAMP, NEE =
# NEE_VUT_REF, GPP = GPP_NT_VUT_REF, RECO = RECO_NT_VUT_REF, ET = LE_F_MDS)`)
# is direct repo evidence -- not an assumption -- that the real YY file has
# these columns UNCHANGED from their half-hourly names: TIMESTAMP (the bare
# 4-digit year, e.g. 2017 -- confirmed by `year = TIMESTAMP` with no date
# parsing), NEE_VUT_REF, GPP_NT_VUT_REF, RECO_NT_VUT_REF, LE_F_MDS. So
# GPP_NT_VUT_REF is CONFIRMED present at annual resolution, not assumed.
#
# PAR-EQUIVALENT COLUMN -- NOT confirmed by that evidence (NEON_FluxVariability.R
# never selects a PAR/radiation column, so its presence or absence in the
# real file is genuinely unknown here). resolve_par_column() (Section 2)
# tries "PPFD_IN" first -- plausible given every OTHER confirmed column kept
# its half-hourly name unchanged, but still unverified for this specific
# column -- and falls back to "SW_IN_F"/"SW_IN" (broadband shortwave, W/m2)
# with an explicit, literature-default W/m2->umol-photons-m-2-s-1 conversion
# factor (2.02, McCree/Britton & Dodd; NOT site-calibrated, same style of
# caveat as the NDVI_min/NDVI_max fPAR defaults) if PPFD_IN isn't present.
# Which path was taken is recorded in the `par_source` output column
# (ppfd_in_annual vs. sw_in_derived) -- an addition beyond the literally
# requested output columns, justified by this project's established
# never-silently-blend-methods convention (growing_season_source,
# real_metadata/fallback_default in the other two LUE scripts).
#
# UNITS CONSISTENCY -- flagged explicitly per task instructions, not just
# assumed correct. FLUXNET2015/AmeriFlux YY aggregation, per the standard
# product convention AND per repo evidence (NEON_FluxVariability.R's own
# comment: "convert LE (W/m2) to ET (mm/yr) if you need absolute ET
# values... LE_F_MDS ... not mm" -- i.e. LE_F_MDS in the YY file is STILL a
# W/m2 flux-density MEAN, not a converted/summed total), reports annual
# ARITHMETIC MEANS of the underlying flux-rate variables, not summed
# totals. On that basis, GPP_NT_VUT_REF and a PPFD_IN-named annual column
# (if present) are both expected to remain in their original per-second
# flux-rate units (umol CO2 m-2 s-1 / umol photons m-2 s-1) -- i.e. the SAME
# units convention as ComputeLUE.R's half-hourly values, just time-averaged
# over the year -- making GPP_annual / (PAR_annual * fPAR) unit-consistent
# in the same mol-CO2/mol-photon LUE convention ComputeLUE.R's regression
# slope uses. This reasoning is evidence-informed, not verified against a
# real downloaded value. check_gpp_plausible_range() (Section 3) is a
# concrete runtime guard against the specific units failure mode the task
# called out -- an annual SUM (e.g. gC/m2/yr, typically ~200-3000) landing
# where an annual MEAN flux rate (umol CO2 m-2 s-1, typically roughly -5 to
# 30) is expected -- flagged in `status`, not silently computed through.
#
# -9999 fill-value handling and Site.ID<->neon_site matching are reused
# verbatim (repo-wide convention, same as ComputeLUE.R/ExtractMODIS.R/
# NEON_FluxVariability.R).
#
# SOURCING, NOT COPY-PASTING (per task instructions -- read ComputeLUE.R in
# full before writing this): import_functions_from() (the bootstrap utility
# itself) is necessarily redefined here rather than imported -- it cannot
# import itself without circularity, the same way every script in this repo
# independently states its own library() calls. Everything it is USED to
# pull in is reused, not reimplemented:
#   - from AnnualSpectralDiversity.R (same set ComputeLUE.R imports):
#     bad_band_ranges, read_neon_h5_tile, mask_bad_bands, nearest_band,
#     get_tower_reflectance, compute_ndvi_raster, extract_year_from_path,
#     discover_site_years.
#   - from FieldDiversity.R: get_flight_acquisition_date.
#   - from ComputeLUE.R ITSELF (the task's suggestion, followed after
#     confirming these are cleanly standalone top-level functions there,
#     not entangled with ComputeLUE.R's half-hourly-specific state):
#     compute_mean_ndvi, compute_fpar, resolve_growing_season. These three
#     reference the free variables `spec_fns`, `ndvi_thresh`, `ndvi_min`,
#     `ndvi_max`, `field_fns`, and `fallback_growing_season_months` --
#     resolved via normal lexical scoping against THIS script's own globals
#     (Section 0/1 below), which are set to the same values/imports
#     ComputeLUE.R uses, for exactly the reasons the task asked to reuse
#     them (same literature-default fPAR calibration limitation, same
#     real-flight-month-vs-fallback growing-season logic).
#
# Output: one row per tower-YEAR, written incrementally after each site-year
# (crash-safe), to out_csv (Section 0).
#
# VALIDATION: synthetic annual-format fixtures only (no real Data/ access in
# this sandbox) -- see SESSION_LOG.md for what was checked (including a
# deliberately wrong-units GPP fixture to confirm check_gpp_plausible_range()
# actually catches it) and what still needs a real-data run on the server.

library(terra)
library(rhdf5)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# ============================================================================
# 0. Config
# ============================================================================
hyperspec_dir    <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
annual_data_dir  <- "./Data/NEON_Ameriflux/AnnualData"   # CONFIRMED real path -- see header
out_csv          <- "D:/projects/moore/SpectralBiodiversity/Data/lue_by_year_annual.csv"

spectral_diversity_script <- "./Code/AnnualSpectralDiversity.R"
field_diversity_script    <- "./Code/DataAnalysis/FieldDiversity.R"
lue_hh_script              <- "./Code/DataAnalysis/ComputeLUE.R"

buffer_m            <- 500    # matches ComputeLUE.R / AnnualSpectralDiversity.R
ndvi_thresh          <- 0.4   # veg-mask threshold, matches ComputeLUE.R
ndvi_min             <- 0.05  # fPAR literature default (NOT site-calibrated -- see ComputeLUE.R header)
ndvi_max             <- 0.95  # fPAR literature default (NOT site-calibrated -- see ComputeLUE.R header)
fallback_growing_season_months <- 6:9   # June-Sept, same convention as ComputeLUE.R

sw_to_ppfd_factor <- 2.02     # umol photons per J of broadband shortwave -- literature default
                               # (McCree 1972 / Britton & Dodd 1976 convention), NOT site-calibrated;
                               # only used if PPFD_IN is absent and SW_IN_F/SW_IN is used instead --
                               # see resolve_par_column() and the header's PAR-EQUIVALENT COLUMN note.

gpp_plausible_range <- c(-5, 30)   # umol CO2 m-2 s-1, plausible ANNUAL MEAN flux-rate range for
                                     # terrestrial GPP -- an annual SUM (e.g. gC/m2/yr, typically
                                     # ~200-3000) would land far outside this; see header's UNITS
                                     # CONSISTENCY note. A value outside this range is flagged in
                                     # `status`, not silently divided through.

towers_df <- read.csv("./Data/NEONsites.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites in NEONsites.csv.\n")

# ============================================================================
# 1. Import (not copy) shared helpers -- see header's SOURCING note
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

spec_fns <- import_functions_from(
  spectral_diversity_script,
  c("bad_band_ranges", "read_neon_h5_tile", "mask_bad_bands", "nearest_band",
    "get_tower_reflectance", "compute_ndvi_raster", "extract_year_from_path",
    "discover_site_years")
)

field_fns <- import_functions_from(
  field_diversity_script,
  c("get_flight_acquisition_date")
)

# compute_mean_ndvi/compute_fpar/resolve_growing_season, imported from
# ComputeLUE.R itself rather than re-derived -- see header's SOURCING note
# for why this is safe (their free variables resolve against this script's
# own spec_fns/field_fns/ndvi_thresh/ndvi_min/ndvi_max/fallback_growing_season_months,
# defined above/below with matching values) and why import_functions_from()
# is not itself importable here (it's the bootstrap; no circularity avoidance
# possible, same as every script's own library() calls).
lue_fns <- import_functions_from(
  lue_hh_script,
  c("compute_mean_ndvi", "compute_fpar", "resolve_growing_season")
)

# ============================================================================
# 2. Annual AmeriFlux (YY) file loading + PAR-equivalent column resolution
# ============================================================================
load_annual_flux_data <- function() {
  if (!dir.exists(annual_data_dir)) {
    stop("Required input directory not found: ", annual_data_dir,
         " -- this is the CONFIRMED real annual AmeriFlux directory (see",
         " NEON_FluxVariability.R / SESSION_LOG.md); if it has moved, update",
         " annual_data_dir above.")
  }
  files <- list.files(annual_data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No CSV files found under ", annual_data_dir, " -- nothing to process.")
  }
  raw <- files %>%
    setNames(files) %>%
    map_dfr(read.csv, na.strings = "-9999", stringsAsFactors = FALSE, .id = "source_file")
  raw$tower_id <- str_extract(basename(raw$source_file), "US-[A-Za-z0-9]+")
  as_tibble(raw)
}

# PPFD_IN tried first (plausible given every other confirmed YY column kept
# its HH name unchanged), SW_IN_F/SW_IN as a documented fallback requiring
# sw_to_ppfd_factor -- see header's PAR-EQUIVALENT COLUMN note. Resolved ONCE
# against the file schema (not per-row): all YY files share one column set.
resolve_par_column <- function(df) {
  if ("PPFD_IN" %in% names(df)) {
    return(list(col = "PPFD_IN", source = "ppfd_in_annual", conversion = 1))
  }
  sw_hit <- intersect(c("SW_IN_F", "SW_IN"), names(df))
  if (length(sw_hit) > 0) {
    return(list(col = sw_hit[1], source = "sw_in_derived", conversion = sw_to_ppfd_factor))
  }
  NULL
}

check_gpp_plausible_range <- function(gpp_val) {
  if (is.na(gpp_val)) return(TRUE)
  gpp_val >= gpp_plausible_range[1] && gpp_val <= gpp_plausible_range[2]
}

# ============================================================================
# 3. LUE ratio calculation (Section header note: see METHODOLOGICAL
#    DISTINCTION FROM ComputeLUE.R at the top of this file)
# ============================================================================
compute_annual_lue <- function(gpp_annual, par_annual, fpar_val) {
  apar_annual <- par_annual * fpar_val
  lue <- if (is.na(apar_annual) || apar_annual == 0) NA_real_ else gpp_annual / apar_annual
  list(apar_annual = apar_annual, lue = lue)
}

# ============================================================================
# 4. Main loop: one row per tower-YEAR
# ============================================================================
make_row <- function(tower_id, neon_site, year, gpp_annual = NA_real_, par_annual = NA_real_,
                      par_source = NA_character_, ndvi = NA_real_, fpar = NA_real_,
                      apar_annual = NA_real_, lue = NA_real_,
                      growing_season_source = NA_character_, status) {
  tibble(tower_id = tower_id, neon_site = neon_site, year = year,
         gpp_annual = gpp_annual, par_annual = par_annual, par_source = par_source,
         ndvi = ndvi, fpar = fpar, apar_annual = apar_annual, lue = lue,
         growing_season_source = growing_season_source, status = status)
}

run_lue_annual_pipeline <- function() {
  annual_flux <- load_annual_flux_data()
  annual_flux$year <- suppressWarnings(as.integer(annual_flux$TIMESTAMP))

  if (!("GPP_NT_VUT_REF" %in% names(annual_flux))) {
    stop("GPP_NT_VUT_REF not found in the annual AmeriFlux files -- this column was",
         " expected to be present unchanged from the half-hourly naming (confirmed",
         " elsewhere by NEON_FluxVariability.R's own successful use of it); inspect",
         " the real file header and update this script rather than guessing.")
  }
  par_info <- resolve_par_column(annual_flux)
  if (is.null(par_info)) {
    stop("No PAR-equivalent column found in the annual AmeriFlux files (tried PPFD_IN,",
         " SW_IN_F, SW_IN) -- inspect the real file header (names(annual_flux)) and",
         " add the real column name to resolve_par_column() rather than guessing.")
  }
  cat("Annual PAR-equivalent column resolved to '", par_info$col, "' (par_source = '",
      par_info$source, "')", if (par_info$source == "sw_in_derived")
        paste0(" -- converting via sw_to_ppfd_factor = ", sw_to_ppfd_factor) else "", "\n", sep = "")

  site_years_files <- setNames(map(towers_df$Site.ID, spec_fns$discover_site_years), towers_df$Site.ID)

  site_year_jobs <- list()
  for (tid in names(site_years_files)) {
    df <- site_years_files[[tid]]
    if (nrow(df) == 0) next
    for (yr in sort(unique(df$year))) {
      site_year_jobs[[length(site_year_jobs) + 1]] <- list(
        tower_id = tid, year = yr, files = df$file[df$year == yr]
      )
    }
  }

  cat("Total site-year jobs to process:", length(site_year_jobs), "\n")
  if (length(site_year_jobs) == 0) {
    stop("No site-years with hyperspectral data found under ", hyperspec_dir, " -- nothing to process.")
  }

  results <- tibble(tower_id = character(), neon_site = character(), year = character(),
                     gpp_annual = double(), par_annual = double(), par_source = character(),
                     ndvi = double(), fpar = double(), apar_annual = double(), lue = double(),
                     growing_season_source = character(), status = character())

  pb <- txtProgressBar(min = 0, max = length(site_year_jobs), style = 3)

  for (j in seq_along(site_year_jobs)) {
    job <- site_year_jobs[[j]]
    tower_id  <- job$tower_id
    yr        <- job$year
    neon_site <- towers_df$neon_site[towers_df$Site.ID == tower_id][1]
    target_year <- suppressWarnings(as.integer(yr))

    cat("\n==== [", j, "/", length(site_year_jobs), "] ", tower_id, " - ", yr, " ====\n", sep = "")

    row <- tryCatch({
      ndvi_val <- lue_fns$compute_mean_ndvi(tower_id, job$files, buffer_m)

      if (is.na(ndvi_val)) {
        make_row(tower_id, neon_site, yr, status = "no vegetation pixels for NDVI")
      } else {
        fpar_val <- lue_fns$compute_fpar(ndvi_val)
        gs <- lue_fns$resolve_growing_season(job$files)

        annual_row <- annual_flux %>% filter(tower_id == !!tower_id, year == target_year)

        if (nrow(annual_row) == 0) {
          make_row(tower_id, neon_site, yr, ndvi = ndvi_val, fpar = fpar_val,
                    growing_season_source = gs$source,
                    status = "no annual AmeriFlux data for this tower-year")
        } else {
          gpp_val <- annual_row$GPP_NT_VUT_REF[1]
          par_raw <- annual_row[[par_info$col]][1]
          par_val <- if (is.na(par_raw)) NA_real_ else par_raw * par_info$conversion

          if (is.na(gpp_val) || is.na(par_val)) {
            make_row(tower_id, neon_site, yr, gpp_annual = gpp_val, par_annual = par_val,
                      par_source = par_info$source, ndvi = ndvi_val, fpar = fpar_val,
                      growing_season_source = gs$source,
                      status = "missing annual GPP and/or PAR-equivalent value (fill/-9999)")
          } else {
            calc <- compute_annual_lue(gpp_val, par_val, fpar_val)
            units_ok <- check_gpp_plausible_range(gpp_val)
            status <- if (is.na(calc$lue)) {
              "zero or missing APAR_annual, cannot compute LUE"
            } else if (!units_ok) {
              paste0("success (WARNING: GPP_NT_VUT_REF = ", round(gpp_val, 2),
                     " is outside the plausible annual-mean flux-rate range [",
                     gpp_plausible_range[1], ", ", gpp_plausible_range[2],
                     "] -- may indicate a units mismatch, e.g. a SUMMED annual total",
                     " rather than a mean flux rate; verify against the real file header)")
            } else {
              "success"
            }
            make_row(tower_id, neon_site, yr, gpp_annual = gpp_val, par_annual = par_val,
                      par_source = par_info$source, ndvi = ndvi_val, fpar = fpar_val,
                      apar_annual = calc$apar_annual, lue = calc$lue,
                      growing_season_source = gs$source, status = status)
          }
        }
      }
    }, error = function(e) {
      cat("  unexpected error:", conditionMessage(e), "\n")
      make_row(tower_id, neon_site, yr, status = paste("error:", conditionMessage(e)))
    })

    results <- bind_rows(results, row)
    # save incrementally so a crash/interruption partway through doesn't lose
    # completed site-years (matches ComputeLUE.R / AnnualSpectralDiversity.R)
    write.csv(results, out_csv, row.names = FALSE)

    cat("  status:", row$status, "\n")
    setTxtProgressBar(pb, j)
  }

  close(pb)
  cat("\nDone. Summary written to:", out_csv, "\n")
  results
}

results_annual <- run_lue_annual_pipeline()
