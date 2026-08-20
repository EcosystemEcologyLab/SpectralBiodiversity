# ComputeLUE.R
#
# Light Use Efficiency (LUE, epsilon) from AmeriFlux/FLUXNET half-hourly (HH)
# flux tower data, for the 45 NEON/AmeriFlux sites in this project. Standalone
# analysis, fully independent of the spectral-diversity pipeline
# (AnnualSpectralDiversity.R) and the footprint pipeline
# (NEON_ExtractFootprints.R / AnnualFootprintExtent.R) -- different script,
# different output path, no shared state -- but follows their conventions
# (site discovery from NEONsites.csv, crash-safe incremental write.csv, one
# row per tower-YEAR).
#
# ----------------------------------------------------------------------------
# BEFORE-WRITING INVESTIGATION (see SESSION_LOG.md for the full report):
#
# - No script in this repo reads AmeriFlux HH files today, and this sandbox
#   has no access to the real D:/X: data drives, so the real HH filename and
#   directory location could NOT be confirmed against an actual file. File
#   discovery below therefore uses a FLEXIBLE pattern (site ID + "FLUXNET" +
#   "_HH_" all present in the basename) rather than a rigid version-number
#   regex, so it tolerates whatever exact FULLSET/FLUXMET middle token and
#   v#-#/_r# suffix the real files use. `hh_data_dir` (Section 0) is a
#   best-guess sibling of NEON_FluxVariability.R's
#   "./Data/NEON_Ameriflux/AnnualData" -- CONFIRM/ADJUST THIS PATH ON THE
#   SERVER before trusting a run that reports "no AmeriFlux HH data found"
#   for every site.
# - Site ID matching (US-xAB style Site.ID <-> 4-letter neon_site) and the
#   -9999 fill-value convention are both established elsewhere in this repo
#   (NEON_ExtractFootprints.R, AnnualFootprintExtent.R, AnnualSpectralDiversity.R,
#   NEON_FluxVariability.R) and are reused verbatim here.
# - NDVI per site-year is NOT already available anywhere in this repo's
#   output -- spectral_diversity_by_year.csv stores raoq_ndvi (a DIVERSITY
#   metric of NDVI), never a plain mean-NDVI value -- so it's computed fresh
#   here from the same hyperspectral archive AnnualSpectralDiversity.R uses.
# - QC-flag filtering has no established convention anywhere in this repo
#   (nothing currently filters on it). This script applies the standard
#   AmeriFlux/FLUXNET2015 practice of keeping NEE_VUT_REF_QC <= qc_max_allowed
#   (0 = measured, 1 = good-quality gapfill) as a project-specific choice, not
#   an inherited one -- see qc_max_allowed in Section 0. If NEE_VUT_REF_QC is
#   absent from the real file, QC filtering is skipped (not treated as fatal).
# - TIMESTAMP_START/END format (YYYYMMDDHHMM) is assumed per the standard
#   FLUXNET2015/AmeriFlux HH product and could not be confirmed against a
#   real file; read_hh_flux() stop()s with a clear message if the expected
#   columns are absent, rather than silently proceeding on a wrong schema.
#
# COLLECTION-MONTH RECOVERY: read_neon_h5_tile() (AnnualSpectralDiversity.R)
# carries no date field, but Code/DataAnalysis/FieldDiversity.R already solved
# this problem -- get_flight_acquisition_date() recovers the REAL flight date
# (not just year) from the ATCOR_Input_file text log embedded under
# Metadata/Logs/<flightline>/, confirmed against a real uploaded H5 tile in a
# prior session. That function is imported (Section 2) and its date's MONTH
# is used as the site-year's growing-season window when resolvable -- i.e.
# the flux data analyzed is restricted to the calendar month the canopy was
# actually imaged, so the single NDVI/fPAR snapshot is temporally consistent
# with the GPP/APAR half-hours being regressed. When no real date is
# recoverable (H5 read/parse failure, or no hyperspec files for that
# tower-year), a literature default growing-season window (Section 0) is used
# instead. Every output row's `growing_season_source` column says which
# happened -- never silently blended.
#
# SOURCING, NOT COPY-PASTING: read_neon_h5_tile, mask_bad_bands, nearest_band,
# get_tower_reflectance, compute_ndvi_raster, extract_year_from_path, and
# discover_site_years are reused from AnnualSpectralDiversity.R, and
# get_flight_acquisition_date from FieldDiversity.R -- via import_functions_from()
# (Section 2), which parse()s each origin script and eval()s ONLY the named
# top-level `<-` bindings into a dedicated environment. Plain source() of
# either origin script is NOT safe here: both are monolithic scripts whose top
# level immediately executes a full site-year pipeline (H5 discovery, main
# loop, its own write.csv) as a side effect the moment they're loaded -- they
# are not function libraries. import_functions_from() keeps the actual logic
# single-sourced in the origin files (a later edit there is picked up
# automatically here) without triggering either pipeline.
#
# METHODOLOGY:
#   LUE (epsilon) = the SLOPE of lm(GPP_NT_VUT_REF ~ APAR) across filtered
#   midday half-hourly records for a site-year, NOT a ratio of period means --
#   a ratio-of-means is biased whenever APAR varies within the averaging
#   window (which it always does across a multi-week growing-season window),
#   and the regression form gives an R^2 per site-year as a free QA
#   diagnostic: a poor fit flags a site-year where the light-response
#   assumption may be violated (e.g. water-stress decoupling of light and
#   carbon uptake), surfaced in `lue_r2` rather than hidden.
#   Fit INCLUDES an intercept (lm(GPP ~ APAR), not lm(GPP ~ 0 + APAR)) --
#   forcing the fit through the origin is a stronger assumption than the
#   light-response literature generally requires; flagged here as a modeling
#   choice, not settled by the data.
#   fPAR = (NDVI - ndvi_min) / (ndvi_max - ndvi_min), clipped to [0, 1], with
#   ndvi_min = 0.05 / ndvi_max = 0.95 as LITERATURE DEFAULTS -- these are NOT
#   empirically calibrated for these specific sites, unlike Wang et al.
#   (2016)'s ground-calibrated approach. Known limitation, worth revisiting
#   with real site data later.
#   APAR per half-hour = PPFD_IN * fPAR, since fPAR is one constant value per
#   site-year (a single NDVI snapshot), not itself time-varying within the
#   window.
#   Midday filtering uses LOCAL SOLAR time (Section 3), not local clock time
#   or timezone -- AmeriFlux TIMESTAMP is local STANDARD time (no DST), so
#   solar time is computed as clock time + a longitude-based offset from the
#   nominal 15-degree-wide standard meridian nearest the site + the equation
#   of time (Cooper's approximation). The nominal-meridian assumption is
#   itself an approximation of each site's real civil standard-time offset
#   (which can deviate from the nearest 15-degree line for political-boundary
#   reasons) -- flagged here as a known limitation, not verified per-site
#   against real UTC-offset metadata (none of which is present in
#   NEONsites.csv).
#
# Output: one row per tower-YEAR, written incrementally after each site-year
# (crash-safe), to out_csv (Section 0).
#
# VALIDATION: run in this sandbox against synthetic HH CSVs with a KNOWN true
# LUE slope (no access to real Data/ here) -- see SESSION_LOG.md for what was
# checked and what still needs a real-data run on the server, per this
# project's established working pattern.
#
# ----------------------------------------------------------------------------
# MODIS-DERIVED LUE ADDITION (same session as Code/DataAnalysis/ExtractMODIS.R):
#
# Adds a second, satellite-derived LUE estimate per site-year (`lue_modis`,
# `lue_modis_r2`) alongside the original NDVI-derived one, plus reference
# columns (`fpar_modis`, `gpp_modis_ref`, `n_modis_composites_used`), reading
# ExtractMODIS.R's `modis_fpar_raw.csv` / `modis_gpp_raw.csv` outputs
# (8-day-composite MOD15A2H Fpar_500m / MOD17A2H Gpp_500m, already QC-filtered
# by ExtractMODIS.R). The original `fpar` output column is renamed
# `ndvi_fpar` so the two fPAR sources are unambiguous in the output header.
#
# Apples-to-apples design: LUE_modis = slope of lm(GPP_NT_VUT_REF ~ APAR_modis)
# over the EXACT SAME filtered half-hourly rows (`filtered`, Section 6/7) as
# the NDVI-derived regression -- same tower GPP, same PPFD_IN, same midday/
# growing-season/QC filtering -- only APAR_modis = PPFD_IN * fpar_modis
# differs (fpar_modis in place of the NDVI-derived fPAR). fit_lue_regression()
# is reused unchanged for both fits, not duplicated.
# fpar_modis itself is the mean of QC-passed 8-day composites falling inside
# the SAME growing-season window (gs$months / target_year) already resolved
# for the NDVI-derived version (Section 4, real-flight-month vs. fallback) --
# reused via aggregate_modis_window(), not a separately reimplemented window.
# gpp_modis_ref is reported for comparison only -- it is NEVER used as GPP in
# either LUE regression; tower GPP_NT_VUT_REF remains the sole GPP input to
# both fits.
#
# MODIS inputs are OPTIONAL, not required, at this script's top-level check:
# the NDVI-derived LUE pipeline has no real dependency on MODIS data and
# should not be blocked by ExtractMODIS.R not having been run yet (it's a
# separate, network-dependent, untested-in-this-sandbox script -- see its own
# header). If `modis_fpar_raw.csv`/`modis_gpp_raw.csv` are missing,
# load_modis_csv() cat()s a clear one-time notice and every MODIS-derived
# column is NA for the whole run, rather than stop()ing the entire pipeline.
# A MODIS-side failure for one specific site-year (e.g. a corrupt row) is
# also caught per-site-year and only NAs that site-year's MODIS columns --
# it never discards an already-successful NDVI-derived result for that row.
#
# VALIDATION (MODIS addition): synthetic modis_fpar_raw.csv/modis_gpp_raw.csv
# fixtures with a KNOWN true LUE_modis slope, exercising the join/aggregation
# logic and confirming fit_lue_regression() recovers the known slope when fed
# APAR_modis instead of APAR -- see SESSION_LOG.md. Real MODIS values
# (produced by actually running ExtractMODIS.R on the server) are NOT
# validated here, same caveat as ExtractMODIS.R itself.

library(terra)
library(rhdf5)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# ============================================================================
# 0. Config
# ============================================================================
hyperspec_dir <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
# UNCONFIRMED -- best-guess sibling of NEON_FluxVariability.R's
# "./Data/NEON_Ameriflux/AnnualData"; adjust to the real HH directory before
# trusting a run that reports "no AmeriFlux HH data found" for every site.
hh_data_dir   <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Ameriflux/HHData"
out_csv       <- "D:/projects/moore/SpectralBiodiversity/Data/lue_by_year.csv"

# ExtractMODIS.R outputs -- OPTIONAL (see MODIS-DERIVED LUE ADDITION header
# note above): missing files degrade to all-NA MODIS columns, they don't stop().
modis_fpar_csv <- "D:/projects/moore/SpectralBiodiversity/Data/modis_fpar_raw.csv"
modis_gpp_csv  <- "D:/projects/moore/SpectralBiodiversity/Data/modis_gpp_raw.csv"

spectral_diversity_script <- "./Code/AnnualSpectralDiversity.R"
field_diversity_script    <- "./Code/DataAnalysis/FieldDiversity.R"

buffer_m            <- 500    # matches AnnualSpectralDiversity.R's single site-year buffer
ndvi_thresh          <- 0.4   # veg-mask threshold, matches AnnualSpectralDiversity.R
ndvi_min             <- 0.05  # fPAR literature default (NOT site-calibrated -- see header)
ndvi_max             <- 0.95  # fPAR literature default (NOT site-calibrated -- see header)
ppfd_min_threshold   <- 500   # umol m-2 s-1, midday-regression PPFD floor, adjustable
midday_solar_start   <- 10    # local solar time window (hours), adjustable
midday_solar_end     <- 14
qc_max_allowed       <- 1     # keep NEE_VUT_REF_QC <= this (0=measured,1=good gapfill);
                               # no prior repo convention -- see header
fallback_growing_season_months <- 6:9   # June-Sept, used only when no real flight date resolves
min_halfhours_for_regression   <- 10    # minimum n for a meaningful 2-parameter lm fit; judgment call

towers_df <- read.csv("./Data/NEONsites.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites in NEONsites.csv.\n")

# MODIS reference data (optional -- see header note). Loaded once, not
# per-tower like the HH flux cache, since these are small aggregated 8-day
# time series rather than large half-hourly archives.
load_modis_csv <- function(path) {
  if (!file.exists(path)) {
    cat("MODIS input not found:", path, "-- run ExtractMODIS.R on the server first.",
        "MODIS-derived LUE columns will be NA for this run.\n")
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$date <- as.Date(df$date)
  as_tibble(df)
}

modis_fpar_df <- load_modis_csv(modis_fpar_csv)
modis_gpp_df  <- load_modis_csv(modis_gpp_csv)

# ============================================================================
# 1. Import (not copy) shared helpers from other pipeline scripts
# ============================================================================
# Parses `script_path` and eval()s ONLY the named top-level `<-` bindings
# (function defs or plain config constants, e.g. bad_band_ranges) into a
# dedicated environment whose parent is the caller's environment -- so a free
# variable inside an imported function (e.g. get_tower_reflectance() reading
# the caller's `towers_df`) resolves normally via lexical scoping, without
# executing anything else in the origin script.
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

# ============================================================================
# 2. Local solar time (Section header note: see METHODOLOGY above for the
#    nominal-meridian approximation this relies on)
# ============================================================================
# Equation of time (minutes), Cooper's approximation -- standard solar-
# engineering formula, day-of-year based.
equation_of_time_min <- function(doy) {
  B <- 2 * pi * (doy - 81) / 364
  9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B)
}

# timestamp: POSIXct, treated as a naive local-standard-time wall clock (the
# tz attribute is irrelevant here -- only clock hour/minute/day-of-year are
# read off it). lon: site longitude, decimal degrees. Returns fractional
# local solar hour (e.g. 13.25 = 13:15 solar time).
local_solar_hour <- function(timestamp, lon) {
  doy <- as.numeric(format(timestamp, "%j"))
  clock_hour <- as.numeric(format(timestamp, "%H")) + as.numeric(format(timestamp, "%M")) / 60
  standard_meridian <- 15 * round(lon / 15)
  offset_min <- 4 * (lon - standard_meridian) + equation_of_time_min(doy)
  clock_hour + offset_min / 60
}

# ============================================================================
# 3. NDVI / fPAR per site-year (reuses imported spectral-diversity helpers)
# ============================================================================
compute_mean_ndvi <- function(tower_id, h5_files, buffer_m) {
  data <- spec_fns$get_tower_reflectance(tower_id, h5_files, buffer_m)
  if (is.null(data)) return(NA_real_)
  ndvi <- spec_fns$compute_ndvi_raster(data$raster, data$wavelengths)
  # terra::mask() only excludes cells where the mask is NA (default
  # maskvalues = NA) -- same fix as AnnualSpectralDiversity.R's veg_mask.
  veg_mask <- ifel(ndvi > ndvi_thresh, 1, NA)
  masked <- mask(ndvi, veg_mask)
  val <- terra::global(masked, fun = "mean", na.rm = TRUE)[1, 1]
  if (is.nan(val)) NA_real_ else as.numeric(val)
}

compute_fpar <- function(ndvi_val) {
  fpar <- (ndvi_val - ndvi_min) / (ndvi_max - ndvi_min)
  pmin(pmax(fpar, 0), 1)
}

# ============================================================================
# 4. Growing-season window resolution (real flight month vs. fallback)
# ============================================================================
resolve_growing_season <- function(h5_files) {
  if (length(h5_files) == 0) {
    return(list(months = fallback_growing_season_months, source = "fallback_default"))
  }
  flight_date <- tryCatch(field_fns$get_flight_acquisition_date(h5_files[1]),
                           error = function(e) NA)
  if (!is.na(flight_date)) {
    return(list(months = as.integer(format(flight_date, "%m")), source = "real_metadata"))
  }
  list(months = fallback_growing_season_months, source = "fallback_default")
}

# ============================================================================
# 5. AmeriFlux HH file discovery + reading
# ============================================================================
# Flexible pattern -- deliberately NOT a rigid version-number regex, since the
# real filename/version convention couldn't be confirmed (see header).
find_hh_files <- function(tower_id) {
  if (!dir.exists(hh_data_dir)) return(character(0))
  all_csv <- list.files(hh_data_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  keep <- str_detect(basename(all_csv), fixed(paste0("AMF_", tower_id))) &
    str_detect(basename(all_csv), fixed("FLUXNET")) &
    str_detect(basename(all_csv), fixed("_HH_"))
  all_csv[keep]
}

read_hh_flux <- function(path) {
  df <- read.csv(path, na.strings = "-9999", stringsAsFactors = FALSE)
  required_cols <- c("TIMESTAMP_START", "GPP_NT_VUT_REF", "PPFD_IN")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("read_hh_flux(): missing required column(s) [", paste(missing_cols, collapse = ", "),
         "] in ", path, " -- real file header differs from the assumed FLUXNET2015 HH",
         " schema; inspect the file and update column names here rather than guessing.")
  }
  df$TIMESTAMP_START <- as.POSIXct(as.character(df$TIMESTAMP_START), format = "%Y%m%d%H%M", tz = "UTC")
  df <- as_tibble(df)
  if ("NEE_VUT_REF_QC" %in% names(df)) {
    # QC value or NA both dropped unless confirmed <= qc_max_allowed -- a
    # missing QC flag is not treated as "good".
    df <- df %>% filter(NEE_VUT_REF_QC <= qc_max_allowed)
  }
  df
}

# per-tower_id cache: a site's HH file(s) can span many years and are read
# once, reused across every hyperspec year on disk for that tower, instead of
# re-reading a potentially large CSV once per site-year.
flux_cache <- new.env()
get_flux_data <- function(tower_id) {
  if (!exists(tower_id, envir = flux_cache, inherits = FALSE)) {
    files <- find_hh_files(tower_id)
    combined <- if (length(files) == 0) {
      NULL
    } else {
      dfs <- compact(map(files, possibly(read_hh_flux, otherwise = NULL)))
      if (length(dfs) == 0) NULL else {
        bind_rows(dfs) %>% distinct(TIMESTAMP_START, .keep_all = TRUE) %>% arrange(TIMESTAMP_START)
      }
    }
    assign(tower_id, combined, envir = flux_cache)
  }
  get(tower_id, envir = flux_cache, inherits = FALSE)
}

# ============================================================================
# 6. LUE regression for one site-year's filtered half-hourly records
# ============================================================================
filter_midday_growing_season <- function(flux_df, lon, target_year, months) {
  flux_df %>%
    mutate(row_year = as.integer(format(TIMESTAMP_START, "%Y")),
           month    = as.integer(format(TIMESTAMP_START, "%m")),
           solar_hr = local_solar_hour(TIMESTAMP_START, lon)) %>%
    filter(row_year == target_year,
           month %in% months,
           solar_hr >= midday_solar_start, solar_hr <= midday_solar_end,
           PPFD_IN >= ppfd_min_threshold,
           !is.na(GPP_NT_VUT_REF), !is.na(PPFD_IN))
}

fit_lue_regression <- function(filtered, fpar_val) {
  filtered$APAR <- filtered$PPFD_IN * fpar_val
  fit <- lm(GPP_NT_VUT_REF ~ APAR, data = filtered)
  s <- summary(fit)
  list(lue = unname(coef(fit)["APAR"]),
       intercept = unname(coef(fit)["(Intercept)"]),
       r2 = s$r.squared)
}

# ============================================================================
# 6a. MODIS-derived LUE addon (see MODIS-DERIVED LUE ADDITION header note)
# ============================================================================
# Mean of QC-passed 8-day composites for one tower-year within the SAME
# growing-season window (months/target_year) already resolved for the
# NDVI-derived version -- not a separately reimplemented window.
aggregate_modis_window <- function(modis_df, tower_id_val, target_year, months, value_col) {
  if (is.null(modis_df)) return(list(mean_value = NA_real_, n_used = 0L))
  d <- modis_df %>%
    filter(tower_id == tower_id_val,
           as.integer(format(date, "%Y")) == target_year,
           as.integer(format(date, "%m")) %in% months,
           qc_passed)
  if (nrow(d) == 0) return(list(mean_value = NA_real_, n_used = 0L))
  list(mean_value = mean(d[[value_col]], na.rm = TRUE), n_used = nrow(d))
}

# `filtered`: the SAME filtered half-hourly rows already used for the
# NDVI-derived fit (identical GPP/PPFD_IN/midday/growing-season/QC
# filtering) -- only the fPAR source changes, per the apples-to-apples design
# note above. Reuses fit_lue_regression() rather than duplicating the lm() call.
compute_modis_lue_addon <- function(filtered, tower_id_val, target_year, months) {
  fpar_agg <- aggregate_modis_window(modis_fpar_df, tower_id_val, target_year, months, "fpar_modis")
  gpp_agg  <- aggregate_modis_window(modis_gpp_df,  tower_id_val, target_year, months, "gpp_modis")

  if (is.na(fpar_agg$mean_value)) {
    return(list(fpar_modis = NA_real_, lue_modis = NA_real_, lue_modis_r2 = NA_real_,
                gpp_modis_ref = gpp_agg$mean_value, n_modis_composites_used = fpar_agg$n_used))
  }
  fit_modis <- fit_lue_regression(filtered, fpar_agg$mean_value)
  list(fpar_modis = fpar_agg$mean_value, lue_modis = fit_modis$lue, lue_modis_r2 = fit_modis$r2,
       gpp_modis_ref = gpp_agg$mean_value, n_modis_composites_used = fpar_agg$n_used)
}

# ============================================================================
# 7. Main loop: one row per tower-YEAR
# ============================================================================
make_row <- function(tower_id, neon_site, year, lue = NA_real_, lue_intercept = NA_real_,
                      lue_r2 = NA_real_, ndvi = NA_real_, ndvi_fpar = NA_real_,
                      fpar_modis = NA_real_, lue_modis = NA_real_, lue_modis_r2 = NA_real_,
                      gpp_modis_ref = NA_real_, n_modis_composites_used = NA_integer_,
                      n_halfhours_used = NA_integer_, growing_season_source = NA_character_,
                      status) {
  tibble(tower_id = tower_id, neon_site = neon_site, year = year,
         lue = lue, lue_intercept = lue_intercept, lue_r2 = lue_r2,
         ndvi = ndvi, ndvi_fpar = ndvi_fpar,
         fpar_modis = fpar_modis, lue_modis = lue_modis, lue_modis_r2 = lue_modis_r2,
         gpp_modis_ref = gpp_modis_ref,
         n_modis_composites_used = as.integer(n_modis_composites_used),
         n_halfhours_used = as.integer(n_halfhours_used),
         growing_season_source = growing_season_source, status = status)
}

run_lue_pipeline <- function() {
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
                     lue = double(), lue_intercept = double(), lue_r2 = double(),
                     ndvi = double(), ndvi_fpar = double(),
                     fpar_modis = double(), lue_modis = double(), lue_modis_r2 = double(),
                     gpp_modis_ref = double(), n_modis_composites_used = integer(),
                     n_halfhours_used = integer(),
                     growing_season_source = character(), status = character())

  pb <- txtProgressBar(min = 0, max = length(site_year_jobs), style = 3)

  for (j in seq_along(site_year_jobs)) {
    job <- site_year_jobs[[j]]
    tower_id  <- job$tower_id
    yr        <- job$year
    neon_site <- towers_df$neon_site[towers_df$Site.ID == tower_id][1]
    lon       <- towers_df$Lon[towers_df$Site.ID == tower_id][1]

    cat("\n==== [", j, "/", length(site_year_jobs), "] ", tower_id, " - ", yr, " ====\n", sep = "")

    row <- tryCatch({
      ndvi_val <- compute_mean_ndvi(tower_id, job$files, buffer_m)

      if (is.na(ndvi_val)) {
        make_row(tower_id, neon_site, yr, status = "no vegetation pixels for NDVI")
      } else {
        fpar_val <- compute_fpar(ndvi_val)
        gs <- resolve_growing_season(job$files)
        flux_df <- get_flux_data(tower_id)

        if (is.null(flux_df)) {
          make_row(tower_id, neon_site, yr, ndvi = ndvi_val, ndvi_fpar = fpar_val,
                    growing_season_source = gs$source, status = "no AmeriFlux HH data found")
        } else {
          target_year <- suppressWarnings(as.integer(yr))
          filtered <- filter_midday_growing_season(flux_df, lon, target_year, gs$months)

          if (nrow(filtered) < min_halfhours_for_regression) {
            make_row(tower_id, neon_site, yr, ndvi = ndvi_val, ndvi_fpar = fpar_val,
                      n_halfhours_used = nrow(filtered), growing_season_source = gs$source,
                      status = paste0("insufficient half-hourly data after filtering (n=",
                                       nrow(filtered), ")"))
          } else {
            fit <- fit_lue_regression(filtered, fpar_val)

            # MODIS addon failing (missing inputs, a bad row, etc.) must never
            # discard an already-successful NDVI-derived result -- caught and
            # NA'd on its own, separate from the outer tryCatch.
            modis_add <- tryCatch(
              compute_modis_lue_addon(filtered, tower_id, target_year, gs$months),
              error = function(e) {
                cat("  MODIS-derived LUE failed for this site-year:", conditionMessage(e),
                    "-- NDVI-derived result unaffected.\n")
                list(fpar_modis = NA_real_, lue_modis = NA_real_, lue_modis_r2 = NA_real_,
                     gpp_modis_ref = NA_real_, n_modis_composites_used = 0L)
              }
            )

            make_row(tower_id, neon_site, yr, lue = fit$lue, lue_intercept = fit$intercept,
                      lue_r2 = fit$r2, ndvi = ndvi_val, ndvi_fpar = fpar_val,
                      fpar_modis = modis_add$fpar_modis, lue_modis = modis_add$lue_modis,
                      lue_modis_r2 = modis_add$lue_modis_r2, gpp_modis_ref = modis_add$gpp_modis_ref,
                      n_modis_composites_used = modis_add$n_modis_composites_used,
                      n_halfhours_used = nrow(filtered), growing_season_source = gs$source,
                      status = "success")
          }
        }
      }
    }, error = function(e) {
      cat("  unexpected error:", conditionMessage(e), "\n")
      make_row(tower_id, neon_site, yr, status = paste("error:", conditionMessage(e)))
    })

    results <- bind_rows(results, row)
    # save incrementally so a crash/interruption partway through doesn't lose
    # completed site-years (matches AnnualSpectralDiversity.R / AnnualFootprintExtent.R)
    write.csv(results, out_csv, row.names = FALSE)

    cat("  status:", row$status, "\n")
    setTxtProgressBar(pb, j)
  }

  close(pb)
  cat("\nDone. Summary written to:", out_csv, "\n")
  results
}

results <- run_lue_pipeline()
