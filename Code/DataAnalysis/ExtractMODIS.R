# ExtractMODIS.R
#
# Pulls MOD15A2H FPAR (Fpar_500m) and MOD17A2H GPP (Gpp_500m), 8-day
# composite, at each of the 45 NEON/AmeriFlux towers' Lat/Lon, via NASA's
# AppEEARS point-sample API. Independent, satellite-derived cross-check
# against the NDVI-derived fPAR and tower GPP used in ComputeLUE.R -- this
# script only DOWNLOADS and QC-filters the raw 8-day time series; aggregation
# to one value per site-year and the MODIS-fPAR LUE comparison happen in
# ComputeLUE.R (Section on MODIS-derived LUE, added alongside the original
# NDVI-derived pipeline in this same session).
#
# ==============================================================================
# CRITICAL CAVEAT, READ BEFORE RUNNING: the AppEEARS submit/poll/download
# section of this script (everything that makes an HTTP request) was written
# and NEVER executed -- this sandbox has no internet access, so it could not
# call NASA's login, task-submission, task-status, or bundle-download
# endpoints even once. Every endpoint URL, JSON request-body shape, response
# field name (task_id, status values, bundle file_id/file_name), and the
# real AppEEARS point-CSV column-naming convention are written from
# documented AppEEARS API/product conventions, NOT verified against a live
# response. Per this project's established real-data-validation discipline
# (see SESSION_LOG.md -- e.g. FieldDiversity.R's H5-metadata investigation,
# AnnualSpectralDiversity.R's synthetic-only validation), UNTESTED CODE IS
# NOT CLAIMED TO WORK. You will need to run this on the server, watch it fail
# at whatever point the real API actually differs from what's assumed here,
# and fix that point specifically -- start with `cat()`-ing intermediate
# objects (the login response, the task-submit response, `names(bundle)`,
# `names(raw_combined)`) rather than assuming the whole thing is broken if
# one step doesn't match.
#
# What IS validated in this sandbox (synthetic data, see SESSION_LOG.md):
# the QC-decoding bit logic, the generic column-matching parser
# (parse_appeears_csv / find_col / find_col_exact), and the final
# tower_id/neon_site/date/value/qc_flag/qc_passed output shape -- i.e.
# everything downstream of "a CSV that looks like AppEEARS' documented point-
# sample output exists on disk". That parsing logic is written to be robust
# to small naming differences (see find_col()'s substring match rather than
# an exact hardcoded column name) specifically because the real column names
# could not be confirmed.
# ==============================================================================
#
# CREDENTIALS: AppEEARS requires an Earthdata Login (free account at
# https://urs.earthdata.nasa.gov/). Read from EARTHDATA_USERNAME /
# EARTHDATA_PASSWORD environment variables -- NEVER hardcoded here. Set them
# in your ~/.Renviron (already gitignored in this repo -- confirmed against
# .gitignore before choosing this: `.Renviron` is listed alongside
# `.Rproj.user`/`.Rhistory`/`.RData`/`.Ruserdata`, so no new gitignore entry
# was needed) or in the shell environment before running this script.
#
# PRODUCTS/LAYERS (collection .061, the current MODIS collection as of this
# writing -- confirm this is still current when you run this):
#   MOD15A2H.061  Fpar_500m       -- primary FPAR layer
#   MOD15A2H.061  FparLai_QC      -- its QC band (bit-packed, see decode_fpar_qc())
#   MOD17A2H.061  Gpp_500m        -- secondary GPP cross-check (NOT used as
#                                    primary GPP anywhere -- tower
#                                    GPP_NT_VUT_REF remains primary in
#                                    ComputeLUE.R)
#   MOD17A2H.061  Psn_QC_500m     -- its QC band (see decode_gpp_qc())
# MOD15A2H (Terra-only) is used rather than the combined MCD15A2H per the
# task's stated preference order; swap `fpar_product` below to
# "MCD15A2H.061" for Terra+Aqua combined (denser 8-day compositing, shorter
# period of record) if preferred later -- not done here since it wasn't the
# first-listed option and the tradeoff wasn't discussed.
#
# Submitted as ONE combined AppEEARS point task covering all 4 layers and all
# 45 site coordinates at once (AppEEARS point tasks accept multiple
# products/layers and multiple coordinates in a single task), not one HTTP
# request per site or per product -- matches the "batch point request, not a
# loop of single-site calls" instruction. Output is still split into the two
# requested files (modis_fpar_raw.csv, modis_gpp_raw.csv) after parsing.
#
# Date range is NOT hardcoded -- read from spectral_diversity_by_year.csv
# (AnnualSpectralDiversity.R's output) so it always covers whatever
# site-years actually exist, whenever this is re-run.

library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(httr2)

# ============================================================================
# 0. Config
# ============================================================================
spectral_diversity_csv <- "D:/projects/moore/SpectralBiodiversity/Data/spectral_diversity_by_year.csv"
appeears_download_dir  <- "D:/projects/moore/SpectralBiodiversity/Data/AppEEARS_raw"
fpar_out_csv            <- "D:/projects/moore/SpectralBiodiversity/Data/modis_fpar_raw.csv"
gpp_out_csv             <- "D:/projects/moore/SpectralBiodiversity/Data/modis_gpp_raw.csv"

appeears_base_url <- "https://appeears.earthdatacloud.nasa.gov/api"

fpar_product  <- "MOD15A2H.061"
fpar_layer    <- "Fpar_500m"
fpar_qc_layer <- "FparLai_QC"
fpar_scale_factor <- 0.01    # MOD15A2H documented scale factor -- UNCONFIRMED against a real file, see caveat
fpar_valid_range  <- c(0, 100)   # raw (pre-scale) valid range; outside this treated as fill regardless of QC

gpp_product  <- "MOD17A2H.061"
gpp_layer    <- "Gpp_500m"
gpp_qc_layer <- "Psn_QC_500m"
gpp_scale_factor <- 0.0001   # MOD17A2H documented scale factor -- UNCONFIRMED against a real file, see caveat
gpp_valid_range  <- c(0, 30000)  # raw (pre-scale) valid range; outside this treated as fill regardless of QC

poll_interval_sec <- 30      # how often to re-check task status
max_wait_min       <- 60     # give up (not error out silently) after this long; point tasks are usually much faster

towers_df <- read.csv("./Data/NEONsites.csv", fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites in NEONsites.csv.\n")

coord_xwalk <- towers_df %>% transmute(tower_id = Site.ID, neon_site)

# ============================================================================
# 1. Date range -- derived from spectral_diversity_by_year.csv, never hardcoded
# ============================================================================
if (!file.exists(spectral_diversity_csv)) {
  stop("Required input not found: ", spectral_diversity_csv,
       " -- run AnnualSpectralDiversity.R first. ExtractMODIS.R needs its",
       " site-year range to know what date span to request from AppEEARS,",
       " rather than a hardcoded year list.")
}
spec_div <- read.csv(spectral_diversity_csv)
years_present <- sort(unique(suppressWarnings(as.integer(spec_div$year))))
years_present <- years_present[!is.na(years_present)]
if (length(years_present) == 0) {
  stop("Could not parse any usable year from ", spectral_diversity_csv, "$year -- inspect the file.")
}

# AppEEARS point-task date format, per documented API convention: MM-DD-YYYY
# -- UNCONFIRMED against a live response, see header caveat.
appeears_start_date <- sprintf("01-01-%d", min(years_present))
appeears_end_date   <- sprintf("12-31-%d", max(years_present))
cat("Site-years span", min(years_present), "-", max(years_present),
    "(from spectral_diversity_by_year.csv) -- requesting MODIS data for",
    appeears_start_date, "to", appeears_end_date, "\n")

# ============================================================================
# 2. QC decoding (documented MOD15A2H/MOD17A2H bit-packed QC convention --
#    stable, published product-spec bit layout, independent of the API
#    transport -- but still not verified against a real downloaded value
#    since no live file was available. Confirm on first real run.)
# ============================================================================
# MOD15A2H FparLai_QC, 8-bit:
#   bit 0       MODLAND_QC   0 = good quality, 1 = other quality
#   bits 5-7    CloudState   000 = clear, anything else = some cloud influence
# Kept only when BOTH good MODLAND quality AND clear-sky.
decode_fpar_qc <- function(qc) {
  qc <- as.integer(qc)
  modland_good <- bitwAnd(qc, 1L) == 0L
  cloud_state  <- bitwAnd(bitwShiftR(qc, 5L), 7L)
  clear_sky    <- cloud_state == 0L
  modland_good & clear_sky & !is.na(qc)
}

# MOD17A2H Psn_QC_500m, 8-bit -- same MODLAND_QC bit-0 convention as FPAR/LAI
# products; a fuller multi-bit interpretation (analogous to FparLai_QC's
# cloud-state bits) may exist for this product's QC byte too, but wasn't
# confirmed, so only the documented bit-0 check is applied here. Revisit if a
# real Psn_QC_500m value doesn't behave as expected.
decode_gpp_qc <- function(qc) {
  qc <- as.integer(qc)
  bitwAnd(qc, 1L) == 0L & !is.na(qc)
}

# ============================================================================
# 3. Generic AppEEARS point-CSV parser (VALIDATED against a synthetic fixture
#    -- see SESSION_LOG.md)
# ============================================================================
# Matches a column by requiring the layer name as a SUBSTRING of the column
# name, rather than an exact hardcoded "<Product>_<Version>_<Layer>" string
# -- deliberately tolerant of whatever exact prefixing convention the real
# AppEEARS output uses, since it could not be confirmed here.
find_col <- function(df, layer_name) {
  matches <- names(df)[str_detect(names(df), fixed(layer_name))]
  if (length(matches) == 0) {
    stop("Could not find a column containing '", layer_name, "' in the AppEEARS",
         " results. Real column-naming convention differs from assumed --",
         " inspect names(raw_combined) and adjust. Available columns: ",
         paste(names(df), collapse = ", "))
  }
  if (length(matches) > 1) {
    warning("Multiple columns match '", layer_name, "': ", paste(matches, collapse = ", "),
            " -- using '", matches[1], "'. Confirm this is the right one.")
  }
  matches[1]
}

find_col_exact <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop("Could not find any of [", paste(candidates, collapse = ", "),
         "] among the AppEEARS results columns: ", paste(names(df), collapse = ", "))
  }
  hit[1]
}

parse_appeears_csv <- function(raw_df, coord_xwalk, value_layer, qc_layer, qc_decoder,
                                scale_factor = 1, valid_range = c(-Inf, Inf)) {
  # AppEEARS echoes back the "id" submitted per coordinate (tower_id here) --
  # assumed column name "ID"; "Date" assumed for the composite date.
  id_col    <- find_col_exact(raw_df, c("ID", "id"))
  date_col  <- find_col_exact(raw_df, c("Date", "date"))
  value_col <- find_col(raw_df, value_layer)
  qc_col    <- find_col(raw_df, qc_layer)

  raw_df %>%
    transmute(
      tower_id  = as.character(.data[[id_col]]),
      date      = as.Date(.data[[date_col]]),
      raw_value = suppressWarnings(as.numeric(.data[[value_col]])),
      qc_flag   = suppressWarnings(as.integer(.data[[qc_col]]))
    ) %>%
    mutate(
      qc_passed = qc_decoder(qc_flag) & !is.na(raw_value) &
        raw_value >= valid_range[1] & raw_value <= valid_range[2],
      value = ifelse(qc_passed, raw_value * scale_factor, NA_real_)
    ) %>%
    inner_join(coord_xwalk, by = "tower_id") %>%
    select(tower_id, neon_site, date, value, qc_flag, qc_passed)
}

# ============================================================================
# 4. AppEEARS API calls (UNTESTED -- see header caveat)
# ============================================================================
appeears_login <- function(username, password) {
  resp <- request(paste0(appeears_base_url, "/login")) %>%
    req_method("POST") %>%
    req_auth_basic(username, password) %>%
    req_retry(max_tries = 3) %>%
    req_perform()
  token <- resp_body_json(resp)$token
  if (is.null(token)) stop("AppEEARS login response did not contain a 'token' field -- ",
                            "response shape differs from assumed; inspect resp_body_json(resp).")
  token
}

build_task_body <- function(task_name, start_date, end_date, coord_xwalk_full) {
  layers <- list(
    list(product = fpar_product, layer = fpar_layer),
    list(product = fpar_product, layer = fpar_qc_layer),
    list(product = gpp_product,  layer = gpp_layer),
    list(product = gpp_product,  layer = gpp_qc_layer)
  )
  coordinates <- pmap(
    list(coord_xwalk_full$Site.ID, coord_xwalk_full$Lat, coord_xwalk_full$Lon),
    function(id, lat, lon) list(id = id, latitude = lat, longitude = lon, category = "site")
  )
  list(
    task_type = "point",
    task_name = task_name,
    params = list(
      dates       = list(list(startDate = start_date, endDate = end_date)),
      layers      = layers,
      coordinates = coordinates
    )
  )
}

submit_appeears_task <- function(token, task_body) {
  resp <- request(paste0(appeears_base_url, "/task")) %>%
    req_auth_bearer_token(token) %>%
    req_body_json(task_body) %>%
    req_retry(max_tries = 3) %>%
    req_perform()
  task_id <- resp_body_json(resp)$task_id
  if (is.null(task_id)) stop("AppEEARS task-submit response did not contain a 'task_id' field -- ",
                              "response shape differs from assumed; inspect resp_body_json(resp).")
  task_id
}

poll_appeears_task <- function(token, task_id) {
  deadline <- Sys.time() + max_wait_min * 60
  repeat {
    resp <- request(paste0(appeears_base_url, "/task/", task_id)) %>%
      req_auth_bearer_token(token) %>%
      req_retry(max_tries = 3) %>%
      req_perform()
    status <- resp_body_json(resp)$status
    cat("  AppEEARS task", task_id, "status:", status, "\n")
    if (identical(status, "done")) return(invisible(TRUE))
    if (status %in% c("error", "expired")) {
      stop("AppEEARS task ", task_id, " ended with status '", status, "' -- check ",
           "https://appeears.earthdatacloud.nasa.gov/view/", task_id, " for details.")
    }
    if (Sys.time() > deadline) {
      stop("AppEEARS task ", task_id, " did not complete within ", max_wait_min,
           " minutes -- it may still be processing; check ",
           "https://appeears.earthdatacloud.nasa.gov/view/", task_id,
           " and re-run download_appeears_results() directly once done, rather",
           " than assuming failure.")
    }
    Sys.sleep(poll_interval_sec)
  }
}

download_appeears_results <- function(token, task_id, dest_dir) {
  resp <- request(paste0(appeears_base_url, "/bundle/", task_id)) %>%
    req_auth_bearer_token(token) %>%
    req_retry(max_tries = 3) %>%
    req_perform()
  bundle_files <- resp_body_json(resp)$files
  if (is.null(bundle_files) || length(bundle_files) == 0) {
    stop("AppEEARS bundle listing for task ", task_id, " was empty or in an",
         " unexpected shape -- inspect resp_body_json(resp) manually.")
  }
  # The per-point results CSV vs. metadata/granule-list CSVs in the bundle
  # are distinguished here by name pattern -- UNCONFIRMED, inspect
  # `bundle_files` on first real run and adjust this filter if it picks up
  # the wrong file(s) or misses the real one.
  is_results_csv <- function(f) {
    grepl("\\.csv$", f$file_name, ignore.case = TRUE) &&
      !grepl("granule-list|metadata|log", f$file_name, ignore.case = TRUE)
  }
  csv_entries <- keep(bundle_files, is_results_csv)
  if (length(csv_entries) == 0) {
    stop("No results CSV identified in the AppEEARS bundle for task ", task_id,
         " -- file names present: ", paste(map_chr(bundle_files, "file_name"), collapse = ", "),
         ". Adjust is_results_csv() to match the real naming.")
  }

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  map_chr(csv_entries, function(entry) {
    out_path <- file.path(dest_dir, entry$file_name)
    file_resp <- request(paste0(appeears_base_url, "/bundle/", task_id, "/", entry$file_id)) %>%
      req_auth_bearer_token(token) %>%
      req_retry(max_tries = 3) %>%
      req_perform()
    writeBin(resp_body_raw(file_resp), out_path)
    cat("  downloaded", out_path, "\n")
    out_path
  })
}

# ============================================================================
# 5. Main orchestration (UNTESTED end-to-end -- see header caveat)
# ============================================================================
run_extract_modis <- function() {
  earthdata_username <- Sys.getenv("EARTHDATA_USERNAME")
  earthdata_password <- Sys.getenv("EARTHDATA_PASSWORD")
  if (!nzchar(earthdata_username) || !nzchar(earthdata_password)) {
    stop("EARTHDATA_USERNAME and EARTHDATA_PASSWORD must be set (e.g. in ~/.Renviron,",
         " already gitignored in this repo) -- Earthdata Login credentials, never",
         " hardcoded. Register at https://urs.earthdata.nasa.gov/ if needed.")
  }

  cat("Logging in to AppEEARS...\n")
  token <- appeears_login(earthdata_username, earthdata_password)

  task_name <- paste0("SpectralBiodiversity_MODIS_", min(years_present), "_", max(years_present))
  task_body <- build_task_body(task_name, appeears_start_date, appeears_end_date, towers_df)

  cat("Submitting one combined point task for", nrow(towers_df), "sites,",
      "4 layers (", fpar_layer, ",", fpar_qc_layer, ",", gpp_layer, ",", gpp_qc_layer, ")...\n")
  task_id <- submit_appeears_task(token, task_body)
  cat("Submitted task_id:", task_id, "\n")

  cat("Polling for completion (every", poll_interval_sec, "sec, up to", max_wait_min, "min)...\n")
  poll_appeears_task(token, task_id)

  cat("Downloading results bundle...\n")
  raw_csv_paths <- download_appeears_results(token, task_id, appeears_download_dir)

  raw_combined <- map_dfr(raw_csv_paths, read.csv, stringsAsFactors = FALSE)
  cat("Combined raw AppEEARS results:", nrow(raw_combined), "rows,",
      ncol(raw_combined), "columns:", paste(names(raw_combined), collapse = ", "), "\n")

  fpar_result <- parse_appeears_csv(raw_combined, coord_xwalk, fpar_layer, fpar_qc_layer,
                                     decode_fpar_qc, fpar_scale_factor, fpar_valid_range) %>%
    rename(fpar_modis = value)
  write.csv(fpar_result, fpar_out_csv, row.names = FALSE)
  cat("Wrote", nrow(fpar_result), "rows to", fpar_out_csv,
      "(", sum(fpar_result$qc_passed), "QC-passed )\n")

  gpp_result <- parse_appeears_csv(raw_combined, coord_xwalk, gpp_layer, gpp_qc_layer,
                                    decode_gpp_qc, gpp_scale_factor, gpp_valid_range) %>%
    rename(gpp_modis = value)
  write.csv(gpp_result, gpp_out_csv, row.names = FALSE)
  cat("Wrote", nrow(gpp_result), "rows to", gpp_out_csv,
      "(", sum(gpp_result$qc_passed), "QC-passed )\n")

  list(fpar = fpar_result, gpp = gpp_result)
}

results_modis <- run_extract_modis()
