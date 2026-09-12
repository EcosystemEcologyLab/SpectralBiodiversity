# FieldDiversity.R
#
# Field-measured (ground-plot) floristic diversity, computed independently of
# the hyperspectral pipelines (ComputeSpecBiodiv.R / AnnualSpectralDiversity.R)
# but using their same hyperspectral tile archive as a source of real flight
# acquisition dates, purely for temporal matching -- this script never reads
# reflectance values, only H5 metadata.
#
# Two axes are varied independently, giving up to 4 output "modes" per
# site-year:
#   plot_scope:      "tower" (tower plots only) vs. "all" (tower + distributed)
#   temporal_scope:   "per_bout" (every field bout, no selection) vs.
#                     "peak_flight" (the one bout nearest the AOP flight date,
#                     only where a flight date exists and a bout falls within
#                     30 days of it)
#
# Metrics per (tower_id, year, plot_scope, temporal_scope, bout):
#   1. floristic_shannon_mean / floristic_shannon_gamma -- abundance-weighted
#      Hill-Shannon (q=1) from the six-subplot percent-cover data, via hillR.
#   2. floristic_richness -- gamma (site-level) presence richness, unioning
#      the 1m2 percent-cover taxa, the nested 10/100m2 presence taxa, and the
#      semicolon-delimited additionalSpecies taxa.
#
# ============================================================================
# STEP 1 FINDING (H5 metadata investigation, see Section 2 below for the code
# that reproduces this at runtime against whatever tile it actually finds):
#
# Confirmed against a real DP3.30006 reflectance tile (NEON_D16_ABBY_DP3_...
# _reflectance.h5, ABBY, 2017) via h5ls(). There is NO acquisition-date
# attribute anywhere on the file, site group, /Reflectance group, or
# /Reflectance/Reflectance_Data dataset -- only sensor/geometry/scale-factor
# metadata lives at those levels. The flight date is NOT encoded in the
# Metadata/Logs/<flightline> group NAME either -- that name is just the
# flight line's HHMMSS time-of-day (e.g. "183828"), no date.
#
# The date is recoverable only from the TEXT CONTENT of the
# Metadata/Logs/<flightline>/ATCOR_Input_file string dataset, two ways:
#   (a) a "<YYYYMMDD>_<flightline_id>" token embedded in the flight-line
#       identifier string it contains (e.g. "..._20170622_183828"), OR
#   (b) its own literal first line, "22\06\2017          Date (dd/mm/year)"
#       -- note the delimiter is a literal backslash, not a slash.
# Both were verified to independently resolve to the same date (2017-06-22)
# on the real file. get_flight_acquisition_date() below tries (a) first (tied
# to the exact flightline_id opened, so it can't grab an unrelated date token
# elsewhere in the text) and falls back to (b); it stop()s with a clear
# message if neither parses, rather than guessing.
#
# All flight-line subgroups within one tile were the same date (one flight
# day, several passes) -- consistent with "day-to-day variation within a
# flight campaign doesn't matter", so reading just the first flightline of
# the first tile per tower-year is sufficient.
# ============================================================================

# Needed on this server for terra/sf's PROJ lookup; set defensively even
# though this script itself never touches a raster (rhdf5/tabular only).
Sys.setenv(PROJ_LIB = "C:/Program Files/R/R-4.4.1/library/terra/proj")

library(rhdf5)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(hillR)
library(readr)

# ============================================================================
# 0. Setup / paths
# ============================================================================
div_1m2_path         <- "./Data/NEON_FieldData/div_1m2Data.csv"
div_nested_path       <- "./Data/NEON_FieldData/div_10m2Data100m2Data.csv"
neonsites_path        <- "./Data/NEONsites.csv"
hyperspec_dir         <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
out_csv               <- "D:/projects/moore/SpectralBiodiversity/Data/field_diversity_long.csv"

# NEON DP1.10098.001 Vegetation Structure -- used only to classify species as
# canopy-exposed vs. understory for the *_canopy metrics (Section 4b). Exact
# filenames TBD once downloaded; these are the expected basic-package table
# names per NEON's own convention.
vst_apparent_path    <- "./Data/NEON_FieldData/vst_apparentindividual.csv"
vst_mapping_path     <- "./Data/NEON_FieldData/vst_mappingandtagging.csv"

# Herbaceous/non-woody companion table -- vst_apparentindividual only covers
# trees/shrubs, so without this the *_canopy metrics can only ever exclude
# confirmed-shaded woody species, never understory herbs/forbs/graminoids
# (see Section 4c). Structure unconfirmed -- investigated at runtime there.
vst_nonwoody_path    <- "./Data/NEON_FieldData/vst_non-woody.csv"

flight_match_tolerance_days <- 30

# ---- CANOPY DESIGN DECISION, flagged for review (see Section 4b) ----
# How species that never appear in vst_ at all -- vst_ only measures
# trees/shrubs above a size threshold, so many herbs/forbs/graminoids will
# never be in it -- are treated in the canopy-filtered metrics.
# Default here is "include": "unmeasured" reflects a genuine data gap in
# vst_'s sampling design (it doesn't sample below its size threshold), not
# confirmed evidence that a species is understory/non-canopy. Defaulting to
# exclude would silently discard a potentially large fraction of real
# species from floristic_richness_canopy on that inference alone, rather
# than on a measurement. "exclude" remains available as a stricter
# alternative -- toggle and re-run rather than treating either as settled.
unmeasured_species_treatment <- "include"  # "include" or "exclude" -- REVIEW

required_inputs <- c(div_1m2_path, div_nested_path, neonsites_path,
                      vst_apparent_path, vst_mapping_path, vst_nonwoody_path)
missing_inputs  <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Required input(s) not found:\n  ", paste(missing_inputs, collapse = "\n  "),
       "\n\nThe vst_ paths are NEON DP1.10098.001 (Vegetation Structure) basic-package",
       " tables, needed for the canopy-filtered floristic_*_canopy metrics. Download",
       " DP1.10098.001 for the relevant site(s)/year(s) and place",
       " vst_apparentindividual.csv and vst_mappingandtagging.csv under",
       " ./Data/NEON_FieldData/ (same convention as div_1m2Data.csv). Do not fabricate",
       " or approximate this data -- if it is genuinely not available yet, the",
       " canopy metrics cannot be computed.")
}
if (!dir.exists(hyperspec_dir)) {
  stop("Hyperspectral tile directory not found: ", hyperspec_dir)
}

# ============================================================================
# 1. Site crosswalk (siteID <-> tower_id), same convention as
#    ComputeSpecBiodiv.R / AnnualSpectralDiversity.R
# ============================================================================
towers_df <- read.csv(neonsites_path, fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

site_xwalk <- towers_df %>% transmute(tower_id = Site.ID, neon_site)

cat("Found", nrow(site_xwalk), "sites in NEONsites.csv.\n")

# ============================================================================
# 2. H5 metadata helpers -- flight acquisition date (see STEP 1 note above)
# ============================================================================
# byTileAOP() writes tiles under nested product/year/FullSite/domain/
# <year>_<SITECODE>_<visit>/... directories -- year is read off that folder
# name, not the h5 filename (same convention as AnnualSpectralDiversity.R's
# extract_year_from_path()/discover_site_years(), reused here verbatim).
extract_year_from_path <- function(path) {
  parts <- str_split(path, "[/\\\\]")[[1]]
  m <- str_match(parts, "^(\\d{4})_[A-Za-z]{4}_\\d+$")
  yr <- m[!is.na(m[, 2]), 2]
  if (length(yr) >= 1) return(yr[1])
  yr2 <- parts[str_detect(parts, "^(19|20)\\d{2}$")]
  if (length(yr2) >= 1) return(yr2[1])
  NA_character_
}

discover_site_year_tiles <- function(tower_id) {
  tile_dir <- file.path(hyperspec_dir, tower_id)
  if (!dir.exists(tile_dir)) return(tibble(file = character(), year = character()))
  h5_files <- list.files(tile_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(h5_files) == 0) return(tibble(file = character(), year = character()))
  tibble(file = h5_files, year = map_chr(h5_files, extract_year_from_path)) %>%
    filter(!is.na(year))
}

# Reads the acquisition date from the FIRST flight-line log group of a given
# tile -- one tile is sufficient per tower-year (see STEP 1 note); this does
# not aggregate/median across tiles or flightlines.
get_flight_acquisition_date <- function(h5_path) {
  fid <- H5Fopen(h5_path)
  on.exit(H5Fclose(fid), add = TRUE)

  listing <- h5ls(fid)
  refl_group <- listing %>% filter(str_detect(name, "^Reflectance$")) %>% pull(group) %>% unique()
  if (length(refl_group) == 0) stop("No Reflectance group found in ", h5_path)
  site_group <- str_split(refl_group[1], "/")[[1]][2]

  logs_path <- paste0("/", site_group, "/Reflectance/Metadata/Logs")
  flightlines <- listing %>% filter(group == logs_path) %>% pull(name)
  if (length(flightlines) == 0) {
    stop("No flight-line log groups found under ", logs_path, " in ", h5_path)
  }
  flightline_id <- flightlines[1]

  input_txt <- h5read(fid, paste0(logs_path, "/", flightline_id, "/ATCOR_Input_file"))

  # (a) primary: "<YYYYMMDD>_<flightline_id>" token tied to this exact
  # flightline, e.g. "..._20170622_183828" when flightline_id == "183828"
  date_match <- str_match(input_txt, paste0("(\\d{8})_", flightline_id))
  if (!is.na(date_match[1, 2])) {
    return(as.Date(date_match[1, 2], format = "%Y%m%d"))
  }

  # (b) fallback: the file's own literal header line,
  # "22\06\2017          Date (dd/mm/year)" (backslash-delimited, not "/")
  date_line <- str_match(input_txt, "(\\d{2})\\\\(\\d{2})\\\\(\\d{4})\\s+Date \\(dd/mm/year\\)")
  if (!is.na(date_line[1, 2])) {
    return(as.Date(paste(date_line[1, 4], date_line[1, 3], date_line[1, 2], sep = "-")))
  }

  stop("Could not parse flight acquisition date from ATCOR_Input_file in ", h5_path,
       " (flightline ", flightline_id, ") -- structure differs from the confirmed",
       " convention documented at the top of this script; do not guess, investigate first.")
}

# ---- run the STEP 1 investigation live, against whatever tile is actually
# on disk right now, so its real structure is visible in the console before
# any date-matching logic built on it is trusted. ----
cat("\n==== STEP 1: H5 metadata structure investigation ====\n")
probe_tiles <- list.files(hyperspec_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
if (length(probe_tiles) == 0) {
  cat("No hyperspectral tiles found under", hyperspec_dir, "-- skipping investigation.\n")
} else {
  probe_file <- probe_tiles[1]
  cat("Probe file:", probe_file, "\n\n")
  probe_listing <- h5ls(probe_file, all = TRUE)
  meta_rows <- probe_listing %>% filter(str_detect(group, "Metadata"))
  print(as.data.frame(meta_rows[, c("group", "name", "otype", "num_attrs")]), row.names = FALSE)

  probe_date <- tryCatch(get_flight_acquisition_date(probe_file),
                          error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NA })
  cat("\nResolved acquisition date:", format(probe_date), "\n")
  cat("(No date attribute exists at file/site/Reflectance/Reflectance_Data level --\n",
      " date is recovered from the ATCOR_Input_file text blob under\n",
      " Metadata/Logs/<flightline>/, per the STEP 1 note above.)\n", sep = "")
}

# ============================================================================
# 3. STEP 2 -- one representative flight date per tower_id-year on disk
# ============================================================================
get_representative_flight_dates <- function() {
  all_tower_ids <- unique(site_xwalk$tower_id)
  rows <- list()
  for (tid in all_tower_ids) {
    tiles <- discover_site_year_tiles(tid)
    if (nrow(tiles) == 0) next
    for (yr in sort(unique(tiles$year))) {
      first_file <- tiles$file[tiles$year == yr][1]
      fdate <- tryCatch(get_flight_acquisition_date(first_file), error = function(e) {
        cat("  [flight date] FAILED for", tid, yr, ":", conditionMessage(e), "\n")
        as.Date(NA)
      })
      rows[[length(rows) + 1]] <- tibble(tower_id = tid, year = as.integer(yr),
                                          flight_date = fdate, source_file = first_file)
    }
  }
  if (length(rows) == 0) return(tibble(tower_id = character(), year = integer(),
                                        flight_date = as.Date(character()), source_file = character()))
  bind_rows(rows)
}

cat("\n==== STEP 2: resolving representative flight dates per tower-year ====\n")
flight_dates <- get_representative_flight_dates()
cat("Resolved", sum(!is.na(flight_dates$flight_date)), "of", nrow(flight_dates), "tower-year flight dates.\n")

# ============================================================================
# 4. Field data -- load, filter, crosswalk
# ============================================================================
div_1m2 <- read.csv(div_1m2_path, fileEncoding = "UTF-8-BOM") %>%
  filter(divDataType == "plantSpecies") %>%
  mutate(endDate = as.Date(substr(endDate, 1, 10)),
         year = as.integer(format(endDate, "%Y"))) %>%
  inner_join(site_xwalk, by = c("siteID" = "neon_site"))

div_nested <- read.csv(div_nested_path, fileEncoding = "UTF-8-BOM") %>%
  mutate(endDate = as.Date(substr(endDate, 1, 10)),
         year = as.integer(format(endDate, "%Y"))) %>%
  inner_join(site_xwalk, by = c("siteID" = "neon_site"))

cat("\nLoaded", nrow(div_1m2), "1m2 plantSpecies rows and", nrow(div_nested),
    "nested presence rows, crosswalked to", length(unique(c(div_1m2$tower_id, div_nested$tower_id))),
    "tower_ids.\n")

# ============================================================================
# 4b. Vegetation Structure (vst_) data -- canopy-exposure species
#     classification, for the floristic_*_canopy metrics.
#
# Same discipline as the STEP 1 H5 investigation above: nothing below assumes
# a column name, category value, or ID scheme without confirming it against
# the real file first via find_one_column()/table(); anything that doesn't
# match documented expectation stop()s with the actual structure printed,
# rather than guessing.
# ============================================================================
find_one_column <- function(df, pattern, df_name) {
  candidates <- names(df)[str_detect(names(df), regex(pattern, ignore_case = TRUE))]
  if (length(candidates) == 0) {
    stop("No column matching /", pattern, "/ found in ", df_name, ". Columns present: ",
         paste(names(df), collapse = ", "),
         ". Structure differs from what this script expects -- investigate before proceeding.")
  }
  if (length(candidates) > 1) {
    stop("Multiple columns match /", pattern, "/ in ", df_name, ": ", paste(candidates, collapse = ", "),
         ". Ambiguous -- name the correct one explicitly rather than guessing.")
  }
  candidates[1]
}

# Like find_one_column(), but for exploratory checks where the column
# genuinely may not exist -- returns NA instead of stop()ing on zero or
# multiple matches, so a candidate-field search can report "not found"
# rather than aborting the script.
find_optional_column <- function(df, pattern) {
  candidates <- names(df)[str_detect(names(df), regex(pattern, ignore_case = TRUE))]
  if (length(candidates) != 1) return(NA_character_)
  candidates
}

vst_apparent <- read.csv(vst_apparent_path, fileEncoding = "UTF-8-BOM")

# ---- vst_mappingandtagging.csv: a real run found base read.csv() silently
# truncating this file (158,898 rows written at download time -> only 3,189
# rows loaded here), with "invalid input found on input connection" / "EOF
# within quoted string" -- the classic signature of an unescaped stray `"`
# inside one of this table's several free-text columns (remarks,
# identificationQualifier, morphospeciesIDRemarks, etc.), which makes
# read.csv()'s quote-balancing parser treat everything after that point as
# still inside one open quoted field, corrupting/dropping the remainder of
# the file. A line with an ODD count of `"` characters is the fingerprint of
# exactly this failure mode -- scanned for here (not hand-patched into the
# source file, per this project's write-safety convention of never editing a
# raw downloaded data file) so the actual offending line is visible rather
# than "somewhere in there".
mapping_raw_lines  <- readLines(vst_mapping_path, warn = FALSE)
mapping_quote_counts <- lengths(regmatches(mapping_raw_lines, gregexpr('"', mapping_raw_lines)))
mapping_odd_quote_lines <- which(mapping_quote_counts %% 2 == 1)
if (length(mapping_odd_quote_lines) > 0) {
  cat("\nvst_mappingandtagging.csv: line(s) with an ODD '\"' count (likely source of the",
      "CSV-quoting parse failure) at line number(s):\n")
  print(head(mapping_odd_quote_lines, 20))
} else {
  cat("\nvst_mappingandtagging.csv: no odd-quote-count line found by the raw scan -- the",
      "malformed-quote diagnosis may not be the exact mechanism here; check the",
      "readr::read_csv() problems() printout below instead.\n")
}

# Checked (not assumed) whether vst_apparentindividual.csv -- same download,
# same product -- shows the same fingerprint before leaving its read.csv()
# call unchanged.
apparent_raw_lines  <- readLines(vst_apparent_path, warn = FALSE)
apparent_quote_counts <- lengths(regmatches(apparent_raw_lines, gregexpr('"', apparent_raw_lines)))
apparent_odd_quote_lines <- which(apparent_quote_counts %% 2 == 1)
if (length(apparent_odd_quote_lines) > 0) {
  cat("vst_apparentindividual.csv ALSO has", length(apparent_odd_quote_lines),
      "odd-quote-count line(s) -- this file may have the same read.csv() truncation problem;",
      "investigate before trusting its row count (nrow =", nrow(vst_apparent), ") too.\n")
} else {
  cat("vst_apparentindividual.csv: no odd-quote-count lines found -- no evidence of the same",
      "parsing problem, read.csv() left unchanged for this file (nrow =", nrow(vst_apparent), ").\n")
}

# Switched to readr::read_csv() for THIS file only -- it tokenizes fields
# far more defensively around a malformed embedded quote (reports the
# problem via problems() and keeps parsing, rather than base read.csv()'s
# all-or-nothing quote-balancing that silently drops everything after the
# first bad quote). div_1m2Data.csv/div_10m2Data100m2Data.csv and
# vst_apparentindividual.csv are untouched -- no evidence found above that
# they share this problem.
vst_mapping <- read_csv(vst_mapping_path, show_col_types = FALSE, progress = FALSE)
mapping_parse_problems <- problems(vst_mapping)
if (nrow(mapping_parse_problems) > 0) {
  cat("\nvst_mappingandtagging.csv: read_csv() reported", nrow(mapping_parse_problems),
      "parsing problem(s):\n")
  print(mapping_parse_problems)
}
cat("\nLoaded", nrow(vst_mapping), "rows from vst_mappingandtagging.csv via read_csv()",
    "-- compare against the row count NEON_Download_VegStructure.R's own \"Wrote ...\"",
    "message reported at download time (158,898 on the run that surfaced this bug) to",
    "confirm the fix worked. If this is still far below that, the fix above did not fully",
    "resolve it -- investigate further rather than assuming success.\n")

cat("\n==== STEP 1B: vst_ (Vegetation Structure) structure investigation ====\n")
cat("vst_apparentindividual: ", nrow(vst_apparent), " rows, columns:\n", sep = "")
print(names(vst_apparent))
cat("\nSample rows:\n")
print(head(vst_apparent, 5))
cat("\nvst_mappingandtagging: ", nrow(vst_mapping), " rows, columns:\n", sep = "")
print(names(vst_mapping))
print(head(vst_mapping, 5))

individual_col_apparent <- find_one_column(vst_apparent, "^individualid$", "vst_apparentindividual")
site_col_apparent       <- find_one_column(vst_apparent, "^siteid$",       "vst_apparentindividual")
plot_col_apparent       <- find_one_column(vst_apparent, "^plotid$",       "vst_apparentindividual")
canopy_col              <- find_one_column(vst_apparent, "canopy.?position", "vst_apparentindividual")

individual_col_mapping  <- find_one_column(vst_mapping, "^individualid$", "vst_mappingandtagging")
taxon_col_mapping       <- find_one_column(vst_mapping, "^taxonid$",      "vst_mappingandtagging")
date_col_apparent       <- find_one_column(vst_apparent, "^date$", "vst_apparentindividual")
date_col_mapping        <- find_one_column(vst_mapping,  "^date$", "vst_mappingandtagging")

cat("\nCanopy-position column identified as '", canopy_col, "'. Unique values:\n", sep = "")
print(table(vst_apparent[[canopy_col]], useNA = "always"))

# ---- Duplicate individualID investigation, BEFORE joining. A real run
# found the join fanning out: 502,012 apparentindividual rows joined
# against 158,898 mapping rows produced 505,124 matched rows (MORE than the
# apparentindividual side alone), with dplyr's own "unexpected many-to-many
# relationship" warning. Investigated rather than assumed which side is
# actually non-unique.
#
# vst_apparentindividual: the SAME individualID appearing across multiple
# rows (different visit dates) is EXPECTED, not a defect -- every visit's
# canopyPosition is independent evidence and must be preserved, not
# deduplicated. The real anomaly to check for is a duplicate
# (individualID, date) PAIR -- the same individual measured twice on the
# same visit.
apparent_dupe_visits <- vst_apparent %>%
  count(.data[[individual_col_apparent]], .data[[date_col_apparent]]) %>%
  filter(n > 1)
cat("\nvst_apparentindividual:", n_distinct(vst_apparent[[individual_col_apparent]]),
    "distinct individualIDs across", nrow(vst_apparent), "rows (repeat-visit structure is",
    "expected and NOT deduplicated). Duplicate (individualID, date) pairs -- a genuine",
    "anomaly, distinct from expected repeat-visit structure:", nrow(apparent_dupe_visits), "\n")
if (nrow(apparent_dupe_visits) > 0) {
  cat("Sample duplicate (individualID, date) pairs in vst_apparentindividual:\n")
  print(head(apparent_dupe_visits, 10))
}

# vst_mappingandtagging: SHOULD be one static identity row per individualID
# (assigned once at tagging) -- any individualID appearing more than once
# here is the real candidate source of the join fan-out.
mapping_dupe_counts <- vst_mapping %>% count(.data[[individual_col_mapping]]) %>% filter(n > 1)
mapping_dupe_ids <- mapping_dupe_counts[[individual_col_mapping]]

cat("\nvst_mappingandtagging:", length(mapping_dupe_ids), "of",
    n_distinct(vst_mapping[[individual_col_mapping]]),
    "distinct individualIDs have MORE THAN ONE row (", nrow(vst_mapping), "total rows).\n")

conflict_ids <- character(0)
if (length(mapping_dupe_ids) > 0) {
  cat("Full rows for the first few duplicated individualIDs (inspect the actual pattern",
      "before deduplicating):\n")
  print(vst_mapping %>% filter(.data[[individual_col_mapping]] %in% head(mapping_dupe_ids, 5)) %>%
          arrange(.data[[individual_col_mapping]]))

  # ---- taxonID conflicts within duplicate individualID groups. A real run
  # found 172 such individualIDs; inspection showed duplicate rows sharing
  # consistent physical metadata (plot, stem distance/azimuth) but
  # differing taxonID across dates spanning YEARS (e.g. 2014 vs. 2025, 2022
  # vs. 2025). Read as a field identification corrected at a later revisit
  # -- an expected pattern in decade-spanning ecological monitoring, not
  # data corruption.
  #
  # DECISION: resolved by the SAME most-recent-date rule applied to every
  # duplicated individualID below (see dedup block) -- these 172 cases need
  # no special-casing, since "the latest record is authoritative" already
  # resolves both the general re-tagging duplication and the specific
  # taxonID disagreements identically. Before the superseded taxonID is
  # discarded for these specific individuals, an audit trail is written
  # (below, after the dedup) documenting exactly what was overridden --
  # this is NOT a silent resolution, it's a documented one.
  taxon_conflicts <- vst_mapping %>%
    filter(.data[[individual_col_mapping]] %in% mapping_dupe_ids) %>%
    group_by(.data[[individual_col_mapping]]) %>%
    summarise(n_taxa = n_distinct(.data[[taxon_col_mapping]]), .groups = "drop") %>%
    filter(n_taxa > 1)
  conflict_ids <- taxon_conflicts[[individual_col_mapping]]

  if (length(conflict_ids) > 0) {
    cat("\n", length(conflict_ids), " individualID(s) have duplicate rows disagreeing on ",
        "taxonID -- resolved by the most-recent-date rule below (see the audit trail written ",
        "to ./Data/NEON_FieldData/vst_taxonID_reassignments.csv for exactly what was ",
        "overridden).\n", sep = "")
  } else {
    cat("\nNo taxonID conflicts found within duplicate individualID groups -- duplicates appear",
        "to be benign re-tagging/re-mapping history with a consistent species call.\n")
  }
}

# ---- Deduplicate vst_mappingandtagging to one row per individualID,
# keeping the most recent `date` -- ONE consistent rule applied to ALL
# duplicated individualIDs, whether or not they show a taxonID conflict, so
# the 172 conflicting cases are not special-cased separately from the
# general re-tagging duplication. Applied BEFORE the join, not as a
# post-hoc filter on already-fanned-out rows. vst_apparentindividual is NOT
# deduplicated (its repeat-visit rows are legitimate, independent
# canopyPosition evidence, confirmed above).
vst_mapping_dated <- vst_mapping %>%
  mutate(.dedupe_date = suppressWarnings(as.Date(.data[[date_col_mapping]])))

# Exact ties on the max date within one individualID can't be broken by
# date alone -- reported here rather than silently, in case it ever
# actually occurs; the first such row encountered after arrange() is kept
# deterministically regardless.
exact_date_ties <- vst_mapping_dated %>%
  group_by(.data[[individual_col_mapping]]) %>%
  filter(.dedupe_date == max(.dedupe_date, na.rm = TRUE)) %>%
  summarise(n_at_max_date = n(), .groups = "drop") %>%
  filter(n_at_max_date > 1)
if (nrow(exact_date_ties) > 0) {
  cat("\n!!", nrow(exact_date_ties), "individualID(s) have an exact TIE on the maximum date --",
      "date alone cannot fully resolve these; the first such row encountered after arrange()",
      "is kept. Inspect if this matters for any of the taxonID-conflict cases specifically:\n")
  print(head(exact_date_ties, 10))
} else {
  cat("\nNo exact ties on the maximum date within any duplicated individualID -- the",
      "most-recent-date rule resolves every case unambiguously.\n")
}

# ---- Audit trail for taxonID reassignments -- written BEFORE the
# superseded rows are discarded by the dedup below. Documentation/QA only;
# not read anywhere downstream in the canopy-filtering logic itself.
if (length(conflict_ids) > 0) {
  conflict_rows <- vst_mapping_dated %>%
    filter(.data[[individual_col_mapping]] %in% conflict_ids) %>%
    arrange(.data[[individual_col_mapping]], desc(.dedupe_date))

  kept_rows <- conflict_rows %>%
    group_by(.data[[individual_col_mapping]]) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(individualID = .data[[individual_col_mapping]],
              new_taxonID  = .data[[taxon_col_mapping]],
              new_date     = .dedupe_date)

  superseded_rows <- conflict_rows %>%
    group_by(.data[[individual_col_mapping]]) %>%
    slice(-1) %>%
    ungroup() %>%
    transmute(individualID = .data[[individual_col_mapping]],
              old_taxonID  = .data[[taxon_col_mapping]],
              old_date     = .dedupe_date)

  reassignment_audit <- superseded_rows %>%
    inner_join(kept_rows, by = "individualID") %>%
    filter(old_taxonID != new_taxonID) %>%
    select(individualID, old_taxonID, old_date, new_taxonID, new_date)

  audit_path <- "./Data/NEON_FieldData/vst_taxonID_reassignments.csv"
  write.csv(reassignment_audit, audit_path, row.names = FALSE)
  cat("\nWrote", nrow(reassignment_audit), "taxonID reassignment(s) for", length(conflict_ids),
      "individualID(s) to", audit_path, "(a conflicting individualID with more than 2 total",
      "rows produces one audit row per superseded record that actually disagrees with the",
      "kept taxonID).\n")
}

vst_mapping_deduped <- vst_mapping_dated %>%
  arrange(.data[[individual_col_mapping]], desc(.dedupe_date)) %>%
  distinct(across(all_of(individual_col_mapping)), .keep_all = TRUE) %>%
  select(-.dedupe_date)

cat("\nDeduplicated vst_mappingandtagging:", nrow(vst_mapping), "->", nrow(vst_mapping_deduped),
    "rows (one per individualID;", length(mapping_dupe_ids), "individualID(s) were",
    "deduplicated, including", length(conflict_ids), "with a taxonID reassignment).\n")

# ---- The 38,758 duplicate (individualID, date) pairs found in
# vst_apparentindividual (see apparent_dupe_visits above) are NOT
# deduplicated -- repeat-visit structure is preserved. But same-visit
# duplicates deserve the same identity-conflict check just applied to
# vst_mappingandtagging: canopyPosition is EXPECTED to vary/disagree
# between duplicate rows (independent readings, and the any-exposed-
# individual rule is designed around exactly that), but a disagreement on
# an identity-relevant field (growthForm) at the SAME visit would mean
# individualID isn't uniquely identifying at that resolution -- a genuine
# data-quality concern, not something the any()-based classification can
# paper over. Flagged the same way as the taxonID conflict (stop(), not a
# silent resolution) if found.
growthform_col_apparent <- find_optional_column(vst_apparent, "^growthform$")

if (nrow(apparent_dupe_visits) == 0) {
  cat("\nvst_apparentindividual: no duplicate (individualID, date) pairs exist, so there is",
      "nothing to identity-check here.\n")
} else if (is.na(growthform_col_apparent)) {
  cat("\nvst_apparentindividual: no growthForm-like column found to cross-check duplicate",
      "(individualID, date) groups against -- could not confirm identity consistency beyond",
      "canopyPosition itself. canopyPosition disagreement within a duplicate group is treated",
      "as expected (independent readings), not an identity conflict, consistent with why these",
      "rows are preserved rather than deduplicated.\n")
} else {
  apparent_identity_conflicts <- vst_apparent %>%
    inner_join(apparent_dupe_visits %>% select(all_of(c(individual_col_apparent, date_col_apparent))),
               by = c(individual_col_apparent, date_col_apparent)) %>%
    group_by(.data[[individual_col_apparent]], .data[[date_col_apparent]]) %>%
    summarise(n_growthforms = n_distinct(.data[[growthform_col_apparent]]), .groups = "drop") %>%
    filter(n_growthforms > 1)

  if (nrow(apparent_identity_conflicts) > 0) {
    cat("\n!! ", nrow(apparent_identity_conflicts), " duplicate (individualID, date) group(s) in ",
        "vst_apparentindividual disagree on '", growthform_col_apparent, "':\n", sep = "")
    print(head(apparent_identity_conflicts, 10))
    stop(nrow(apparent_identity_conflicts), " duplicate (individualID, date) group(s) in ",
         "vst_apparentindividual disagree on '", growthform_col_apparent, "' -- a genuine ",
         "identity conflict at the SAME visit, not independent multi-year evidence like the ",
         "vst_mappingandtagging case. Investigate before trusting the any-exposed-individual ",
         "classification for these individuals; not resolved automatically.")
  }
  cat("\nvst_apparentindividual: no duplicate (individualID, date) group disagrees on '",
      growthform_col_apparent, "' -- the any-exposed-individual canopy classification rule is",
      " confirmed sufficient without further changes here. Same-visit duplicate rows describe",
      " the same individual consistently; multiple canopyPosition votes (agreeing or",
      " disagreeing on EXPOSURE, which is expected to genuinely vary by reading) don't reflect",
      " an identity problem the way a taxonID mismatch would.\n", sep = "")
}

# ---- vst_apparentindividual (repeated per-visit measurements, incl.
# canopyPosition) carries no taxonID of its own -- species identity lives on
# vst_mappingandtagging (assigned once per individual at tagging). Joining by
# individualID is the documented NEON DP1.10098.001 linkage; if that
# assumption is wrong for the real files, the join below will surface it as
# zero matched rows rather than silently producing an empty/garbage lookup.
# relationship = "many-to-one" makes the expected cardinality explicit --
# vst_apparentindividual is legitimately many rows per individualID (visit
# history), vst_mappingandtagging is now exactly one per individualID after
# deduplication -- so any future regression in that uniqueness (e.g. a
# re-download reintroducing duplicates) fails loudly here instead of
# silently fanning out again.
vst_joined <- vst_apparent %>%
  transmute(individualID = .data[[individual_col_apparent]],
            siteID       = .data[[site_col_apparent]],
            plotID       = .data[[plot_col_apparent]],
            canopyPosition = .data[[canopy_col]]) %>%
  inner_join(
    vst_mapping_deduped %>% transmute(individualID = .data[[individual_col_mapping]],
                                       taxonID      = .data[[taxon_col_mapping]]),
    by = "individualID",
    relationship = "many-to-one"
  )

if (nrow(vst_joined) == 0) {
  stop("Joining vst_apparentindividual to vst_mappingandtagging by individualID produced ",
       "zero rows. This assumes the two tables share an individualID key per NEON ",
       "DP1.10098.001's documented schema -- that assumption doesn't hold for these files. ",
       "Investigate the real join key before proceeding; do not guess.")
}

# Belt-and-suspenders check alongside dplyr's own relationship = "many-to-one"
# guard: the join must never produce more rows than went in on the "many"
# side once the "one" side is genuinely unique.
if (nrow(vst_joined) > nrow(vst_apparent)) {
  stop("Join still fans out after deduplication: ", nrow(vst_joined), " matched rows exceeds ",
       "the ", nrow(vst_apparent), " vst_apparentindividual rows going in. relationship = ",
       "\"many-to-one\" should have caught this as an error -- investigate further before ",
       "trusting downstream canopy classification.")
}

join_match_rate <- nrow(vst_joined) / nrow(vst_apparent)
cat("\nJoined vst_apparentindividual x vst_mappingandtagging by individualID:",
    nrow(vst_joined), "of", nrow(vst_apparent), "apparentindividual rows matched",
    "(", round(100 * join_match_rate, 1), "%).\n")
if (join_match_rate < 0.5) {
  cat("!! Match rate is still under 50% -- a run that previously saw 14,734/502,012 (2.9%)",
      "matched against a truncated vst_mappingandtagging.csv (see the read_csv() fix above)",
      "should climb sharply once the parsing fix is in effect. If it's STILL this low after",
      "that fix, that points to a separate, real issue (e.g. individualIDs genuinely absent",
      "from the mapping table) -- investigate rather than assuming the parsing fix alone",
      "explains the remainder.\n")
}

# ---- Attempt to confirm exposed/shaded placement against NEON's
# authoritative categoricalCodes_10098.csv definitions, rather than the
# label text alone. A real run confirmed this file's actual column
# structure is name/pubCode/description/startDate/endDate -- NOT the
# fieldName/name/definition layout originally guessed here, which meant
# the lookup below could never have matched anything on the first attempt.
# Corrected to filter on `name` directly (the category label itself, e.g.
# "Mostly shaded") against the known canopyPosition vocabulary, since this
# structure has no fieldName column to scope the search to one field first.
# categorical_codes/cc_name_col/cc_desc_col are kept in scope for reuse by
# the vst_non-woody section (4c) below, so its own exposure-proxy field (if
# any) can be looked up the same way. Still best-effort: a missing file or
# unexpected layout falls back to the documented judgment call below rather
# than blocking the whole script on a reference-table lookup.
categorical_codes_path <- "./Data/NEON_FieldData/categoricalCodes_10098.csv"
canopy_categories_of_interest <- c("Full shade", "Mostly shaded", "Partially shaded",
                                    "Open grown", "Full sun")
categorical_codes <- NULL
cc_name_col <- NA_character_
cc_desc_col <- NA_character_

if (file.exists(categorical_codes_path)) {
  categorical_codes <- read.csv(categorical_codes_path, fileEncoding = "UTF-8-BOM")
  cc_name_col <- find_optional_column(categorical_codes, "^name$")
  cc_desc_col <- find_optional_column(categorical_codes, "^description$")
  cc_pub_col  <- find_optional_column(categorical_codes, "^pubcode$")

  if (!is.na(cc_name_col) && !is.na(cc_desc_col)) {
    canopy_position_definitions <- categorical_codes %>%
      filter(.data[[cc_name_col]] %in% canopy_categories_of_interest) %>%
      transmute(category = .data[[cc_name_col]],
                pubCode = if (!is.na(cc_pub_col)) .data[[cc_pub_col]] else NA,
                definition = .data[[cc_desc_col]])
    cat("\n==== categoricalCodes_10098.csv: canopyPosition definitions",
        "(matched by category label against the known vocabulary) ====\n")
    if (nrow(canopy_position_definitions) > 0) {
      print(as.data.frame(canopy_position_definitions))
      cat("REVIEW the definitions printed above against the exposed_categories/",
          "shaded_categories placement below -- especially \"Mostly shaded\", added on",
          "an ordinal judgment call (see comment) because this codebase could not read",
          "this file's real contents at the time that call was written. Correct the two",
          "vectors below if the printed definitions disagree with any of the five",
          "placements.\n")
    } else {
      cat("(no rows matched any known canopyPosition category label -- this file may not",
          "cover this field, or use different label text. Falling back to the ordinal",
          "judgment call below.)\n")
    }
  } else {
    cat("\ncategoricalCodes_10098.csv exists but doesn't have the expected name/description",
        "columns (found:", paste(names(categorical_codes), collapse = ", "),
        ") -- cannot confirm category definitions against it programmatically here;",
        "inspect it by hand if the judgment call below needs revisiting.\n")
  }
} else {
  cat("\ncategoricalCodes_10098.csv not found at", categorical_codes_path, "-- cannot confirm",
      "canopyPosition category definitions against an authoritative source in this run;",
      "using the ordinal judgment call below instead.\n")
}

# ---- exposed vs. shaded category mapping.
# "Full sun"/"Open grown"/"Partially shaded" (exposed) and "Full shade"
# (shaded) were the original PROVISIONAL placements, based on NEON's
# documented vocabulary, not a confirmed definition-file read. A real run
# then surfaced a fifth category, "Mostly shaded" (39,490 of 502,012 rows),
# which correctly stop()ed the original version of this script rather than
# silently guessing.
#
# "Mostly shaded" is placed on the SHADED side here by ORDINAL JUDGMENT
# CALL, NOT a confirmed categoricalCodes_10098.csv definition -- this
# codebase has no access to real NEON data/files to read that definition
# text directly (see the investigation block immediately above, which DOES
# read and print the real definitions when this script is actually run;
# check that printout and correct this placement if it disagrees). The
# judgment: NEON's canopyPosition vocabulary orders by degree of light
# exposure, roughly Full shade < Mostly shaded < Partially shaded <
# Partially shaded < Full sun / Open grown; "Partially shaded" (majority
# sun -- already placed exposed) and "Mostly shaded" (majority shade) sit
# on opposite sides of the 50% line implied by their own labels, so
# "Mostly shaded" goes with "Full shade" rather than with "Partially
# shaded". All five categories' placements, including the four assumed
# before this fix, should be re-verified against the printed definitions
# above rather than trusted as previously confirmed.
exposed_categories <- c("Full sun", "Open grown", "Partially shaded")
shaded_categories  <- c("Full shade", "Mostly shaded")

observed_categories <- unique(na.omit(vst_joined$canopyPosition))
unmapped_categories <- setdiff(observed_categories, c(exposed_categories, shaded_categories))
if (length(unmapped_categories) > 0) {
  stop("Unmapped canopyPosition categories found in real data: ",
       paste(unmapped_categories, collapse = ", "),
       ". Extend exposed_categories/shaded_categories above with an explicit decision ",
       "about which side of the exposed/shaded line each belongs on -- do not guess.")
}

# ============================================================================
# 4c. vst_non-woody -- investigate for a usable canopy/light-exposure proxy
# to extend the *_canopy metrics beyond woody species. vst_apparentindividual
# only covers trees/shrubs, so without this the canopy filter can only ever
# exclude confirmed-shaded TREES, never understory herbs/forbs/graminoids --
# a real motivating case for this whole addition. Nothing below assumes this
# table's structure, columns, or category vocabulary; everything is
# investigated live and printed before any integration decision is made. If
# no usable field is found, this section reports that clearly and leaves
# vst_joined (and therefore the classification below) untouched -- it does
# NOT fabricate an exposure signal from an unrelated field (cover, height,
# growth form) just to close the coverage gap.
# ============================================================================
vst_nonwoody <- read.csv(vst_nonwoody_path, fileEncoding = "UTF-8-BOM")

cat("\n==== STEP 1C: vst_non-woody structure investigation ====\n")
cat("vst_non-woody: ", nrow(vst_nonwoody), " rows, columns:\n", sep = "")
print(names(vst_nonwoody))
cat("\nSample rows:\n")
print(head(vst_nonwoody, 5))

nonwoody_taxon_col  <- find_optional_column(vst_nonwoody, "^taxonid$")
nonwoody_site_col   <- find_optional_column(vst_nonwoody, "^siteid$")
nonwoody_plot_col   <- find_optional_column(vst_nonwoody, "^plotid$")
nonwoody_canopy_col <- find_optional_column(vst_nonwoody, "canopy.?position")
nonwoody_other_candidates <- setdiff(
  names(vst_nonwoody)[str_detect(names(vst_nonwoody),
    regex("cover|height|growthform|light|exposure|shade|sun", ignore_case = TRUE))],
  na.omit(nonwoody_canopy_col)
)

cat("\nvst_non-woody carries taxonID directly: ", !is.na(nonwoody_taxon_col),
    "; siteID directly: ", !is.na(nonwoody_site_col),
    "; plotID directly: ", !is.na(nonwoody_plot_col), ".\n", sep = "")
cat(if (!is.na(nonwoody_taxon_col)) {
  "No separate identity-mapping join needed -- taxonID is carried directly per row (non-woody\nsurveys apparently don't tag/map individual plants the way vst_apparentindividual does).\n"
} else {
  "No direct taxonID column found -- this table would need an identity-mapping join (mirroring\nvst_apparentindividual x vst_mappingandtagging) before species could be classified from it;\nnot attempted since no direct taxonID path was found.\n"
})

cat("canopyPosition-analogous field found in vst_non-woody: ",
    ifelse(is.na(nonwoody_canopy_col), "NONE", nonwoody_canopy_col), "\n", sep = "")
if (length(nonwoody_other_candidates) > 0) {
  cat("Other candidate exposure-proxy-ish column(s) present (NOT integrated automatically --",
      "converting cover/height/growth-form into exposed/shaded evidence would require a new,",
      "undiscussed thresholding rule this script does not fabricate):",
      paste(nonwoody_other_candidates, collapse = ", "), "\n")
}

nonwoody_usable <- !is.na(nonwoody_canopy_col) && !is.na(nonwoody_taxon_col) &&
                   !is.na(nonwoody_site_col) && !is.na(nonwoody_plot_col)

if (nonwoody_usable) {
  cat("\nvst_non-woody has a usable canopyPosition-analogous field ('", nonwoody_canopy_col,
      "') with direct taxonID/siteID/plotID -- integrating into the unified canopy-evidence",
      " table below. Unique values:\n", sep = "")
  print(table(vst_nonwoody[[nonwoody_canopy_col]], useNA = "always"))

  # Re-check categoricalCodes_10098.csv for THIS field's real categories,
  # reusing the categorical_codes/cc_name_col/cc_desc_col objects already
  # loaded above (same confirmed name/description structure).
  nonwoody_observed_categories <- unique(na.omit(vst_nonwoody[[nonwoody_canopy_col]]))
  if (!is.null(categorical_codes) && !is.na(cc_name_col) && !is.na(cc_desc_col)) {
    nonwoody_definitions <- categorical_codes %>%
      filter(.data[[cc_name_col]] %in% nonwoody_observed_categories) %>%
      transmute(category = .data[[cc_name_col]], definition = .data[[cc_desc_col]])
    cat("\n==== categoricalCodes_10098.csv: vst_non-woody '", nonwoody_canopy_col,
        "' definitions ====\n", sep = "")
    if (nrow(nonwoody_definitions) > 0) {
      print(as.data.frame(nonwoody_definitions))
      cat("CONFIRMED against categoricalCodes_10098.csv -- not a judgment call, unlike the",
          "original canopyPosition placements above.\n")
    } else {
      cat("(no matching definitions found -- category placement below, if any, is a judgment",
          "call, not a confirmed definition.)\n")
    }
  } else {
    cat("\ncategoricalCodes_10098.csv unavailable/unusable -- category placement below, if any,",
        "is a judgment call, not a confirmed definition.\n")
  }

  # This field reuses the SAME exposed_categories/shaded_categories vocabulary
  # already established for woody canopyPosition -- only valid because NEON's
  # canopyPosition vocabulary is shared vocabulary, not because non-woody
  # categories are assumed to match without checking. Any category found here
  # that isn't already in one of those two vectors stops rather than guesses.
  nonwoody_unmapped <- setdiff(nonwoody_observed_categories, c(exposed_categories, shaded_categories))
  if (length(nonwoody_unmapped) > 0) {
    stop("vst_non-woody's '", nonwoody_canopy_col, "' field contains categories not covered by ",
         "exposed_categories/shaded_categories: ", paste(nonwoody_unmapped, collapse = ", "),
         ". Extend those vectors (informed by the categoricalCodes_10098.csv definitions ",
         "printed above, if found) before proceeding -- do not guess.")
  }

  nonwoody_evidence <- vst_nonwoody %>%
    transmute(individualID   = NA_character_,
              siteID         = .data[[nonwoody_site_col]],
              plotID         = .data[[nonwoody_plot_col]],
              taxonID        = .data[[nonwoody_taxon_col]],
              canopyPosition = .data[[nonwoody_canopy_col]])

  vst_joined <- vst_joined %>% mutate(source = "vst_apparentindividual") %>%
    bind_rows(nonwoody_evidence %>% mutate(source = "vst_non-woody"))

  cat("\nCombined canopy-evidence table:", sum(vst_joined$source == "vst_apparentindividual"),
      "woody rows +", sum(vst_joined$source == "vst_non-woody"), "non-woody rows =",
      nrow(vst_joined), "total.\n")

  # ---- Flag species with evidence from BOTH sources -- not silently
  # resolved. build_site_canopy_lookup() below already pools ALL rows for a
  # taxonID (regardless of source) before voting exposed/understory/
  # unmeasured, so "exposed wins if either source shows it" falls out of
  # the existing any-exposed-individual rule with no extra code needed;
  # this just reports how often it actually matters.
  both_source_taxa <- vst_joined %>%
    filter(!is.na(taxonID)) %>%
    group_by(taxonID) %>%
    summarise(n_sources = n_distinct(source), .groups = "drop") %>%
    filter(n_sources > 1)
  cat(nrow(both_source_taxa), "taxonID(s) have canopy evidence from BOTH vst_apparentindividual",
      "and vst_non-woody -- resolved by the existing any-exposed-individual rule (exposed wins",
      "if either source shows it for that species), not a separate cross-source resolution",
      "step.\n")
} else {
  cat("\nvst_non-woody does NOT have a usable canopyPosition-analogous field together with the",
      "taxonID/siteID/plotID needed to integrate it automatically -- the herb/forb/graminoid",
      "canopy-filter coverage gap REMAINS OPEN. Not fabricating an exposure signal from an",
      "unrelated field (cover/height/growth form) per explicit instruction; vst_joined is left",
      "unchanged (woody evidence only). A different approach would be needed to close this gap",
      "for non-woody species.\n")
}

# ---- STEP 3: plot-level vs. site-level linkage -- investigate the actual
# plotID overlap rather than assuming plot-level granularity is available.
vst_plot_ids <- unique(vst_joined$plotID)
div_plot_ids <- unique(c(div_1m2$plotID, div_nested$plotID))
plot_overlap <- intersect(vst_plot_ids, div_plot_ids)
plot_overlap_frac <- if (length(vst_plot_ids) == 0) 0 else length(plot_overlap) / length(vst_plot_ids)

cat("\nvst_ plotID / diversity-data plotID overlap:", length(plot_overlap), "of",
    length(vst_plot_ids), "vst_ plots (", round(100 * plot_overlap_frac, 1), "%).\n")

# Threshold is a flagged judgment call, not a NEON-documented rule: >50% of
# vst_ plots recognized as diversity plots is treated as "linkage usable",
# below that the schemes are assumed unrelated enough to fall back to
# site-level pooling.
canopy_linkage_granularity <- if (plot_overlap_frac > 0.5) "plot" else "site"
cat("Canopy linkage granularity achieved: '", canopy_linkage_granularity, "'.\n", sep = "")

# ---- NA canopyPosition diagnostic. A real run found 338,712 of 502,012
# vst_apparentindividual rows (67%) have canopyPosition == NA (apparently
# not recorded at every visit for every individual). FIXED below: NA rows
# are now excluded from exposed/shaded evidence entirely -- a data gap is
# not confirmed non-exposure, matching the unmeasured_species_treatment
# philosophy already established for species entirely absent from vst_.
# The count printed here is exactly the set of taxon(+plot) groups this fix
# moves from the old (incorrect) "understory" classification to the
# corrected "unmeasured" one -- every one of these groups has zero
# individuals with a usable (non-NA) canopyPosition, so build_site_
# canopy_lookup() below now omits them from its lookup table entirely,
# which is what makes classify_canopy_status() fall through to its
# "unmeasured" default for them.
na_canopy_group_cols <- if (canopy_linkage_granularity == "plot") c("plotID", "taxonID") else "taxonID"
na_canopy_diagnostic <- vst_joined %>%
  group_by(across(all_of(na_canopy_group_cols))) %>%
  summarise(all_na = all(is.na(canopyPosition)), .groups = "drop")
cat("\nNA canopyPosition diagnostic:", sum(is.na(vst_joined$canopyPosition)), "of",
    nrow(vst_joined), "vst_joined rows have NA canopyPosition;", sum(na_canopy_diagnostic$all_na),
    "of", nrow(na_canopy_diagnostic), canopy_linkage_granularity, "-level taxon group(s) have",
    "canopyPosition == NA for EVERY recorded individual -- these move from 'understory' (the",
    "prior, incorrect behavior) to 'unmeasured' (the corrected behavior) as a result of the",
    "NA-handling fix in build_site_canopy_lookup() below.\n")

# ---- STEP 2: per-species canopy-status classification, "any exposed
# individual" rule (permissive -- flagged choice, see task note). A stricter
# alternative (e.g. >50% of a species' measured individuals exposed, where N
# is large enough to be meaningful) is a legitimate alternative not
# implemented here. Computed ONCE per site, reused across every
# bout/plot_scope/temporal_scope combo for that site (it doesn't vary by any
# of those axes).
#
# NA canopyPosition individuals are excluded from exposed/shaded evidence
# entirely (neither "this species is exposed" nor "this species is
# understory" evidence) -- a species is "exposed" if ANY individual with a
# non-NA canopyPosition is in an exposed category; "understory" if it has
# at least one non-NA individual and none of them are exposed; and if EVERY
# individual has NA canopyPosition (no usable classification at all), the
# group is dropped from the lookup table here so classify_canopy_status()'s
# existing missing-key fallback classifies it "unmeasured" -- the same
# bucket used for species entirely absent from vst_ data, since both cases
# are the same thing: no measured evidence either way, not confirmed
# non-exposure.
build_site_canopy_lookup <- function(nsite) {
  va_site <- vst_joined %>% filter(siteID == nsite)
  if (nrow(va_site) == 0) return(tibble(key = character(), canopy_status = character()))

  group_cols <- if (canopy_linkage_granularity == "plot") c("plotID", "taxonID") else "taxonID"
  va_site %>%
    mutate(is_exposed = canopyPosition %in% exposed_categories) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      canopy_status = case_when(
        all(is.na(canopyPosition)) ~ NA_character_,
        any(is_exposed)            ~ "exposed",
        TRUE                       ~ "understory"
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(canopy_status)) %>%
    mutate(key = if (canopy_linkage_granularity == "plot") paste(plotID, taxonID, sep = "\r") else taxonID) %>%
    select(key, canopy_status)
}

all_neon_sites <- unique(site_xwalk$neon_site)
canopy_lookup_by_site <- set_names(map(all_neon_sites, build_site_canopy_lookup), all_neon_sites)

exposed_statuses <- if (unmeasured_species_treatment == "include") c("exposed", "unmeasured") else "exposed"
cat("\nUnmeasured-species treatment: '", unmeasured_species_treatment,
    "' -- unmeasured species are ", if (unmeasured_species_treatment == "include") "INCLUDED in" else "EXCLUDED from",
    " the canopy-filtered metrics.\n", sep = "")

# taxonIDs/plotIDs must be same-length vectors from one siteID; returns
# "exposed" / "understory" / "unmeasured" per element.
classify_canopy_status <- function(taxonIDs, plotIDs, nsite) {
  lookup <- canopy_lookup_by_site[[nsite]]
  if (is.null(lookup) || nrow(lookup) == 0) return(rep("unmeasured", length(taxonIDs)))

  key <- if (canopy_linkage_granularity == "plot") paste(plotIDs, taxonIDs, sep = "\r") else taxonIDs
  status <- lookup$canopy_status[match(key, lookup$key)]
  status[is.na(status)] <- "unmeasured"
  status
}

# ============================================================================
# 5. Per-(tower_id, year, plot_scope, bout) metric computation
# ============================================================================
compute_combo_metrics <- function(tid, yr, plot_scope, bout) {
  d1 <- div_1m2   %>% filter(tower_id == tid, year == yr, boutNumber == bout)
  dn <- div_nested %>% filter(tower_id == tid, year == yr, boutNumber == bout)
  if (plot_scope == "tower") {
    d1 <- d1 %>% filter(plotType == "tower")
    dn <- dn %>% filter(plotType == "tower")
  }

  # ---- metric 1: abundance-weighted Hill-Shannon (q=1) from 1m2 % cover ----
  # mean cover per species per plot, across the SIX 1m2 subplots (subplots
  # where a species wasn't recorded count as 0 cover, hence /6, not /n rows)
  plot_cover <- d1 %>%
    group_by(plotID, taxonID) %>%
    summarise(mean_cover = sum(percentCover, na.rm = TRUE) / 6, .groups = "drop")

  n_plots <- length(unique(plot_cover$plotID))

  if (n_plots == 0) {
    floristic_shannon_mean  <- NA_real_
    floristic_shannon_gamma <- NA_real_
  } else {
    comm <- plot_cover %>%
      pivot_wider(names_from = taxonID, values_from = mean_cover, values_fill = 0) %>%
      column_to_rownames("plotID") %>%
      as.matrix()

    floristic_shannon_mean  <- mean(hill_taxa(comm, q = 1))
    floristic_shannon_gamma <- hill_taxa_parti(comm, q = 1, show_warning = FALSE)$TD_gamma
  }

  # ---- metric 2: gamma richness from combined presence data, per plot too
  # (per-plot breakdown is kept only for the sanity-floor check) ----
  per_plot_taxa <- bind_rows(
    select(d1, plotID, taxonID),
    select(dn, plotID, taxonID)
  ) %>% filter(!is.na(taxonID), nzchar(taxonID))

  per_plot_additional <- dn %>%
    filter(!is.na(additionalSpecies), nzchar(additionalSpecies)) %>%
    select(plotID, additionalSpecies) %>%
    mutate(taxonID = str_split(additionalSpecies, ";")) %>%
    unnest(taxonID) %>%
    mutate(taxonID = str_trim(taxonID)) %>%
    filter(nzchar(taxonID)) %>%
    select(plotID, taxonID)

  per_plot_all <- bind_rows(per_plot_taxa, per_plot_additional) %>% distinct(plotID, taxonID)

  floristic_richness <- length(unique(per_plot_all$taxonID))
  max_single_plot_richness <- if (nrow(per_plot_all) == 0) 0L else {
    per_plot_all %>% count(plotID) %>% pull(n) %>% max()
  }

  # ---- canopy-filtered metrics: same unions/matrix as above, restricted to
  # taxa classified canopy-exposed via vst_ (Section 4b). nsite_for_canopy
  # matches the siteID key the Section 4b lookup was built on (NEON site
  # code, pre-crosswalk -- retained on d1/dn from the Section 4 join). ----
  nsite_for_canopy <- if (nrow(d1) > 0) d1$siteID[1] else if (nrow(dn) > 0) dn$siteID[1] else NA_character_

  if (nrow(per_plot_all) == 0 || is.na(nsite_for_canopy)) {
    per_plot_all_canopy <- per_plot_all[0, ]
  } else {
    per_plot_all_canopy <- per_plot_all %>%
      mutate(canopy_status = classify_canopy_status(taxonID, plotID, nsite_for_canopy)) %>%
      filter(canopy_status %in% exposed_statuses) %>%
      select(plotID, taxonID)
  }

  floristic_richness_canopy <- length(unique(per_plot_all_canopy$taxonID))
  max_single_plot_richness_canopy <- if (nrow(per_plot_all_canopy) == 0) 0L else {
    per_plot_all_canopy %>% count(plotID) %>% pull(n) %>% max()
  }

  if (n_plots == 0 || is.na(nsite_for_canopy)) {
    floristic_shannon_mean_canopy  <- NA_real_
    floristic_shannon_gamma_canopy <- NA_real_
  } else {
    canopy_taxa <- unique(per_plot_all_canopy$taxonID)
    canopy_cols <- intersect(colnames(comm), canopy_taxa)
    if (length(canopy_cols) == 0) {
      floristic_shannon_mean_canopy  <- NA_real_
      floristic_shannon_gamma_canopy <- NA_real_
    } else {
      comm_canopy <- comm[, canopy_cols, drop = FALSE]
      floristic_shannon_mean_canopy  <- mean(hill_taxa(comm_canopy, q = 1))
      floristic_shannon_gamma_canopy <- hill_taxa_parti(comm_canopy, q = 1, show_warning = FALSE)$TD_gamma
    }
  }

  # ---- bout end date, for peak_flight matching ----
  end_dates <- if (nrow(d1) > 0) d1$endDate else dn$endDate
  bout_end_date <- if (length(end_dates) == 0) as.Date(NA) else median(end_dates, na.rm = TRUE)

  list(n_plots = n_plots,
       floristic_shannon_mean = floristic_shannon_mean,
       floristic_shannon_gamma = floristic_shannon_gamma,
       floristic_richness = floristic_richness,
       max_single_plot_richness = max_single_plot_richness,
       floristic_shannon_mean_canopy = floristic_shannon_mean_canopy,
       floristic_shannon_gamma_canopy = floristic_shannon_gamma_canopy,
       floristic_richness_canopy = floristic_richness_canopy,
       max_single_plot_richness_canopy = max_single_plot_richness_canopy,
       bout_end_date = bout_end_date,
       d1 = d1, dn = dn)
}

# ============================================================================
# 6. Build the full site x year x plot_scope grid and compute every row
# ============================================================================
plot_scopes <- c("tower", "all")
all_years   <- sort(unique(c(div_1m2$year, div_nested$year)))

grid <- expand_grid(site_xwalk %>% distinct(tower_id, neon_site), year = all_years, plot_scope = plot_scopes)

cat("\n==== Computing field diversity for", nrow(grid), "tower-year-plot_scope combos ====\n")

result_rows <- list()

for (i in seq_len(nrow(grid))) {
  tid <- grid$tower_id[i]; nsite <- grid$neon_site[i]
  yr  <- grid$year[i];     ps    <- grid$plot_scope[i]

  d1_scope <- div_1m2 %>% filter(tower_id == tid, year == yr)
  dn_scope <- div_nested %>% filter(tower_id == tid, year == yr)
  if (ps == "tower") {
    d1_scope <- d1_scope %>% filter(plotType == "tower")
    dn_scope <- dn_scope %>% filter(plotType == "tower")
  }

  bouts <- sort(unique(c(d1_scope$boutNumber, dn_scope$boutNumber)))

  if (length(bouts) == 0) {
    result_rows[[length(result_rows) + 1]] <- tibble(
      tower_id = tid, neon_site = nsite, year = yr, plot_scope = ps,
      temporal_scope = NA_character_, bout = NA_integer_, n_plots = NA_integer_,
      flight_date_matched = as.Date(NA), bout_end_date = as.Date(NA),
      floristic_shannon_mean = NA_real_, floristic_shannon_gamma = NA_real_,
      floristic_richness = NA_integer_,
      floristic_shannon_mean_canopy = NA_real_, floristic_shannon_gamma_canopy = NA_real_,
      floristic_richness_canopy = NA_integer_, status = "no data")
    next
  }

  bout_rows <- map(bouts, function(b) {
    m <- compute_combo_metrics(tid, yr, ps, b)
    tibble(tower_id = tid, neon_site = nsite, year = yr, plot_scope = ps,
           temporal_scope = "per_bout", bout = b, n_plots = m$n_plots,
           flight_date_matched = as.Date(NA), bout_end_date = m$bout_end_date,
           floristic_shannon_mean = m$floristic_shannon_mean,
           floristic_shannon_gamma = m$floristic_shannon_gamma,
           floristic_richness = m$floristic_richness,
           floristic_shannon_mean_canopy = m$floristic_shannon_mean_canopy,
           floristic_shannon_gamma_canopy = m$floristic_shannon_gamma_canopy,
           floristic_richness_canopy = m$floristic_richness_canopy, status = "ok")
  }) %>% bind_rows()

  result_rows[[length(result_rows) + 1]] <- bout_rows

  # ---- peak_flight: nearest bout to the representative flight date, within
  # a 30-day tolerance; absent entirely (not NA-filled) if no flight date or
  # no bout qualifies. ----
  fd_row <- flight_dates %>% filter(tower_id == tid, year == yr)
  if (nrow(fd_row) == 1 && !is.na(fd_row$flight_date[1])) {
    fdate <- fd_row$flight_date[1]
    day_diff <- abs(as.numeric(bout_rows$bout_end_date - fdate))
    if (!all(is.na(day_diff))) {
      min_diff <- min(day_diff, na.rm = TRUE)
      if (min_diff <= flight_match_tolerance_days) {
        candidate_idx <- which(day_diff == min_diff)
        best_idx <- candidate_idx[which.min(bout_rows$bout_end_date[candidate_idx])]
        peak_row <- bout_rows[best_idx, ]
        peak_row$temporal_scope <- "peak_flight"
        peak_row$flight_date_matched <- fdate
        result_rows[[length(result_rows) + 1]] <- peak_row
      }
    }
  }
}

field_diversity_long <- bind_rows(result_rows) %>%
  select(tower_id, neon_site, year, plot_scope, temporal_scope, bout, n_plots,
         flight_date_matched, bout_end_date, floristic_shannon_mean,
         floristic_shannon_gamma, floristic_richness,
         floristic_shannon_mean_canopy, floristic_shannon_gamma_canopy,
         floristic_richness_canopy, status)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(field_diversity_long, out_csv, row.names = FALSE)
cat("\nWrote", nrow(field_diversity_long), "rows to", out_csv, "\n")

# ============================================================================
# 7. VALIDATE
# ============================================================================
cat("\n==== Peak-flight match inventory ====\n")
site_years_with_data  <- grid %>% distinct(tower_id, year)
site_years_with_flight <- flight_dates %>% filter(!is.na(flight_date)) %>% distinct(tower_id, year)
site_years_matched     <- field_diversity_long %>% filter(temporal_scope == "peak_flight") %>%
  distinct(tower_id, year)

for (i in seq_len(nrow(site_years_with_data))) {
  tid <- site_years_with_data$tower_id[i]; yr <- site_years_with_data$year[i]
  has_flight  <- nrow(semi_join(site_years_with_flight, tibble(tower_id = tid, year = yr), by = c("tower_id", "year"))) > 0
  is_matched  <- nrow(semi_join(site_years_matched, tibble(tower_id = tid, year = yr), by = c("tower_id", "year"))) > 0
  reason <- if (is_matched) "MATCHED" else if (!has_flight) "no flight date (no hyperspec coverage that year)" else "no bout within 30 days"
  cat(" ", tid, yr, "->", reason, "\n")
}

no_flight_date <- flight_dates %>% filter(is.na(flight_date))
if (nrow(no_flight_date) > 0) {
  cat("\nTower-years with hyperspec tiles but flight date could not be resolved:\n")
  print(as.data.frame(no_flight_date))
}

# ---- spot check: one real site-year, tower vs. all, at peak_flight --------
spot_check_key <- field_diversity_long %>%
  filter(temporal_scope == "peak_flight") %>%
  count(tower_id, year) %>% filter(n == 2) %>% slice(1)

if (nrow(spot_check_key) == 0) {
  cat("\n==== Spot check: no site-year has a peak_flight match for BOTH plot_scope",
      "values -- nothing to compare side by side. ====\n")
} else {
  sc_tid <- spot_check_key$tower_id[1]; sc_yr <- spot_check_key$year[1]
  cat("\n==== Spot check:", sc_tid, sc_yr, "(temporal_scope = peak_flight) ====\n")

  sc_rows <- field_diversity_long %>%
    filter(tower_id == sc_tid, year == sc_yr, temporal_scope == "peak_flight")
  cat("Matched flight date:", format(sc_rows$flight_date_matched[1]), "\n")
  print(as.data.frame(sc_rows))

  for (ps in plot_scopes) {
    sc_row <- sc_rows %>% filter(plot_scope == ps)
    sc_bout <- sc_row$bout[1]
    cat("\n--- plot_scope =", ps, ", selected bout =", sc_bout,
        ", bout_end_date =", format(sc_row$bout_end_date[1]), "---\n")
    m <- compute_combo_metrics(sc_tid, sc_yr, ps, sc_bout)
    cat("raw d1 (1m2 percentCover) rows:", nrow(m$d1), "\n")
    print(head(m$d1 %>% select(plotID, subplotID, taxonID, percentCover), 10))
    cat("raw dn (nested presence) rows:", nrow(m$dn), "\n")
    print(head(m$dn %>% select(plotID, subplotID, taxonID, additionalSpecies), 10))
    cat("n_plots:", m$n_plots, " floristic_shannon_mean:", m$floristic_shannon_mean,
        " floristic_shannon_gamma:", m$floristic_shannon_gamma,
        " floristic_richness:", m$floristic_richness,
        " max_single_plot_richness:", m$max_single_plot_richness, "\n")
    cat("floristic_shannon_mean_canopy:", m$floristic_shannon_mean_canopy,
        " floristic_shannon_gamma_canopy:", m$floristic_shannon_gamma_canopy,
        " floristic_richness_canopy:", m$floristic_richness_canopy,
        " max_single_plot_richness_canopy:", m$max_single_plot_richness_canopy, "\n")
    if (m$floristic_richness < m$max_single_plot_richness) {
      cat("  !! SANITY FLOOR VIOLATED: gamma richness below a single plot's richness.\n")
    } else {
      cat("  OK: floristic_richness >= max_single_plot_richness (sanity floor holds).\n")
    }
    if (m$floristic_richness_canopy < m$max_single_plot_richness_canopy) {
      cat("  !! SANITY FLOOR VIOLATED (canopy): gamma richness below a single plot's richness.\n")
    } else {
      cat("  OK: floristic_richness_canopy >= max_single_plot_richness_canopy (sanity floor holds).\n")
    }
  }

  richness_tower <- sc_rows$floristic_richness[sc_rows$plot_scope == "tower"]
  richness_all   <- sc_rows$floristic_richness[sc_rows$plot_scope == "all"]
  if (length(richness_tower) == 1 && length(richness_all) == 1) {
    cat("\nfloristic_richness(all) =", richness_all, ">= floristic_richness(tower) =", richness_tower,
        "->", if (richness_all >= richness_tower) "OK" else "!! VIOLATED", "\n")
  }
}

# ---- sanity floor, all computed rows -------------------------------------
cat("\n==== Sanity floor check across all computed per_bout/peak_flight rows ====\n")
violations <- 0
computed_rows <- field_diversity_long %>% filter(status == "ok")
for (i in seq_len(nrow(computed_rows))) {
  r <- computed_rows[i, ]
  m <- compute_combo_metrics(r$tower_id, r$year, r$plot_scope, r$bout)
  if (r$floristic_richness < m$max_single_plot_richness) {
    violations <- violations + 1
    cat("  !! VIOLATION:", r$tower_id, r$year, r$plot_scope, r$temporal_scope, "bout", r$bout,
        "-- richness", r$floristic_richness, "< max single-plot richness", m$max_single_plot_richness, "\n")
  }
  if (r$floristic_richness_canopy < m$max_single_plot_richness_canopy) {
    violations <- violations + 1
    cat("  !! VIOLATION (canopy):", r$tower_id, r$year, r$plot_scope, r$temporal_scope, "bout", r$bout,
        "-- richness_canopy", r$floristic_richness_canopy, "< max single-plot richness_canopy",
        m$max_single_plot_richness_canopy, "\n")
  }
}
cat(if (violations == 0) "All rows pass the sanity floor.\n" else paste(violations, "violation(s) found -- see above.\n"))
