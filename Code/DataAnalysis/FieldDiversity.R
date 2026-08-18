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

# ============================================================================
# 0. Setup / paths
# ============================================================================
div_1m2_path         <- "./Data/NEON_FieldData/div_1m2Data.csv"
div_nested_path       <- "./Data/NEON_FieldData/div_10m2Data100m2Data.csv"
neonsites_path        <- "./Data/NEONsites.csv"
hyperspec_dir         <- "D:/projects/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
out_csv               <- "D:/projects/moore/SpectralBiodiversity/Data/field_diversity_long.csv"

flight_match_tolerance_days <- 30

required_inputs <- c(div_1m2_path, div_nested_path, neonsites_path)
missing_inputs  <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Required input(s) not found:\n  ", paste(missing_inputs, collapse = "\n  "))
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

  # ---- bout end date, for peak_flight matching ----
  end_dates <- if (nrow(d1) > 0) d1$endDate else dn$endDate
  bout_end_date <- if (length(end_dates) == 0) as.Date(NA) else median(end_dates, na.rm = TRUE)

  list(n_plots = n_plots,
       floristic_shannon_mean = floristic_shannon_mean,
       floristic_shannon_gamma = floristic_shannon_gamma,
       floristic_richness = floristic_richness,
       max_single_plot_richness = max_single_plot_richness,
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
      floristic_richness = NA_integer_, status = "no data")
    next
  }

  bout_rows <- map(bouts, function(b) {
    m <- compute_combo_metrics(tid, yr, ps, b)
    tibble(tower_id = tid, neon_site = nsite, year = yr, plot_scope = ps,
           temporal_scope = "per_bout", bout = b, n_plots = m$n_plots,
           flight_date_matched = as.Date(NA), bout_end_date = m$bout_end_date,
           floristic_shannon_mean = m$floristic_shannon_mean,
           floristic_shannon_gamma = m$floristic_shannon_gamma,
           floristic_richness = m$floristic_richness, status = "ok")
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
         floristic_shannon_gamma, floristic_richness, status)

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
    if (m$floristic_richness < m$max_single_plot_richness) {
      cat("  !! SANITY FLOOR VIOLATED: gamma richness below a single plot's richness.\n")
    } else {
      cat("  OK: floristic_richness >= max_single_plot_richness (sanity floor holds).\n")
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
}
cat(if (violations == 0) "All rows pass the sanity floor.\n" else paste(violations, "violation(s) found -- see above.\n"))
