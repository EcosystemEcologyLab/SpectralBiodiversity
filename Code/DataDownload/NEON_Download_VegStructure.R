# NEON_Download_VegStructure.R
#
# Downloads NEON Vegetation Structure (DP1.10098.001) for the 45 towers in
# NEONsites.csv, matched to the year range the field diversity data
# (div_1m2Data.csv / div_10m2Data100m2Data.csv, from NEON_Download_FieldVeg.R)
# already covers, and stacks the result into per-table CSVs. Written to
# supply the two inputs FieldDiversity.R's canopy-filtering addition expects:
# ./Data/NEON_FieldData/vst_apparentindividual.csv and
# ./Data/NEON_FieldData/vst_mappingandtagging.csv.
#
# Mirrors NEON_Download_FieldVeg.R's loadByProduct()-then-write-every-table
# structure (same DP1.100xx.001 family, same output directory convention),
# and NEON_Download_Hyperspec.R's / FieldDiversity.R's NEONsites.csv
# site-code crosswalk (str_extract the 4-letter NEON code out of
# "Site Name (CODE)"), rather than re-deriving either convention differently.
#
# ==============================================================================
# CAVEAT, READ BEFORE RUNNING: this sandbox has no internet/server access, so
# loadByProduct() has never actually been called against the live NEON API
# for this product -- nothing below is verified against a real response.
#
# The two table names this script (and FieldDiversity.R) rely on --
# "vst_apparentindividual" and "vst_mappingandtagging" -- are NEON's
# documented DP1.10098.001 table names per its public data dictionary, not
# confirmed against a real returned list the way this repo's convention
# otherwise requires (see FieldDiversity.R's STEP 1 H5 investigation,
# ExtractMODIS.R's AppEEARS caveat). Section 5 below prints
# names(veg_structure_data) unconditionally on every run and stop()s, listing
# the actual tables returned, if either expected name is absent -- inspect
# that printout on the first real run before trusting the written CSV paths
# or FieldDiversity.R's downstream logic. Do not assume success silently.
#
# Also unconfirmed: whether DP1.10098.001 requires any auth at all for a
# "basic" package pull. NEON's public API documentation states product
# downloads work token-free (a token only raises rate limits), which is why
# this script -- like NEON_Download_FieldVeg.R and NEON_Download_Hyperspec.R
# before it -- treats NEON_API_TOKEN as optional, read from the environment
# (set it in .Renviron, already gitignored in this repo, if you have one).
# If a live run 401s/403s, that documented assumption needs revisiting.
# ==============================================================================

# install.packages("neonUtilities")  # uncomment if not yet installed
library(neonUtilities)
library(dplyr)
library(stringr)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
dpID             <- "DP1.10098.001"   # Vegetation structure
neonsites_path   <- "./Data/NEONsites.csv"
div_1m2_path     <- "./Data/NEON_FieldData/div_1m2Data.csv"
div_nested_path  <- "./Data/NEON_FieldData/div_10m2Data100m2Data.csv"
output_dir       <- "./Data/NEON_FieldData"

neon_token <- Sys.getenv("NEON_API_TOKEN")

required_inputs <- c(neonsites_path, div_1m2_path, div_nested_path)
missing_inputs  <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Required input(s) not found:\n  ", paste(missing_inputs, collapse = "\n  "),
       "\n\nNEONsites.csv defines the site list this script downloads for;",
       " div_1m2Data.csv/div_10m2Data100m2Data.csv (from NEON_Download_FieldVeg.R)",
       " define the year range Vegetation Structure is matched to. Run",
       " NEON_Download_FieldVeg.R first if the div_ files are missing.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Site list -- same NEONsites.csv crosswalk convention already used in
# NEON_Download_Hyperspec.R and FieldDiversity.R, not re-derived differently.
# ------------------------------------------------------------
towers_df <- read.csv(neonsites_path, fileEncoding = "UTF-8-BOM") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

site_codes <- unique(towers_df$neon_site)
cat("Found", length(site_codes), "sites in NEONsites.csv.\n")

# ------------------------------------------------------------
# 3. Year range -- derived from the field diversity data's own coverage,
# not hardcoded (same "derive the range from what's actually needed"
# philosophy as ExtractMODIS.R, which reads its date range off
# spectral_diversity_by_year.csv rather than a fixed year span).
# ------------------------------------------------------------
div_1m2    <- read.csv(div_1m2_path, fileEncoding = "UTF-8-BOM")
div_nested <- read.csv(div_nested_path, fileEncoding = "UTF-8-BOM")

end_dates <- as.Date(substr(c(div_1m2$endDate, div_nested$endDate), 1, 10))
end_dates <- end_dates[!is.na(end_dates)]
if (length(end_dates) == 0) {
  stop("No usable endDate values found in div_1m2Data.csv / div_10m2Data100m2Data.csv --",
       " cannot derive a year range for the Vegetation Structure download.")
}

start_year <- format(min(end_dates), "%Y")
end_year   <- format(max(end_dates), "%Y")
startdate  <- paste0(start_year, "-01")
enddate    <- paste0(end_year, "-12")

cat("Field diversity data spans", start_year, "-", end_year,
    "-- matching Vegetation Structure download to that range.\n")

# ------------------------------------------------------------
# 4. Download + stack in one step
# ------------------------------------------------------------
# site = site_codes   -> only the 45 NEONsites.csv towers, not every NEON
#                        site carrying this product
# startdate/enddate   -> matched to div_ data's own coverage (Section 3),
#                        not the full period of record
# package = "basic"   -> core fields, same as NEON_Download_FieldVeg.R
# check.size = FALSE  -> unattended run; set TRUE to confirm size first.
#                        Vegetation Structure is per-individual tabular data
#                        (not sensor data), so 45 sites x a few years is not
#                        expected to be huge, but this is unconfirmed --
#                        flip to TRUE for the first real run if in doubt.

cat("\nDownloading", dpID, "for", length(site_codes), "sites,",
    startdate, "to", enddate, "...\n")

veg_structure_data <- loadByProduct(
  dpID       = dpID,
  site       = site_codes,
  startdate  = startdate,
  enddate    = enddate,
  package    = "basic",
  token      = neon_token,
  check.size = FALSE
)

# ------------------------------------------------------------
# 5. Inspect the real returned structure BEFORE trusting it -- mandatory
# per the CAVEAT above, not an optional sanity check.
# ------------------------------------------------------------
cat("\n==== Returned table names (loadByProduct output, ", dpID, ") ====\n", sep = "")
print(names(veg_structure_data))

required_tables <- c("vst_apparentindividual", "vst_mappingandtagging")
missing_tables  <- setdiff(required_tables, names(veg_structure_data))
if (length(missing_tables) > 0) {
  stop("loadByProduct() for ", dpID, " did not return the expected table(s): ",
       paste(missing_tables, collapse = ", "),
       ". Actual tables returned: ", paste(names(veg_structure_data), collapse = ", "),
       ". FieldDiversity.R's canopy-filtering addition expects these exact table",
       " names -- investigate the real structure printed above and update either",
       " this script or FieldDiversity.R's expected paths accordingly. Do not guess.")
}

# ------------------------------------------------------------
# 6. Write each stacked table to its own CSV
# ------------------------------------------------------------
# Same convention as NEON_Download_FieldVeg.R: write every returned data
# frame, not just the two FieldDiversity.R needs, so reference tables
# (variables, validation, categoricalCodes, readme) are kept too.

for (tbl_name in names(veg_structure_data)) {

  tbl <- veg_structure_data[[tbl_name]]

  # readme/citation entries are character vectors, not data frames -- skip
  # those, everything else is a table
  if (!is.data.frame(tbl)) next

  out_file <- file.path(output_dir, paste0(tbl_name, ".csv"))
  write.csv(tbl, out_file, row.names = FALSE)
  cat("Wrote ", out_file, " (", nrow(tbl), " rows)\n", sep = "")
}

cat("\nConfirmed required outputs for FieldDiversity.R's canopy-filtering addition:\n")
for (tbl_name in required_tables) {
  cat("  ", file.path(output_dir, paste0(tbl_name, ".csv")), "\n", sep = "")
}

# ------------------------------------------------------------
# 7. Also keep the full list in the R environment for immediate use
# ------------------------------------------------------------
# e.g. veg_structure_data$vst_apparentindividual, $vst_mappingandtagging, etc.

list2env(veg_structure_data, envir = .GlobalEnv)

names(veg_structure_data)
