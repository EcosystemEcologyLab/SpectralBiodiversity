# ============================================================
# Download NEON Plant Presence and Percent Cover (DP1.10058.001)
# for all sites, all years, and stack into merged CSV tables.
# ============================================================

# install.packages("neonUtilities")  # uncomment if not yet installed
library(neonUtilities)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

dpID       <- "DP1.10058.001"          # Plant presence and percent cover
output_dir <- "./Data/NEON_FieldData"

neon_token <- Sys.getenv("NEON_API_TOKEN")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Download + stack in one step
# ------------------------------------------------------------
# site = "all"          -> every NEON site with this product
# startdate/enddate = NA -> full period of record
# package = "basic"     -> core fields (no expanded QA fields)
# check.size = FALSE    -> skip interactive y/n size prompt so the
#                          script can run unattended; DP1.10058.001
#                          is vegetation survey data, not sensor data,
#                          so full period-of-record/all-sites is not huge.
#                          Set to TRUE if you'd rather confirm first.

veg_data <- loadByProduct(
  dpID       = dpID,
  site       = "all",
  startdate  = NA,
  enddate    = NA,
  package    = "basic",
  token      = neon_token,
  check.size = FALSE
)

# ------------------------------------------------------------
# Write each stacked table to its own CSV
# ------------------------------------------------------------
# `veg_data` is a named list of data frames -- the actual data
# tables plus a few reference tables (variables, validation,
# categoricalCodes, readme). We write all of them out so the
# reference tables are available for interpreting column codes.

for (tbl_name in names(veg_data)) {
  
  tbl <- veg_data[[tbl_name]]
  
  # readme/citation entries are character vectors, not data
  # frames -- skip those, everything else is a table
  if (!is.data.frame(tbl)) next
  
  out_file <- file.path(output_dir, paste0(tbl_name, ".csv"))
  write.csv(tbl, out_file, row.names = FALSE)
  message("Wrote ", out_file, " (", nrow(tbl), " rows)")
}

# ------------------------------------------------------------
# Also keep the full list in the R environment for immediate use
# ------------------------------------------------------------
# e.g. veg_data$div_1m2Data, veg_data$div_10m2Data100m2Data, etc.
# depending on which sub-tables DP1.10058.001 returns.

list2env(veg_data, envir = .GlobalEnv)

names(veg_data)
