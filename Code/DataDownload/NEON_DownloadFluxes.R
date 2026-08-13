# download_footprint_data.R
# Downloads DP4.00200.001 (Bundled eddy covariance data), expanded package,
# for all NEON sites in NEONsites.csv. The expanded package is required --
# it's the one that contains the half-hourly footprint arrays.

library(neonUtilities)
library(dplyr)
library(stringr)

# ---- 1. Setup ---------------------------------------------------------
out_dir    <- "X:/moore/SpectralBiodiversity/Data/NEON_Flux"
dpID       <- "DP4.00200.001"
NEON_TOKEN <- Sys.getenv("NEON_API_TOKEN")

#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load sites and extract NEON 4-letter site codes ---------------
towers_df <- read.csv("./Data/NEONsites.csv") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites to process.\n")

# ---- 3. Download expanded EC bundle for each site, all available dates -
log <- tibble(tower_id = character(), neon_site = character(), status = character())

for (i in seq_len(nrow(towers_df))) {
  
  tower_id  <- towers_df$Site.ID[i]
  neon_site <- towers_df$neon_site[i]
  site_out_dir <- file.path(out_dir, tower_id)
  dir.create(site_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n==== [", i, "/", nrow(towers_df), "] ", tower_id,
      " (", neon_site, ") ====\n", sep = "")
  
  result <- tryCatch({
    zipsByProduct(
      dpID       = dpID,
      package    = "expanded",   # required for footprint arrays
      site       = neon_site,
      startdate  = NA,           # NA = no lower bound, grab everything available
      enddate    = NA,
      savepath   = site_out_dir,
      check.size = FALSE,
      token      = NEON_TOKEN
    )
    "success"
  }, error = function(e) {
    cat("  -> failed for", tower_id, ":", conditionMessage(e), "\n")
    paste("error:", conditionMessage(e))
  })
  
  log <- add_row(log, tower_id = tower_id, neon_site = neon_site, status = result)
}

write.csv(log, file.path(out_dir, "download_log.csv"), row.names = FALSE)
cat("\nDone. EC bundles saved under:", out_dir, "\n")
cat("See download_log.csv for a per-site summary.\n")