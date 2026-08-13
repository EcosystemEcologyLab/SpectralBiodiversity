library(neonUtilities)
library(sf)
library(dplyr)
library(purrr)
library(stringr)

# ---- 1. Setup ---------------------------------------------------------
out_dir  <- "X:/moore/SpectralBiodiversity/Data/NEON_Hyperspec"
buffer_m <- 530
dpIDs    <- c("DP3.30006.002", "DP3.30006.001")  # prefer bidirectional; fall back to directional
NEON_TOKEN <- Sys.getenv("NEON_API_TOKEN")        # set this in .Renviron if you have one

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load all sites and extract NEON 4-letter site codes -----------
# "Site Name" looks like "NEON Abby Road (ABBY)" -> pull out "ABBY"
towers_df <- read.csv("./Data/NEONsites.csv") %>%
  mutate(neon_site = str_extract(Site.Name, "(?<=\\()[A-Za-z0-9]{4}(?=\\)\\s*$)")) %>%
  filter(!is.na(neon_site))

cat("Found", nrow(towers_df), "sites to process.\n")

towers_sf <- st_as_sf(towers_df, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)

# assign each tower its correct UTM zone (all CONUS/AK/HI/PR sites -> northern hemisphere)
towers_sf <- towers_sf %>%
  mutate(
    utm_zone = floor((Lon + 180) / 6) + 1,
    utm_crs  = paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m +no_defs")
  )

# ---- 3. For each tower: find available years, then download -----------
# keep a running log so you can see what succeeded/failed/had no data
log <- tibble(tower_id = character(), neon_site = character(),
              year = character(), dpID = character(), status = character())

for (i in seq_len(nrow(towers_sf))) {
  
  tower_id  <- towers_sf$Site.ID[i]
  neon_site <- towers_sf$neon_site[i]
  
  # reproject this tower's point into its own UTM zone
  pt_utm <- st_transform(towers_sf[i, ], crs = towers_sf$utm_crs[i])
  coords <- st_coordinates(pt_utm)
  easting  <- coords[1, "X"]
  northing <- coords[1, "Y"]
  
  cat("\n==== [", i, "/", nrow(towers_sf), "] ", tower_id,
      " (", neon_site, ") ====\n", sep = "")
  
  # find which years have hyperspectral data for this site, across both dpIDs
  years_by_dp <- map(dpIDs, function(dp) {
    info <- tryCatch(getProductInfo(dp), error = function(e) NULL)
    if (is.null(info)) return(NULL)
    site_row <- info$siteCodes %>% filter(siteCode == neon_site)
    if (nrow(site_row) == 0) return(NULL)
    months <- unlist(site_row$availableMonths)
    if (length(months) == 0) return(NULL)
    unique(substr(months, 1, 4))
  })
  names(years_by_dp) <- dpIDs
  
  all_years <- sort(unique(unlist(years_by_dp)))
  if (length(all_years) == 0) {
    cat("  No hyperspectral data found for", neon_site, "- skipping.\n")
    log <- add_row(log, tower_id = tower_id, neon_site = neon_site,
                   year = NA, dpID = NA, status = "no data available")
    next
  }
  
  site_out_dir <- file.path(out_dir, tower_id)
  dir.create(site_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (yr in all_years) {
    # prefer .002, fall back to .001 for this year
    dp_to_use <- dpIDs[map_lgl(dpIDs, ~ yr %in% years_by_dp[[.x]])][1]
    
    cat("  Downloading", yr, "using", dp_to_use, "...\n")
    
    result <- tryCatch({
      byTileAOP(
        dpID       = dp_to_use,
        site       = neon_site,
        year       = yr,
        easting    = easting,
        northing   = northing,
        buffer     = buffer_m,
        check.size = FALSE,
        savepath   = site_out_dir,
        token      = NEON_TOKEN
      )
      "success"
    }, error = function(e) {
      cat("    -> failed for", tower_id, yr, ":", conditionMessage(e), "\n")
      paste("error:", conditionMessage(e))
    })
    
    log <- add_row(log, tower_id = tower_id, neon_site = neon_site,
                   year = yr, dpID = dp_to_use, status = result)
  }
}

# ---- 4. Save a summary log so you know what was actually retrieved ----
write.csv(log, file.path(out_dir, "download_log.csv"), row.names = FALSE)
cat("\nDone. Tiles saved under:", out_dir, "\n")
cat("See download_log.csv for a per-site/year summary.\n")