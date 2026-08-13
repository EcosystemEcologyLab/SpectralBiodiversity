# diagnose_placeholder_files.R
# Standalone diagnostic: checks every monthly footprint HDF5 file for a given
# site (or site-year) against NEON's documented "force execution" failure
# mode -- a file with valid HDF5 structure/timestamps but every footprint
# cell filled with NaN (see DP4.00200.001 readme). footRaster() reads these
# without raising an error, so this must be checked by content.
#
# Usage:
#   Rscript diagnose_placeholder_files.R <site_dir> [year]
# Example:
#   Rscript diagnose_placeholder_files.R X:/moore/SpectralBiodiversity/Data/NEON_Flux/ABBY 2025

library(neonUtilities)
library(terra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript diagnose_placeholder_files.R <site_dir> [year]")
site_dir <- args[1]
year_filter <- if (length(args) >= 2) args[2] else NULL

all_files <- list.files(site_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
if (!is.null(year_filter)) {
  all_files <- all_files[grepl(paste0("\\.", year_filter, "-\\d{2}\\."), basename(all_files))]
}
all_files <- sort(all_files)

cat("Found", length(all_files), "file(s) to check.\n\n")

results <- data.frame(file = character(), n_layers = integer(),
                       max_frac_valid = double(), mean_frac_valid = double(),
                       verdict = character(), stringsAsFactors = FALSE)

for (f in all_files) {
  cat("----", basename(f), "----\n")
  out <- tryCatch({
    raw_st <- footRaster(filepath = f, progress = FALSE)
    notna <- terra::global(raw_st, fun = "notNA")$notNA
    n_cell <- terra::ncell(raw_st)
    frac <- notna / n_cell
    list(n_layers = terra::nlyr(raw_st), max_frac = max(frac), mean_frac = mean(frac))
  }, error = function(e) {
    cat("  READ ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(out)) {
    results <- rbind(results, data.frame(file = basename(f), n_layers = NA,
                                          max_frac_valid = NA, mean_frac_valid = NA,
                                          verdict = "read failed"))
    next
  }

  verdict <- if (out$max_frac < 0.01) "PLACEHOLDER (NaN-filled)" else "valid"
  cat("  layers:", out$n_layers,
      " | max non-NA fraction (best layer):", round(out$max_frac, 4),
      " | mean non-NA fraction:", round(out$mean_frac, 4),
      " ->", verdict, "\n\n")

  results <- rbind(results, data.frame(file = basename(f), n_layers = out$n_layers,
                                        max_frac_valid = out$max_frac,
                                        mean_frac_valid = out$mean_frac,
                                        verdict = verdict))
}

cat("\n==== SUMMARY ====\n")
print(results, row.names = FALSE)
cat("\nPlaceholder files found:", sum(results$verdict == "PLACEHOLDER (NaN-filled)", na.rm = TRUE),
    "of", nrow(results), "\n")
