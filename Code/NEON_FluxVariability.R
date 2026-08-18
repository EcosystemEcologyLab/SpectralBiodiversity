# ============================================================
# Interannual variability of NEON flux tower data (NEE, GPP, RECO, ET)
# CV = sd / abs(mean), computed per site across years
# ============================================================

library(tidyverse)

data_dir       <- "./Data/NEON_Ameriflux/AnnualData"
footprint_file <- "./Data/NEONsites_Footprints.csv"
output_dir     <- "./Data/NEON_InterannualVariability"
dir.create(output_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# Site metadata
# ------------------------------------------------------------
sites <- read_csv(footprint_file, show_col_types = FALSE) %>%
  select(site = `Site ID`, site_name = `Site Name`, ecosystem = `Veg Type`)

# ------------------------------------------------------------
# Read + select flux variables (one YY csv per site in data_dir)
# ------------------------------------------------------------
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
length(files)

flux_data <- files %>%
  set_names() %>%
  map_dfr(read_csv, show_col_types = FALSE, na = "-9999", .id = "source_file") %>%
  mutate(site = str_extract(basename(source_file), "US-[A-Za-z0-9]+")) %>%
  select(
    site,
    year = TIMESTAMP,
    NEE  = NEE_VUT_REF,
    GPP  = GPP_NT_VUT_REF,     # use GPP_DT_VUT_REF for daytime partitioning instead
    RECO = RECO_NT_VUT_REF,    # use RECO_DT_VUT_REF for daytime partitioning instead
    ET   = LE_F_MDS            # latent heat flux (W/m2), not mm -- see note below
  )

# Optional: convert LE (W/m2) to ET (mm/yr) if you need absolute ET values
# elsewhere. Not needed for CV, since CV is unaffected by constant scaling.
  flux_data <- flux_data %>%
    mutate(ET = ET * 0.0864 / 2.45 * 365.25)   # W/m2 -> mm/yr, lambda = 2.45 MJ/kg

# ------------------------------------------------------------
# Interannual variability (CV) per site
# ------------------------------------------------------------
variability <- flux_data %>%
  left_join(sites, by = "site") %>%
  group_by(site, site_name, ecosystem) %>%
  summarise(
    n_years   = sum(!is.na(NEE)),
    NEE_mean  = mean(NEE, na.rm = TRUE),  NEE_sd  = sd(NEE, na.rm = TRUE),  NEE_CV  = NEE_sd  / abs(NEE_mean),
    GPP_mean  = mean(GPP, na.rm = TRUE),  GPP_sd  = sd(GPP, na.rm = TRUE),  GPP_CV  = GPP_sd  / abs(GPP_mean),
    RECO_mean = mean(RECO, na.rm = TRUE), RECO_sd = sd(RECO, na.rm = TRUE), RECO_CV = RECO_sd / abs(RECO_mean),
    ET_mean   = mean(ET, na.rm = TRUE),   ET_sd   = sd(ET, na.rm = TRUE),   ET_CV   = ET_sd   / abs(ET_mean),
    .groups = "drop"
  )%>%filter(site != "US-xSC")

write_csv(variability, file.path(output_dir, "NEON_interannual_variability.csv"))

# ------------------------------------------------------------
# Long format for plotting
# ------------------------------------------------------------
variability_long <- variability %>%
  select(site, site_name, ecosystem, ends_with("_CV")) %>%
  pivot_longer(ends_with("_CV"), names_to = "variable", values_to = "CV") %>%
  mutate(variable = str_remove(variable, "_CV"))

# ------------------------------------------------------------
# Plot: CV by site
# ------------------------------------------------------------
p_site <- ggplot(variability_long, aes(x = reorder(site, CV), y = CV, fill = ecosystem)) +
  geom_col() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "NEON site", y = "Interannual variability (CV)", fill = "Ecosystem type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p_site
#ggsave(file.path(output_dir, "variability_by_site.png"), p_site, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------
# Plot: CV by ecosystem
# ------------------------------------------------------------
p_ecosystem <- ggplot(variability_long, aes(x = ecosystem, y = CV)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Ecosystem type", y = "Interannual variability (CV)") +
  theme_bw()
p_ecosystem
#ggsave(file.path(output_dir, "variability_by_ecosystem.png"), p_ecosystem, width = 10, height = 7, dpi = 300)

