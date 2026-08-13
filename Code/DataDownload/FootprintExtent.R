#find max radial footprint extent
#plot a figure of the distrib of ffp radii, annotate with max and median

library(dplyr)
library(ggplot2)

ffp_df <- read.csv("./Data/Chu_FFP_Summary.csv")

TowerMax_df <- ffp_df%>%
  group_by(site_id)%>%
  summarize(ffp_day_max = max(fpt_extent_day))

MedianExtent <- median(TowerMax_df$ffp_day_max, na.rm = TRUE)
MaxExtent <- max(TowerMax_df$ffp_day_max, na.rm = TRUE)
n_sites <- nrow(TowerMax_df)

ggplot(TowerMax_df, aes(x = ffp_day_max)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = MedianExtent, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_vline(xintercept = MaxExtent, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  annotate("text", x = MedianExtent, y = Inf, 
           label = paste0("Median = ", round(MedianExtent, 0), " m"),
           color = "black", angle = 90, vjust = -0.5, hjust = 1.1, size = 3.5) +
  annotate("text", x = MaxExtent, y = Inf, 
           label = paste0("Max = ", round(MaxExtent, 0), " m"),
           color = "black", angle = 90, vjust = -0.5, hjust = 1.1, size = 3.5) +
  labs(
    title = "Maximum Daily Footprint Radii (80% Contributing Area)",
    subtitle = paste0("n = ", n_sites, " sites"),
    x = "Maximum Footprint Extent (m)",
    y = "Number of Sites"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

#===============================================================================
#NEON max extents

NEON_IDs <- c("US-xBR", 
              "US-xHA", 
              "US-xKZ", 
              "US-xSR", 
              "US-xCP",
              "US-xDL",
              "US-xKA",
              "US-xRM",
              "US-xWD")

NEON_max_df <- TowerMax_df%>%
  filter(site_id %in% NEON_IDs)

NEON_max = max(NEON_max_df$ffp_day_max, na.rm = T)
