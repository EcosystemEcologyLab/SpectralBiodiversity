# Make a histogram showing the number of sites with a specific frequency of overflights
# NEON sites with NEON hyperspectral flights

library(ggplot2)

data <- read.csv("./Data/NEONsites_Footprints.csv")

ggplot(data, aes(x = NumberOverflights)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", boundary = -0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = "NEON Overflight Frequency Across NEON Towers",
    subtitle = paste0("n = ", nrow(data), " sites"),
    x = "Number of Overflights",
    y = "Number of Sites"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
