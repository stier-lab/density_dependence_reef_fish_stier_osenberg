# File: 7_global_study_distribution_map.R
# Purpose: Plot global distribution of density-dependence studies on a world map
# Author: Adrian Stier
# Date: 2025-07-11

# ----------------------------------------------------------------------------
# 1. Dependencies & Setup
# ----------------------------------------------------------------------------
library(here)                # project‐relative paths
library(dplyr)               # data manipulation
library(sf)                  # spatial data handling
library(ggplot2)             # plotting
library(viridis)             # color scales
library(rnaturalearth)       # world map data
library(rnaturalearthdata)   # supporting data
library(ggspatial)           # scale bar & north arrow

# ensure output directory exists
dir.create(here::here("figures"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Prepare study‐location data
# ----------------------------------------------------------------------------
spots <- all_dat2 %>%
  group_by(long_deci, lat_deci) %>%
  summarise(n_studies = n_distinct(betanls2_raw_cm), .groups = "drop") %>%
  st_as_sf(coords = c("long_deci", "lat_deci"), crs = 4326)

# ----------------------------------------------------------------------------
# 3. Load world basemap
# ----------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326)

# ----------------------------------------------------------------------------
# 4. Build publication‐quality map
# ----------------------------------------------------------------------------
p_map <- ggplot() +
  # draw land
  geom_sf(data = world,
          fill  = "white",
          color = "grey80",
          size  = 0.2) +
  # overlay study‐site bubbles
  geom_sf(data = spots,
          aes(size = n_studies, fill = n_studies),
          shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  # combined size & color legend
  scale_size_continuous(
    name  = "Number of Studies",
    range = c(2, 10)
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Number of Studies"
  ) +
  # add a scale bar (bottom left)
  annotation_scale(
    location   = "bl",
    width_hint = 0.2,
    text_cex   = 0.7,
    line_width = 0.5
  ) +
  # add a north arrow (top right)
  annotation_north_arrow(
    location    = "tr",
    which_north = "true",
    pad_x       = unit(0.5, "cm"),
    pad_y       = unit(0.5, "cm"),
    style       = north_arrow_fancy_orienteering
  ) +
  # fix aspect, set global extent
  coord_sf(
    xlim   = c(-180, 180),
    ylim   = c(-60, 85),
    expand = FALSE
  ) +
  # axis labels & title
  labs(
    title = "Global Distribution of Density-Dependence Studies",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  # clean, white background theme
  theme_minimal(base_size = 14) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major   = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    legend.position    = "right"
  )

# ----------------------------------------------------------------------------
# 5. Render & Save
# ----------------------------------------------------------------------------
print(p_map)

# standard name
ggsave(
  filename = here::here("figures", "global_study_distribution_map.png"),
  plot     = p_map,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)
