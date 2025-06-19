library(osmdata)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(dplyr)
library(ggspatial)
here::i_am("bf_map.R")

# Get the bounding box for Ouagadougou for a focused query
ouaga_bb <- getbb("Ouagadougou, Burkina Faso")

# Build the query:
# 1. Boundary set to administrative.
# 2. Filter by admin_level; level "8" is common for city boundaries.
# 3. Match the name "Ouagadougou".
ouaga_query <- opq(bbox = ouaga_bb) %>%
  add_osm_feature(key = "boundary", value = "administrative") %>%
  add_osm_feature(key = "name", value = "Ouagadougou")

# Retrieve OSM data as an sf object
ouaga_osm <- osmdata_sf(ouaga_query)

# Extract multipolygons; if none are returned, use polygons instead.
ouaga_boundary <- ouaga_osm$osm_multipolygons
if (is.null(ouaga_boundary) || nrow(ouaga_boundary) == 0) {
  ouaga_boundary <- ouaga_osm$osm_polygons
}

# Ensure valid geometries
ouaga_boundary <- st_make_valid(ouaga_boundary)

bf <- st_read("/Users/lander/Downloads/bfa_adm_igb_20200323_shp/bfa_admbnda_adm3_igb_20200323.shp")
province_df <- bf |>
  filter(ADM2_FR == "Kadiogo")
province_boundary <- st_make_valid(province_df$geometry)

ouaga_df <- bf |>
  filter(ADM3_FR == "Ouagadougou")
ouaga_boundary <- st_make_valid(ouaga_df$geometry)

# Optional: Check what features you obtained (uncomment to print summary)
# print(ouaga_boundary)

# Retrieve Africa boundaries and extract Burkina Faso (for the inset map)
africa <- ne_countries(continent = "Africa", returnclass = "sf")
burkina_faso <- africa %>% filter(name == "Burkina Faso")

# Plot the main map using the refined Ouagadougou boundary
main_map <- ggplot() +
  geom_sf(data = province_boundary, fill = "grey90", color = "grey", size = 1) +
  geom_sf(data = ouaga_boundary, fill = "#FFEABD", color = "grey", size = 1) +
  geom_point(aes(y = 12.372472, x = -1.494194)) +
  geom_text(aes(y = 12.372472, x = -1.494194, label = "1200 Logements"), hjust = 0.5, vjust = 2, size = 3) +
  geom_point(aes(y = 12.433583, x = -1.51075)) +
  geom_text(aes(y = 12.433583, x = -1.51075, label = "Toudbweogo"), hjust = 0.5, vjust = -1, size = 3) +
  geom_text(aes(y = 12.29, x = -1.525), label = expression(underline(bold("OUAGADOUGOU"))), size = 3, hjust = 0.5) +
  annotation_scale(location = "br", width_hint = 0.3, unit_category = "metric", line_width = .5, height = unit(.2, "cm")) +
  annotation_north_arrow(
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    location = "tl", which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid.major = element_line(size = .2)
  )

# Create the inset map of Africa with Burkina Faso highlighted in green
inset_map <- ggplot() +
  geom_sf(data = africa, fill = "gray90", color = "grey") +
  geom_sf(data = burkina_faso, fill = "#C9F4BC", color = "black") +
  geom_sf(data = province_boundary, fill = "red", color = "red") +
  theme_void()

# Combine both maps using cowplot. Adjust x, y, width, height as desired.
final_plot <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.57, y = 0.5, width = 0.5, height = 0.5)

# Display the final plot
print(final_plot)
ggsave("figures/burkina_map.pdf", dpi = 300, width = 5, height = 3.6)
