library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# ========= Prepare the data =========
df_connections <- X_full %>%
  dplyr::select(lon_breed = V1, lat_breed = V2,
                lon_wint  = V3, lat_wint  = V4,
                source) %>%
  na.omit()

# Breeding dataframe
df_breed <- df_connections %>%
  transmute(lon = lon_breed,
            lat = lat_breed,
            season = "breeding",
            source = source)

# Wintering dataframe
df_wint <- df_connections %>%
  transmute(lon = lon_wint,
            lat = lat_wint,
            season = "wintering",
            source = source)

# Combine and create group key
df_long <- bind_rows(df_breed, df_wint) %>%
  mutate(group = paste(season, source, sep = "_"))

# ========= Basemap =========
world <- ne_countries(scale = "medium", returnclass = "sf")

# Manual palette: 2 reds/orange (breeding) + 2 blues (wintering)
cols <- c(
  "breeding_cmr"   = "#b22222",  # dark red
  "breeding_gps"   = "#e6550d",  # orange
  "wintering_cmr"  = "#08306b",  # dark blue
  "wintering_gps"  = "#4292c6"   # light blue
)

# ========= Plot =========
p = ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.3) +
  
  # Connection lines
  geom_segment(data = df_connections,
               aes(x = lon_breed, y = lat_breed,
                   xend = lon_wint, yend = lat_wint),
               color = "grey40", alpha = 0.14) +
  
  # Breeding + wintering points
  geom_point(data = df_long,
             aes(x = lon, y = lat, color = group),
             size = 1.3, alpha = 0.8) +
  
  scale_color_manual(values = cols,
                     name = "Season / Data source",
                     labels = c("Breeding (CMR)", "Breeding (GPS)",
                                "Wintering (CMR)", "Wintering (GPS)")) +
  
  coord_sf(xlim = c(-17, 40), ylim = c(32, 68), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +       # <<==== axes labels
  theme_minimal() +
  theme(legend.position = "bottom")


ggsave("brut.png", p,
       width = 6000, height = 4500, units = "px", dpi = 600, bg = "white")

