#setwd("/Users/guanjie/Desktop/Labwork/Research_Data/")
#setwd("/Users/fs/Downloads/GuanJie/")
#setwd("/Users/guanjie/huanglab.mycoentre@gmail.com - Google Drive/Other computers/My Computer/GuanJie")

if (!require(devtools)) install.packages("devtools")
if (!require(dplyr)) install.packages("dplyr")
if (!require(geodata)) install.packages("geodata")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggspatial)) install.packages('ggspatial')
if (!require(terra)) install.packages('terra')
if (!require(tidyterra)) install.packages('tidyterra')
if (!require(elevatr)) devtools::install_github("jhollist/elevatr")

library(devtools)
library(dplyr)
library(elevatr)
library(geodata)
library(ggplot2)
library(ggspatial)
library(sf)
library(terra)
library(tidyterra)

#############
# FUNCTIONS #
#############

lat_elev_profile <- function(r) {
  n <- nrow(r)
  y <- yFromRow(r, 1:n)
  maxElev <- numeric(n)
  for (i in 1:n) {
    vals <- terra::values(r, row = i, nrows = 1)
    m <- suppressWarnings(max(vals, na.rm = TRUE))
    maxElev[i] <- ifelse(is.finite(m), m, NA)  # replace -Inf with NA
  }
  tibble(latitude = y, maxElev = maxElev)
}

# Function mapping elevation (m) to fake_long
to_fake_long <- function(elev) {
  fake_long_min + (elev / max_elev) * (fake_long_max - fake_long_min)
}

##############
# INPUT DATA #
##############

alt_colors <- colorRampPalette(c("#acd0a5", "#94bf8b", "#a8c68f", "#bdcc96", "#d1d7ab", "#e1e4b5", "#efebc0", "#e8e1b6", "#ded6a3", "#d3ca9d", "#cab982", "#c3a76b", "#b9985a", "#aa8753", "#ac9a7c", "#baae9a", "#cac3b8", "#e0ded8", "#f5f4f2"))

sites <- read.csv2("ero/Collection_data.csv", header = TRUE, sep = ",")
sites <- sites %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  )

tw <- st_transform(st_read("shp/COUNTY_MOI_1140318.shp"), crs = "+proj=longlat +datum=WGS84 +no_defs")
#https://data.moi.gov.tw/MoiOD/System/DownloadFile.aspx?DATA=72874C55-884D-4CEA-B7D6-F60B0BE85AB0
#tw <- gadm("Taiwan", level = 2, version = "latest", resolution = 1, path = getwd())
tw_sf <- st_as_sf(tw)
tw_main <- tw_sf %>%
  filter(!COUNTYENG %in% tw_sf$COUNTYENG[c(1, 14, 18)])
tw_mt <- tw_sf %>%
  filter(COUNTYENG %in% tw_sf$COUNTYENG[1])
tw_km <- tw_sf %>%
  filter(COUNTYENG %in% tw_sf$COUNTYENG[14])
tw_ph <- tw_sf %>%
  filter(COUNTYENG %in% tw_sf$COUNTYENG[18])
tw_isle <- st_crop(tw_main, st_bbox(c(
  xmin = 119.92, xmax = 122.22, 
  ymin = 21.83, ymax = 25.65), crs = st_crs(tw_sf)))
dongsha_isle <- st_crop(tw_main, st_bbox(c(
  xmin = 116.70, xmax = 116.75, 
  ymin = 20.68, ymax = 20.73), crs = st_crs(tw_sf)))
taiwan_box <- c(xmin = 114.3593, xmax = 124.5612, ymin = 10.37135, ymax = 26.38528) 
elv <- rast("output_SRTMGL1.tif") #Taiwan elevation with 1 arc-second resolution
#116.5, 21.4499, 124.5612, 26.38528

#### Run this if you have multiple rasters to combine ####----
elv2 <- rast("output_SRTMGL1-2.tif")
trgt_res <- res(elv2)
trgt_ext <- union(ext(elv), ext(elv2))
trgt_rast <- rast(
  extent = trgt_ext,
  resolution = trgt_res,
  crs = crs(elv)
)
elv_resampled  <- resample(elv,  trgt_rast, method = "bilinear")
elv2_resampled <- resample(elv2, trgt_rast, method = "bilinear")
elv_tw <- cover(elv2_resampled, elv_resampled)
#####----

#############################
# ELEVATION DATA GENERATION #
#############################
#elev_tw <- crop(elv_tw, ext(taiwan_box))
#elv_masked <- mask(elev_tw, tw)
#elvMax_df <- lat_elev_profile(elv_masked)

elev_main <- crop(elv_tw, ext(tw_main_box))
elev_matsu <- crop(elv_tw, ext(tw_mt))
elev_kinmen <- crop(elv_tw, ext(tw_km))
elev_penghu <- crop(elv_tw, ext(tw_ph))
elev_dongsha <- crop(elv_tw, ext(dongsha_box))
elv_main_masked <- mask(elev_main, tw_main)
elv_mt_masked <- mask(elev_matsu, tw_main)
elv_km_masked <- mask(elev_kinmen, tw_main)
elv_ph_masked <- mask(elev_penghu, tw_main)
elv_ds_masked <- mask(elev_dongsha, tw_main)

elvMax_df <- lat_elev_profile(elv_main_masked)
elvMax <- elvMax_df %>%
  mutate(
    latitude = as.numeric(latitude),
    maxElev  = as.numeric(maxElev)
  ) %>%
  filter(!is.na(latitude) & !is.na(maxElev)) %>%
  arrange(latitude)
elvMax_scaled <- elvMax %>%
  mutate(
    fake_long = vert_bound + (maxElev / max_elev) # map elev to long
  )
sites_scaled <- sites %>%
  mutate(
    fake_long = vert_bound + (Elevation / max_elev)
  )

###################################
# ELEVATION PROFILE PLOT SETTINGS #
###################################

vert_bound <- 123  
max_elev   <- 4000 
elev_breaks <- seq(0, 4000, 1000)
lwr_bound_axis <- 122.75
fake_long_min <- min(elvMax_scaled$fake_long, na.rm = TRUE)
fake_long_max <- vert_bound + (max_elev / max_elev)
max_elev <- max(elev_breaks)
elev_pos <- to_fake_long(elev_breaks)
y_axis_line <- 20.35
y_offset <- 0.2
yend <- y_axis_line + y_offset + 5.2

################
# MAP PLOTTING #
################
#ggplot() +
#  geom_spatraster(data = elv_main_masked) + #, maxcell = 3.8e7) +
#  scale_fill_gradientn(
#    colours = alt_colors(100), 
#    na.value = NA, 
#    guide = guide_colorbar(
#      ticks.colour = "#000000"
#    )
#  ) +
#  geom_sf(data = tw_isle, color = "black", fill = NA) +
#  geom_point(
#    data = sites, aes(x = Longitude, y = Latitude), size = 0.8
#  ) +
#  geom_point(
#    data = sites_scaled, aes(x = fake_long, y = Latitude), size = 0.8
#  ) +
#  coord_sf(crs = 4326, xlim = c(119, 124.5), ylim = c(20.0, 26.0), expand = FALSE) +
#  # --- custom elevation axis ---
#  # outer vertical grid line for the elevation profile
#  geom_segment(data = data.frame(x = lwr_bound_axis),
#               aes(x = x,
#                   y = y_axis_line + y_offset, 
#                   yend = yend),   # extend up to plot top (adjust if needed)
#               inherit.aes = FALSE, 
#               color = "#000000", 
#               linewidth = 0.4
#  ) + 
#  #vertical grid lines for the elevation profile
#  geom_segment(data = data.frame(x = elev_pos),
#               aes(x = x, xend = x,
#                   y = y_axis_line + y_offset, 
#                   yend = yend),
#               inherit.aes = FALSE, 
#               color = "#9e9e9e", 
#               linewidth = 0.4,
#               alpha = 0.3
#  ) + 
#  geom_path(data = elvMax_scaled, aes(x = fake_long, y = latitude, group = 1),
#            color = "#3f3f3f", linewidth = 0.4) +
#  # axis line spanning only the fake_long strip
#  geom_segment(aes(x = lwr_bound_axis, xend = fake_long_max,
#                   y = y_axis_line + y_offset, yend = y_axis_line + y_offset),
#               linewidth = 0.4
#  ) +
#  # ticks
#  geom_segment(data = data.frame(x = elev_pos),
#               aes(x = x, xend = x,
#                   y = y_axis_line + y_offset,
#                   yend = y_axis_line + y_offset - 0.05),
#               inherit.aes = FALSE, linewidth = 0.4
#  ) +
#  # labels
#  geom_text(data = data.frame(x = elev_pos, lab = elev_breaks),
#            aes(x = x, y = y_axis_line + y_offset - 0.10, label = lab),
#            inherit.aes = FALSE, size = 2.5
#  ) +
#  # axis title
#  annotate("text",
#           x = mean(c(lwr_bound_axis, fake_long_max)),
#           y = y_axis_line + y_offset - 0.25,
#           label = "Elevation (m)",
#           size = 3.2, fontface = "bold"
#  ) +
#  scale_x_continuous(
#    limits = c(119, 125),
#    breaks = c(119, 120, 121, 122)
#  ) +
#  scale_y_continuous(
#    limits = c(20.3, 26.0),
#    breaks = c(21, 22, 23, 24, 25, 26)
#  ) +
#  annotation_scale(
#    location = "bl",
#    width_hint = 0.3,
#    text_cex = 0.8,
#    line_width = 0.4
#  ) +
#  theme(
#    axis.text = element_text(color = "#000000", face = "bold"),
#    axis.title = element_blank(),
#    legend.background = element_blank(),
#    legend.key.size = unit(0.5, "cm"),
#    legend.position = "inside",
#    legend.position.inside = c(0.00, 1.00),
#    legend.justification = c("left", "top"),
#    legend.box.just = "right",
#    legend.title = element_blank(),
#    panel.background = element_rect(fill = 'transparent'),
#    panel.border = element_rect(fill = 'transparent'),
#    panel.grid.major = element_blank(),
#    plot.background = element_rect(fill = 'transparent', color = NA)
#  )

##############
# Extra Plot #
##############

## === Create the main map ===
#main_map <- ggplot() +
#  geom_spatraster(data = elv_masked) +
#  scale_fill_gradientn(
#    colours = alpha(alt_colors(100), 0.7),
#    na.value = NA,
#    guide = guide_colorbar(ticks.colour = "#000000")
#  ) +
#  geom_sf(data = tw_sf, color = "black", fill = NA) +
#  coord_sf(crs = 4326, xlim = c(117, 124.5), ylim = c(20, 26.8), expand = FALSE) +
#  theme_minimal()
#
## === Create the inset map ===
#dongsha_map <- ggplot() +
#  geom_sf(data = tw_sf, color = "black", fill = "grey90") +
#  geom_sf(data = tw_sf, color = "black", fill = NA) +
#  coord_sf(crs = 4326, xlim = c(116.70, 116.75), ylim = c(20.68, 20.73), expand = FALSE) +
#  theme_void() +
#  theme(
#    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
#  )
#
#kinmen_map <- ggplot() +
#  geom_sf(data = tw_sf, color = "black", fill = "grey90") +
#  geom_sf(data = tw_sf, color = "black", fill = NA) +
#  coord_sf(crs = 4326, xlim = c(118.13, 118.44), ylim = c(24.37, 24.53), expand = FALSE) +
#  theme_void() +
#  theme(
#    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
#  )
#
## === Add inset inside the main plot ===
#main_map +
#  annotation_custom(
#    grob = ggplotGrob(dongsha_map),
#    xmin = 117.1, xmax = 118.1,   # position of inset box inside the map (tweak)
#    ymin = 20.1, ymax = 22.1      # adjust for where you want it placed
#  ) +
#  annotate("rect", xmin = 116, xmax = 117, ymin = 20, ymax = 21,
#           colour = "black", fill = NA, linewidth = 0.4
#  ) +
#  annotation_custom(
#    grob = ggplotGrob(kinmen_map),
#    xmin = 118.2, xmax = 120.2,   # position of inset box inside the map (tweak)
#    ymin = 22.2, ymax = 24.2      # adjust for where you want it placed
#  ) +
#  annotate("rect", xmin = 116, xmax = 117, ymin = 20, ymax = 21,
#           colour = "black", fill = NA, linewidth = 0.4
#  )
