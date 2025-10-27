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
tw_sf <- st_as_sf(tw)
tw_main <- tw_sf %>%
  filter(!COUNTYENG %in% tw_sf$COUNTYENG[c(1, 14, 18)])
taiwan_main_box <- c(
  xmin = 119.92, xmax = 122.22, 
  ymin = 21.83, ymax = 25.65)
tw_isle <- st_crop(tw_main, st_bbox(taiwan_main_box, crs = st_crs(tw_sf)))
elv <- rast("output_SRTMGL1.tif")
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
elev_main <- crop(elv_tw, ext(tw_main_box))
elv_main_masked <- mask(elev_main, tw_main)
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
