# Villaflores Topographic Analysis
# Author: Augusto Salonio
# Date: 12th Jan 2025
# Description: This script performs topographic analysis including DEM processing,
#              landform classification, and terrain metrics calculation

# Load required libraries ------------------------------------------------
library(sp)
library(terra)
library(sf)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(openxlsx)
library(lubridate)
library(parzer)
library(leaflet)
library(landform)
library(geodata)
library(rgeos)
library(maptools)
library(classInt)
library(RColorBrewer)
library(raster)
library(spatialEco)
library(gridExtra)
library(rasterVis)

# Configuration -----------------------------------------------------------
# Set working directory (modify as needed)
base_dir <- "~/Villaflores/R"
setwd(file.path(base_dir, "Puntos"))

# Define data paths
data_paths <- list(
  limits = file.path(base_dir, "data/Límites"),
  dem = file.path(base_dir, "data/DEM/CEM_V3_20170619_R15_E07_TIF"),
  old_limits = file.path(base_dir, "data/Límites/Viejos"),
  output = file.path(base_dir, "output/Topografico")
)

# Create output directory if it doesn't exist
if (!dir.exists(data_paths$output)) {
  dir.create(data_paths$output, recursive = TRUE)
}

# Data Import -------------------------------------------------------------
cat("Loading spatial data...\n")

# Load boundary limits
lim_6372_file <- file.path(data_paths$limits, "lim_6372.shp")
if (file.exists(lim_6372_file)) {
  lim_6372 <- terra::vect(lim_6372_file)
} else {
  # Load original Villaflores limits and reproject
  villaflores_file <- file.path(data_paths$old_limits, "Villaflores.shp")
  lim_villa <- terra::vect(villaflores_file)
  lim_6372 <- terra::project(lim_villa, "EPSG:6372")
}

# Load DEM data
dem_file <- file.path(data_paths$dem, "Chiapas_r15m.tif")
dem_chiapas <- terra::rast(dem_file)

# Reproject DEM to EPSG:6372
cat("Reprojecting DEM to EPSG:6372...\n")
dem_chiapas_6372 <- terra::project(dem_chiapas, "EPSG:6372")

# Crop DEM to study area
dem_6372 <- terra::crop(dem_chiapas_6372, lim_6372)

# Terrain Analysis -------------------------------------------------------
cat("Calculating terrain metrics...\n")

# Calculate slope (in degrees)
pend <- terra::terrain(dem_6372, "slope", neighbors = 8, unit = "degrees")

# Calculate aspect (orientation in degrees)
orient <- terra::terrain(dem_6372, "aspect", neighbors = 8, unit = "degrees")

# Topographic Position Index (TPI) ---------------------------------------
cat("Calculating Topographic Position Index...\n")

# TPI function with customizable window size
calculate_tpi <- function(dem, window_size = 35) {
  # Create weight matrix
  weight_matrix <- matrix(1/(window_size^2 - 1), 
                          nrow = window_size, 
                          ncol = window_size)
  weight_matrix[ceiling(0.5 * length(weight_matrix))] <- 0
  
  # Calculate focal mean
  focal_mean <- terra::focal(dem, weight_matrix)
  
  # Return TPI (elevation - focal mean)
  return(dem - focal_mean)
}

# Calculate TPI with 35x35 window
tpi_35 <- calculate_tpi(dem_6372, window_size = 35)

# Landform Classification ------------------------------------------------
cat("Classifying landforms...\n")

# Calculate standard deviation of TPI
tpi_values <- terra::values(tpi_35, na.rm = TRUE)
tpi_sd <- sd(tpi_values, na.rm = TRUE)

# Define landform classification matrix based on De Reu et al. (2013) and Weiss (2001)
landform_matrix <- matrix(c(
  -Inf, -tpi_sd, 1,        # Valley
  -tpi_sd, -tpi_sd/2, 2,   # Lower Slope
  -tpi_sd/2, 0, 3,         # Flat Area
  0, tpi_sd/2, 4,          # Middle Slope
  tpi_sd/2, tpi_sd, 5,     # Upper Slope
  tpi_sd, Inf, 6           # Ridge
), ncol = 3, byrow = TRUE)

# Classify landforms
landform <- terra::classify(tpi_35, landform_matrix, right = TRUE)

# Convert to categorical raster with labels
landform <- as.factor(landform)
landform_labels <- data.frame(
  ID = 1:6,
  landform_6372 = c('Valley', 'Lower Slope', 'Flat Area', 
                    'Middle Slope', 'Upper Slope', 'Ridge')
)
levels(landform) <- landform_labels

# Percentage Calculations ------------------------------------------------
cat("Calculating percentage layers...\n")

# Function to calculate percentage in moving window
calculate_percentage <- function(raster_layer, target_values, window_size = 3) {
  # Create binary raster (1 for target values, 0 for others)
  binary_raster <- terra::classify(raster_layer, 
                                   cbind(target_values, rep(1, length(target_values))))
  binary_raster[!values(binary_raster) %in% 1] <- 0
  binary_raster <- as.factor(binary_raster)
  
  # Calculate percentage using focal window
  weight_matrix <- matrix(1, nrow = window_size, ncol = window_size)
  percentage <- terra::focal(binary_raster, 
                             w = weight_matrix,
                             fun = function(x) mean(x, na.rm = TRUE) * 100)
  
  return(percentage)
}

# Calculate valley percentage (landform class 1)
percentage_valley <- calculate_percentage(landform, target_values = 1)

# Calculate ridge percentage (landform class 6)
percentage_ridge <- calculate_percentage(landform, target_values = 6)

# Southwestness Analysis --------------------------------------------------
cat("Calculating southwestness...\n")

# Calculate southwestness (cosine of difference from 225 degrees)
southwestness <- cos((orient - 225) * pi / 180)

# Create binary southwest classification (1 = SW, 0 = other)
sw_6372 <- as.factor(southwestness >= 0)

# Calculate southwest percentage
percentage_sw <- calculate_percentage(sw_6372, target_values = TRUE)

# Clip Data to Study Area ------------------------------------------------
cat("Clipping data to study area...\n")

clip_dem_6372 <- terra::crop(dem_6372, lim_6372, mask = TRUE)
clip_pend_6372 <- terra::crop(pend, lim_6372, mask = TRUE)
clip_tpi_6372 <- terra::crop(tpi_35, lim_6372, mask = TRUE)

# Export Results ----------------------------------------------------------
cat("Exporting results...\n")

# List of rasters to export with their filenames
export_list <- list(
  "dem_6372_30.tif" = clip_dem_6372,
  "dem_entero_6372_30.tif" = dem_6372,
  "pend_6372_30.tif" = clip_pend_6372,
  "tpi_6372_30.tif" = clip_tpi_6372,
  "sw_entero_30.tif" = sw_6372,
  "land_entero_30.tif" = landform,
  "percentage_valley.tif" = percentage_valley,
  "percentage_ridge.tif" = percentage_ridge,
  "percentage_sw.tif" = percentage_sw
)

# Export rasters
for (filename in names(export_list)) {
  filepath <- file.path(data_paths$output, filename)
  terra::writeRaster(export_list[[filename]], filepath, overwrite = TRUE)
  cat("Exported:", filename, "\n")
}

# Export vector data
terra::writeVector(lim_6372, 
                   file.path(data_paths$limits, "lim_6372.shp"), 
                   overwrite = TRUE)

# Create Summary Report ---------------------------------------------------
cat("Creating summary report...\n")

summary_stats <- data.frame(
  Layer = c("DEM", "Slope", "TPI", "Valley %", "Ridge %", "Southwest %"),
  Min = c(
    round(terra::minmax(clip_dem_6372)[1], 2),
    round(terra::minmax(clip_pend_6372)[1], 2),
    round(terra::minmax(clip_tpi_6372)[1], 2),
    round(terra::minmax(percentage_valley)[1], 2),
    round(terra::minmax(percentage_ridge)[1], 2),
    round(terra::minmax(percentage_sw)[1], 2)
  ),
  Max = c(
    round(terra::minmax(clip_dem_6372)[2], 2),
    round(terra::minmax(clip_pend_6372)[2], 2),
    round(terra::minmax(clip_tpi_6372)[2], 2),
    round(terra::minmax(percentage_valley)[2], 2),
    round(terra::minmax(percentage_ridge)[2], 2),
    round(terra::minmax(percentage_sw)[2], 2)
  ),
  Units = c("meters", "degrees", "TPI units", "percent", "percent", "percent")
)

print(summary_stats)

# Save summary to CSV
write.csv(summary_stats, 
          file.path(data_paths$output, "summary_statistics.csv"), 
          row.names = FALSE)

cat("Analysis complete! Check the output directory for results.\n")