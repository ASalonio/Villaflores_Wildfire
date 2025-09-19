# Villaflores Forest Fire Ignition Point Analysis
# Analysis of CONAFOR forest fire clusters and heat point data (2012-2024)
# Author: Augusto Salonio
# Date: 28th Apr 2025
# 
# Description:This script identifies the first ignition point (primer punto de calor) 
#             for each forest fire cluster by intersecting heat points with cluster polygons
#             and selecting the earliest timestamp for each cluster.

# Load required libraries ----
library(sp)
library(terra)
library(sf)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(openxlsx)
library(lubridate)
library(ggplot2)

# Configuration ----
# Set paths - modify these for your environment
BASE_PATH <- "~/Villaflores/R"
DATA_PATH <- file.path(BASE_PATH, "data")
OUTPUT_PATH <- file.path(BASE_PATH, "output")

# Create output directories if they don't exist
dir.create(file.path(OUTPUT_PATH, "Puntos de Ignición - Forestales"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_PATH, "Conglomerados - Forestales"), 
           recursive = TRUE, showWarnings = FALSE)

# Set working directory
setwd(file.path(BASE_PATH, "Puntos"))

# Define analysis years
ANALYSIS_YEARS <- 2012:2024

# Helper Functions ----

#' Load and preprocess cluster data for a given year
#' @param year Integer year to process
#' @param data_path Base path to data directory
#' @param villaflores_limits SpatVector of study area limits
#' @return Preprocessed cluster data clipped to study area
load_cluster_data <- function(year, data_path, villaflores_limits) {
  file_path <- file.path(data_path, "Probabilidad de Ignición/Cong/Forestales", 
                         as.character(year), paste0(year, ".shp"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found for year", year, ":", file_path))
    return(NULL)
  }
  
  # Load and preprocess cluster data
  clusters <- terra::vect(file_path) %>%
    tidyterra::rename_with(~c("FEC_C"), c(FECHA)) %>%
    tidyterra::select(ID_CPC, FEC_C, HA)
  
  # Validate and fix geometries
  invalid_count <- sum(!is.valid(clusters))
  if (invalid_count > 0) {
    message(paste("Fixed", invalid_count, "invalid geometries for year", year))
    clusters <- makeValid(clusters)
  }
  
  # Clip to study area
  clipped_clusters <- terra::crop(clusters, villaflores_limits)
  
  # Report spatial relationships
  overlaps <- sum(relate(clipped_clusters, relation = "overlaps"))
  touches <- sum(relate(clipped_clusters, relation = "touches"))
  covered <- sum(relate(clipped_clusters, relation = "coveredby"))
  
  message(sprintf("Year %d: %d clusters (%d overlaps, %d touches, %d covered)", 
                  year, nrow(clipped_clusters), overlaps, touches, covered))
  
  return(clipped_clusters)
}

#' Load and preprocess heat point data for a given year
#' @param year Integer year to process
#' @param data_path Base path to data directory
#' @param villaflores_limits SpatVector of study area limits
#' @param buffer_width Buffer width in meters for heat points
#' @return List containing clipped points and buffered points
load_heat_data <- function(year, data_path, villaflores_limits, buffer_width = 0.5) {
  file_path <- file.path(data_path, "Probabilidad de Ignición/Calor/Forestales", 
                         as.character(year), paste0(year, ".shp"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found for year", year, ":", file_path))
    return(NULL)
  }
  
  # Load heat points
  heat_points <- terra::vect(file_path)
  
  # Create datetime field
  heat_points_dated <- heat_points %>%
    tidyterra::mutate(Date = as.POSIXlt(paste(FECHA, HORA_GMT, sep = " "), tz = "GMT"))
  
  # Add unique ID
  heat_points_dated$ID <- seq.int(nrow(heat_points_dated))
  
  # Validate geometries
  invalid_count <- sum(!is.valid(heat_points_dated))
  if (invalid_count > 0) {
    message(paste("Found", invalid_count, "invalid heat point geometries for year", year))
  }
  
  # Clip to study area
  clipped_points <- terra::intersect(heat_points_dated, villaflores_limits)
  
  # Create buffer for intersection analysis
  # Note: Some years use different buffer widths for better coverage
  special_buffers <- list("2019" = 1000, "2023" = 2000)
  actual_buffer <- ifelse(as.character(year) %in% names(special_buffers),
                          special_buffers[[as.character(year)]], 
                          buffer_width)
  
  buffered_points <- terra::buffer(clipped_points, width = actual_buffer)
  
  message(sprintf("Year %d: %d heat points (buffer: %g m)", 
                  year, nrow(clipped_points), actual_buffer))
  
  return(list(points = clipped_points, buffered = buffered_points))
}

#' Find ignition points by intersecting heat points with clusters
#' @param heat_data List containing heat points and buffered points
#' @param cluster_data Cluster polygons
#' @return SpatVector of ignition points (earliest per cluster)
find_ignition_points <- function(heat_data, cluster_data) {
  if (is.null(heat_data) || is.null(cluster_data)) {
    return(NULL)
  }
  
  # Intersect buffered heat points with clusters
  intersected <- terra::intersect(heat_data$buffered, cluster_data)
  
  if (nrow(intersected) == 0) {
    warning("No intersections found between heat points and clusters")
    return(NULL)
  }
  
  # Find centroids of intersected areas
  centroids <- terra::centroids(intersected)
  
  # Select earliest point per cluster and clean up variables
  ignition_points <- centroids %>%
    select(IDPUNTOCAL, ID, LATLON, Date, ID_CPC, HA) %>%
    group_by(ID_CPC) %>%
    arrange(Date) %>%
    slice_head() %>%
    ungroup() %>%
    mutate(Date = as.character(Date))
  
  return(ignition_points)
}

# Main Analysis ----

# Load study area boundaries
villaflores_path <- file.path(DATA_PATH, "Límites", "Villaflores.shp")
if (!file.exists(villaflores_path)) {
  stop("Villaflores boundary file not found: ", villaflores_path)
}

lim_villa <- terra::vect(villaflores_path)

# Check CRS
message("Study area CRS: ")
print(crs(lim_villa, describe = TRUE))

# Initialize storage lists
all_clusters <- list()
all_ignition_points <- list()

# Process each year
message("Processing years: ", paste(ANALYSIS_YEARS, collapse = ", "))

for (year in ANALYSIS_YEARS) {
  message(paste("\n--- Processing year", year, "---"))
  
  # Load cluster data
  clusters <- load_cluster_data(year, DATA_PATH, lim_villa)
  if (is.null(clusters)) next
  
  # Load heat point data  
  heat_data <- load_heat_data(year, DATA_PATH, lim_villa)
  if (is.null(heat_data)) next
  
  # Find ignition points
  ignition_points <- find_ignition_points(heat_data, clusters)
  if (is.null(ignition_points)) next
  
  # Store results
  all_clusters[[as.character(year)]] <- clusters
  all_ignition_points[[as.character(year)]] <- ignition_points
  
  # Report intersection success
  intersections <- sum(is.related(clusters, heat_data$buffered, "intersects"))
  message(sprintf("Successfully matched %d/%d clusters with heat points", 
                  nrow(ignition_points), nrow(clusters)))
}

# Combine all years ----
if (length(all_ignition_points) > 0) {
  # Combine ignition points
  combined_ignition <- do.call(tidyterra::bind_spat_rows, all_ignition_points)
  
  # Combine clusters
  combined_clusters <- do.call(tidyterra::bind_spat_rows, all_clusters)
  
  # Verify data integrity
  message(sprintf("\nFinal results: %d ignition points, %d clusters", 
                  nrow(combined_ignition), nrow(combined_clusters)))
  
  # Reproject to EPSG:6372 (Mexico-specific projection)
  ignition_projected <- terra::project(combined_ignition, "EPSG:6372")
  clusters_projected <- terra::project(combined_clusters, "EPSG:6372")
  
  # Export results
  ignition_output <- file.path(OUTPUT_PATH, "Puntos de Ignición - Forestales", 
                               "Seleccion_2012-2024.shp")
  clusters_output <- file.path(OUTPUT_PATH, "Conglomerados - Forestales", 
                               "Cong_2012-2024.shp")
  
  writeVector(ignition_projected, ignition_output, overwrite = TRUE)
  writeVector(clusters_projected, clusters_output, overwrite = TRUE)
  
  message(paste("\nResults exported to:"))
  message(paste("- Ignition points:", ignition_output))
  message(paste("- Clusters:", clusters_output))
  
} else {
  stop("No data was successfully processed. Check file paths and data availability.")
}

# Summary Statistics ----
if (exists("combined_ignition") && exists("combined_clusters")) {
  # Generate summary by year
  ignition_df <- as.data.frame(combined_ignition)
  ignition_df$year <- as.numeric(format(as.Date(ignition_df$Date), "%Y"))
  
  yearly_summary <- ignition_df %>%
    group_by(year) %>%
    summarise(
      n_ignitions = n(),
      total_area_ha = sum(HA, na.rm = TRUE),
      avg_area_ha = mean(HA, na.rm = TRUE),
      .groups = 'drop'
    )
  
  message("\n--- YEARLY SUMMARY ---")
  print(yearly_summary)
  
  # Save summary
  summary_path <- file.path(OUTPUT_PATH, "yearly_summary.csv")
  write.csv(yearly_summary, summary_path, row.names = FALSE)
  message(paste("Summary saved to:", summary_path))
}

message("\n--- ANALYSIS COMPLETE ---")