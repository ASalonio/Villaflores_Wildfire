# Villaflores Demographic Analysis
# Author: Augusto Salonio based on Lado Liñares, Marcos (marcos.lado@udc.es)
# Date: 15th Jan 2025
# Description: Processes road networks and census data, performs spatial interpolation
#              of population density and population change for fire risk analysis

# Load required libraries ------------------------------------------------
library(sp)
library(terra)
library(tidyterra)
library(sf)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(openxlsx)
library(lubridate)
library(parzer)
library(leaflet)
library(gstat)
library(rgdal)
library(raster)
library(ithir)
library(plotKML)
library(dismo)
library(rgeos)
library(proj4)
library(spThin)
library(carData)
library(car)
library(corrplot)
library(ggplot2)
library(MASS)

# Configuration -----------------------------------------------------------  
# Set base directory (modify as needed)
base_dir <- "~/Villaflores/R"
setwd(file.path(base_dir, "Puntos"))

# Define data paths
data_paths <- list(
  limits = file.path(base_dir, "data/Límites"),
  roads = file.path(base_dir, "data/Humano/Caminos 2023"),
  census_2020 = file.path(base_dir, "data/Humano/Censo 2020/iter_07_cpv2020/conjunto_de_datos"),
  census_2010 = file.path(base_dir, "data/Humano/Censo 2010/iter_07_cpv2010/conjunto_de_datos"),
  output_roads = file.path(base_dir, "output/Humano/Caminos"),
  output_census = file.path(base_dir, "output/Humano/Censos"),
  output_ok = file.path(base_dir, "output/Humano/Censos/OK"),
  reference_grid = file.path(base_dir, "output/Malla"),
  output_distance = file.path(base_dir, "output/Humano/Censos/Distancia")
)

# Create output directories if they don't exist
output_dirs <- c(data_paths$output_roads, data_paths$output_census, data_paths$output_ok, data_paths$output_distance)
for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Road network parameters based on Mexican Official Standard NOM-086-SCT2-2015
road_classification <- list(
  highway_speed = 60,     # Carreteras (highways) >= 60 km/h
  local_speed_min = 40,   # Vecinal (local roads) 40-60 km/h
  local_speed_max = 60,
  urban_speed = 40        # Zonas habitadas (urban areas) < 40 km/h
)

# Spatial analysis parameters
spatial_config <- list(
  road_resolution = 50,     # meters for road density calculation
  census_resolution = 1000, # meters for population interpolation
  target_crs = "EPSG:6372"
)

# Load Boundary Data ------------------------------------------------------
cat("Loading boundary data...\n")

# Load Villaflores boundary
villaflores_file <- file.path(data_paths$limits, "Villaflores.shp")
if (!file.exists(villaflores_file)) {
  stop("Villaflores boundary file not found: ", villaflores_file)
}
lim_villa <- terra::vect(villaflores_file)

# Load projected boundary (create if doesn't exist)
lim_6372_file <- file.path(data_paths$limits, "lim_6372.shp")
if (file.exists(lim_6372_file)) {
  lim_6372 <- terra::vect(lim_6372_file)
} else {
  lim_6372 <- terra::project(lim_villa, spatial_config$target_crs)
  terra::writeVector(lim_6372, lim_6372_file)
}

# ROAD NETWORK ANALYSIS ---------------------------------------------------
cat("Processing road network data...\n")

# Load national road network
road_file <- file.path(data_paths$roads, "rvrnc23gw.shp")
if (!file.exists(road_file)) {
  stop("Road network file not found: ", road_file)
}
red_nac <- terra::vect(road_file)

# Clip roads to study area
red_villa <- terra::intersect(red_nac, lim_villa)

# Reproject to target CRS
red_6372 <- terra::project(red_villa, spatial_config$target_crs)

# Fix speed variable for processing
red_6372$VELOCIDAD <- as.integer(red_6372$VELOCIDAD)

# Classify roads according to Mexican Official Standard NOM-086-SCT2-2015
cat("Classifying roads by speed limits...\n")

# Function to classify roads
classify_roads <- function(road_network, config) {
  list(
    highways = road_network[road_network$VELOCIDAD >= config$highway_speed, ],
    local = road_network[road_network$VELOCIDAD >= config$local_speed_min & 
                           road_network$VELOCIDAD < config$local_speed_max, ],
    urban = road_network[road_network$VELOCIDAD < config$urban_speed, ]
  )
}

road_classes <- classify_roads(red_6372, road_classification)

# Export road classifications
cat("Exporting road classifications...\n")
road_exports <- list(
  "red_6372.shp" = red_6372,
  "Car_6372.shp" = road_classes$highways,    # Carreteras
  "Vec_6372.shp" = road_classes$local,       # Vecinal
  "Hab_6372.shp" = road_classes$urban        # Habitadas
)

for (filename in names(road_exports)) {
  filepath <- file.path(data_paths$output_roads, filename)
  terra::writeVector(road_exports[[filename]], filepath, overwrite = TRUE)
}

# Calculate Road Density --------------------------------------------------
cat("Calculating road density...\n")

# Function to calculate road density
calculate_road_density <- function(road_vector, boundary, resolution) {
  # Create raster template
  template <- terra::rast(terra::ext(boundary), resolution = resolution)
  
  # Calculate length density (km per cell)
  road_length_raster <- terra::rasterizeGeom(road_vector, template, 
                                             fun = "length", unit = "km")
  
  # Convert cell area to km²
  cell_area_km2 <- (resolution / 1000)^2
  
  # Normalize by cell area to get density (km/km²)
  road_density <- road_length_raster / cell_area_km2
  
  return(road_density)
}

# Calculate Density for each road class
road_densities <- list()
for (road_type in names(road_classes)) {
  if (nrow(road_classes[[road_type]]) > 0) {
    road_densities[[road_type]] <- calculate_road_density(
      road_classes[[road_type]], 
      lim_6372, 
      spatial_config$road_resolution
    )
  }
}

# Export road density rasters
density_exports <- list(
  "car_dens_1.tif" = road_densities$highways,
  "vec_dens_1.tif" = road_densities$local,
  "hab_dens_1.tif" = road_densities$urban
)

for (filename in names(density_exports)) {
  if (!is.null(density_exports[[filename]])) {
    filepath <- file.path(data_paths$output_roads, filename)
    terra::writeRaster(density_exports[[filename]], filepath, overwrite = TRUE)
  }
}

# CENSUS DATA PROCESSING --------------------------------------------------
cat("Processing census data...\n")

# Function to load and process census data
load_census_2020 <- function(file_path) {
  cat("Loading 2020 census data...\n")
  
  census_2020 <- read_csv(file_path) %>%
    select(LOC, NOM_LOC, NOM_MUN, POBTOT, LONGITUD, LATITUD) %>%
    filter(NOM_MUN == "Villaflores" & 
             !LOC %in% c(9998, 9999) & 
             NOM_LOC != "Total del Municipio") %>%
    mutate(
      LOC = as.integer(LOC),
      x = parzer::parse_lon(LONGITUD),
      y = parzer::parse_lat(LATITUD)
    ) %>%
    select(LOC, NOM_LOC, POBTOT, x, y) %>%
    filter(!is.na(x) & !is.na(y)) %>%
    distinct(LOC, x, y, .keep_all = TRUE)  # Remove duplicates
  
  return(census_2020)
}

load_census_2010 <- function(file_path, coords_reference) {
  cat("Loading 2010 census data...\n")
  
  # Load 2010 census
  census_2010_raw <- read.csv(file_path) %>%
    select(loc, nom_loc, nom_mun, pobtot, longitud, latitud) %>%
    filter(nom_mun == "Villaflores") %>%
    mutate(loc = as.integer(loc))
  
  # Get coordinates for localities present in both censuses
  census_2010_with_coords <- census_2010_raw %>%
    filter(!is.na(latitud) & !is.na(longitud)) %>%
    semi_join(coords_reference, by = c("loc" = "LOC")) %>%
    left_join(coords_reference, by = c("loc" = "LOC")) %>%
    select(loc, nom_loc, pobtot, x, y)
  
  # Process localities missing coordinates
  missing_coords <- census_2010_raw %>%
    filter(!is.na(latitud) & !is.na(longitud)) %>%
    anti_join(coords_reference, by = c("loc" = "LOC")) %>%
    mutate(
      # Convert coordinate format (DDMM to DD.MM)
      latitud = as.numeric(gsub("^(\\d{2})", "\\1.", as.character(latitud))),
      longitud = as.numeric(gsub("^(\\d{2})", "\\1.", as.character(longitud))) * -1
    ) %>%
    select(loc, nom_loc, pobtot, x = longitud, y = latitud)
  
  # Combine all 2010 data
  census_2010_complete <- bind_rows(census_2010_with_coords, missing_coords) %>%
    rename(LOC = loc, NOM_LOC = nom_loc, POB_10 = pobtot) %>%
    filter(!is.na(x) & !is.na(y))
  
  return(census_2010_complete)
}

# Load census data
census_2020_file <- file.path(data_paths$census_2020, "conjunto_de_datos_iter_07CSV20.csv")
census_2010_file <- file.path(data_paths$census_2010, "iter_07_cpv2010.csv")

if (!file.exists(census_2020_file)) {
  stop("2020 Census file not found: ", census_2020_file)
}
if (!file.exists(census_2010_file)) {
  stop("2010 Census file not found: ", census_2010_file)
}

cen_20_villa <- load_census_2020(census_2020_file)
coords_ref <- cen_20_villa %>% select(LOC, x, y)
cen_10_villa <- load_census_2010(census_2010_file, coords_ref)

# Combine census data
cat("Combining census datasets...\n")
combined_census <- cen_20_villa %>%
  rename(POB_20 = POBTOT) %>%
  full_join(cen_10_villa, by = c("LOC", "x", "y")) %>%
  mutate(
    POB_10 = replace_na(POB_10, 0),
    POB_20 = replace_na(POB_20, 0),
    mean_pob = as.integer(ceiling((POB_10 + POB_20) / 2)),
    var_pob = as.integer(POB_20 - POB_10)
  ) %>%
  distinct(LOC, x, y, .keep_all = TRUE)

# Convert to spatial vector
cen_v <- terra::vect(combined_census, geom = c("x", "y"), crs = "EPSG:4326")
cen_v <- terra::crop(cen_v, lim_villa)
cen_6372 <- terra::project(cen_v, spatial_config$target_crs)

# Export census data
cat("Exporting processed census data...\n")
terra::writeVector(cen_6372, file.path(data_paths$output_census, "cen_6372.shp"), overwrite = TRUE)

census_df <- terra::as.data.frame(cen_6372, geom = "XY")
write.csv(census_df, file.path(data_paths$output_census, "cen_6372.csv"), row.names = FALSE)

# POPULATION DENSITY ANALYSIS ---------------------------------------------
cat("Calculating population density using Voronoi polygons...\n")

# Function to calculate population density using Voronoi polygons
calculate_population_density <- function(census_points, boundary) {
  # Create Voronoi polygons
  thiessen <- terra::voronoi(census_points)
  thiessen <- terra::crop(thiessen, boundary)
  
  # Calculate area of influence
  thiessen$area_km2 <- terra::expanse(thiessen, unit = "km")
  
  # Join areas back to points
  census_with_area <- census_points
  census_with_area$area_km2 <- thiessen$area_km2[match(census_points$LOC, thiessen$LOC)]
  
  # Calculate Density
  census_with_area$mean_dens <- census_with_area$mean_pob / census_with_area$area_km2
  census_with_area$var_dens <- census_with_area$var_pob / census_with_area$area_km2
  
  return(census_with_area)
}

cen_area <- calculate_population_density(cen_6372, lim_6372)

# SPATIAL INTERPOLATION ---------------------------------------------------
cat("Performing spatial interpolation...\n")

# Function to check and transform data distribution
check_and_transform_data <- function(data, variable_name, apply_log = TRUE) {
  cat("Checking distribution for", variable_name, "...\n")
  
  # Original distribution
  cat("Original data summary:\n")
  print(summary(data[[variable_name]]))
  
  if (apply_log && variable_name == "mean_dens") {
    # Log transformation for density
    data$transformed_var <- log10(data[[variable_name]] + 0.1)  # Add small constant to avoid log(0)
    transformed_name <- paste0("log_", variable_name)
  } else if (apply_log && variable_name == "var_dens") {
    # Handle negative values in population change
    shift_constant <- abs(min(data[[variable_name]], na.rm = TRUE)) + 0.1
    data$transformed_var <- log(data[[variable_name]] + shift_constant)
    transformed_name <- paste0("log_", variable_name)
  } else {
    data$transformed_var <- data[[variable_name]]
    transformed_name <- variable_name
  }
  
  # Check for NA values
  na_count <- sum(is.na(data$transformed_var))
  if (na_count > 0) {
    warning(paste("Found", na_count, "NA values in", transformed_name))
  }
  
  cat("Transformed data summary:\n")
  print(summary(data$transformed_var))
  
  return(list(data = data, var_name = transformed_name))
}

# Function to perform kriging interpolation
perform_kriging <- function(data, formula_var, boundary, resolution, output_prefix) {
  cat("Performing kriging interpolation for", formula_var, "...\n")
  
  # Create output raster template
  extent <- terra::ext(boundary)
  output_raster <- terra::rast(extent, res = resolution, crs = terra::crs(data))
  
  # Convert to data frame for gstat
  data_df <- terra::as.data.frame(data, geom = "XY")
  
  # Fit variogram
  formula_obj <- as.formula(paste(formula_var, "~ 1"))
  vgm_model <- gstat::variogram(formula_obj, locations = ~x + y, data = data_df)
  
  # Fit spherical model
  fitted_model <- gstat::fit.variogram(vgm_model, gstat::vgm("Sph"))
  
  # Create kriging model
  kriging_model <- gstat::gstat(formula = formula_obj, 
                                locations = ~x + y, 
                                data = data_df, 
                                model = fitted_model)
  
  # Perform interpolation
  interpolated <- terra::interpolate(output_raster, kriging_model)
  
  # Cross-validation
  cross_val <- gstat::krige.cv(formula = formula_obj, 
                               locations = ~x + y, 
                               data = data_df, 
                               model = fitted_model, 
                               nfold = nrow(data_df))
  
  # Calculate validation metrics
  correlation <- cor(cross_val$var1.pred, cross_val$observed, use = "complete.obs")
  mspe <- mean(cross_val$residual^2, na.rm = TRUE)
  mean_error <- mean(cross_val$residual, na.rm = TRUE)
  
  cat("Cross-validation results for", formula_var, ":\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("MSPE:", round(mspe, 6), "\n")
  cat("Mean Error:", round(mean_error, 6), "\n")
  
  # Export results
  pred_file <- file.path(data_paths$output_ok, paste0(output_prefix, "_pred_1km.tif"))
  var_file <- file.path(data_paths$output_ok, paste0(output_prefix, "_var_1km.tif"))
  
  terra::writeRaster(interpolated[[1]], pred_file, overwrite = TRUE)
  terra::writeRaster(interpolated[[2]], var_file, overwrite = TRUE)
  
  return(list(
    interpolated = interpolated,
    validation = list(correlation = correlation, mspe = mspe, mean_error = mean_error),
    variogram = list(empirical = vgm_model, fitted = fitted_model)
  ))
}

# Interpolate population density
density_processed <- check_and_transform_data(cen_area, "mean_dens", apply_log = TRUE)
cen_area_processed <- density_processed$data

density_results <- perform_kriging(
  cen_area_processed, 
  "transformed_var", 
  lim_6372, 
  spatial_config$census_resolution,
  "OK_dens"
)

# Interpolate population change
change_processed <- check_and_transform_data(cen_area, "var_dens", apply_log = TRUE)
cen_area_change <- change_processed$data

change_results <- perform_kriging(
  cen_area_change, 
  "transformed_var", 
  lim_6372, 
  spatial_config$census_resolution,
  "OK_change"
)

# Export processed census data with areas
terra::writeVector(cen_area, file.path(data_paths$output_census, "cen_area.shp"), overwrite = TRUE)

# EUCLIDEAN DISTANCE CALCULATION -----------------------------------------
cat("Calculating Euclidean distance to population centers...\n")

# Create raster template for distance calculation
molde_30 <- terra::rast(terra::ext(lim_6372), resolution = 30, crs = spatial_config$target_crs)

# Calculate Euclidean distance
euc_pob_6372 <- terra::distance(molde_30, cen_6372)

# Clip to study area
clip_pob_6372 <- terra::crop(euc_pob_6372, lim_6372, mask = TRUE)

# Export distance raster
terra::writeRaster(clip_pob_6372, file.path(data_paths$output_distance, "euc_pob_30.tif"), overwrite = TRUE)

# CREATE SUMMARY REPORT ---------------------------------------------------
cat("Creating summary report...\n")

summary_report <- data.frame(
  Analysis = c("Road Network", "Census Data", "Population Density", "Population Change", "Euclidean Distance"),
  Status = c("Complete", "Complete", "Complete", "Complete", "Complete"),
  Records = c(
    nrow(red_6372),
    nrow(combined_census),
    nrow(cen_area),
    nrow(cen_area),
    NA
  ),
  Validation_Correlation = c(
    NA,
    NA,
    round(density_results$validation$correlation, 3),
    round(change_results$validation$correlation, 3),
    NA
  ),
  Notes = c(
    paste("Roads classified:", length(road_classes)),
    "Combined 2010-2020 census",
    "Kriging interpolation with log transformation",
    "Kriging interpolation with shift+log transformation",
    "Distance to population centers at 30m resolution"
  )
)

print(summary_report)
write.csv(summary_report, file.path(data_paths$output_census, "analysis_summary.csv"), row.names = FALSE)

cat("\n=== DEMOGRAPHIC ANALYSIS COMPLETE ===\n")
cat("Road density rasters exported to:", data_paths$output_roads, "\n")
cat("Census interpolation results exported to:", data_paths$output_ok, "\n")
cat("Euclidean distance raster exported to:", data_paths$output_distance, "\n")
cat("Check analysis_summary.csv for detailed results.\n")