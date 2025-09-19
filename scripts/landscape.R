# Villaflores Comprehensive Landscape and Fire Risk Analysis
# Author: Augusto Salonio
# Date: August 1, 2025
# Description: Unified analysis of land cover, forest fire ignition points, forest-agriculture interface (including deciduous density), and multi-temporal spectral indices for fire risk assessment in Villaflores.

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
library(rgee)
library(RStoolbox)
library(raster)
library(rgeos)
library(exactextractr)
library(tidyterra) # For terra-tidyverse integration

# Configuration -----------------------------------------------------------
# Set base directory (modify as needed)
base_dir <- "~/Villaflores/R"
setwd(file.path(base_dir, "Puntos"))

# Define data paths
data_paths <- list(
  limits = file.path(base_dir, "data/Límites"),
  ignition = file.path(base_dir, "output/Puntos de Ignición - Forestales"),
  mexico_admin = file.path(base_dir, "data/Límites/gadm41_MEX_shp"),
  dem = file.path(base_dir, "data/DEM"),
  vegetation = file.path(base_dir, "data/Vegetación/Landcover/mex_land_cover_2010v2_30m_tif/MEX_Land_cover_2010v2_30m_TIF/MEX_NALCMS_landcover_2010v2_30m/data"),
  satellite_april = file.path(base_dir, "data/Vegetación/Abril"),
  satellite_march = file.path(base_dir, "data/Vegetación/Marzo"),
  output_veg = file.path(base_dir, "output/Vegetación"),
  output_cover = file.path(base_dir, "output/Vegetación/coberturas"),
  output_interface = file.path(base_dir, "output/Vegetación/Interfaz Deciduo - Agricola"),
  output_stacks = file.path(base_dir, "output/Vegetación/stacks"),
  output_indices = file.path(base_dir, "output/Vegetación/indices"),
  output_wdvi = file.path(base_dir, "output/Vegetación/WDVI sobre Percentil 75"),
  output_msavi = file.path(base_dir, "output/Vegetación/MSAVI Percentil 10")
)

# Create output directories
output_dirs <- unlist(data_paths[grepl("output", names(data_paths))])
for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Configuration parameters
config <- list(
  # Land cover codes
  forest_types = list(
    needleleaf_temperate = 1,
    broadleaf_evergreen = 3,
    broadleaf_deciduous = 4,
    temperate_deciduous = 5,
    mixed_forest = 6
  ),
  fire_prone_types = c(1, 4, 5), # Needleleaf, deciduous, montane
  cropland = 15,
  # Spatial parameters
  target_crs = "EPSG:6372",
  reference_resolution = 30,
  analysis_resolution = 50,
  buffer_distance = 65.1, # ~210m influence zone
  # Satellite sensor configurations
  landsat7_bands = c("B1", "B2", "B3", "B4", "B5", "B7"),
  landsat8_bands = c("B2", "B3", "B4", "B5", "B6", "B7"),
  landsat9_bands = c("B2", "B3", "B4", "B5", "B6", "B7"),
  # TOA reflectance parameters for Landsat 8/9
  toa_params = list(
    multiplier = 0.00002,
    additive = -0.1
  )
)

# Satellite image metadata (sun elevation angles for TOA conversion)
satellite_metadata <- list(
  "2014" = list(sensor = "L8", sun_elevation = 64.28770715),
  "2016" = list(sensor = "L8", sun_elevation = 66.87137424),
  "2019" = list(sensor = "L8", sun_elevation = 66.22029779),
  "2020" = list(sensor = "L8", sun_elevation = 66.50779175),
  "2021" = list(sensor = "L8", sun_elevation = 64.80402218),
  "2023" = list(sensor = "L9", sun_elevation = 66.59567898),
  "2024" = list(sensor = "L9", sun_elevation = 64.87346657)
)

# Load Boundary and Input Data --------------------------------------------
cat("Loading boundary and input data...\n")

# Administrative boundaries
tryCatch({
  lim_6372 <- terra::vect(file.path(data_paths$limits, "lim_6372.shp"))
  villa <- terra::vect(file.path(data_paths$limits, "Villaflores.shp"))
  lim_mex <- terra::vect(file.path(data_paths$mexico_admin, "gadm41_MEX_2.shp"))
}, error = function(e) {
  stop("Failed to load boundary data: ", e$message)
})

# DEM data
tryCatch({
  dem_villa_EPSG <- terra::rast(file.path(data_paths$dem, "dem_villa_EPSG.tif"))
}, error = function(e) {
  warning("DEM file not found: ", e$message)
  dem_villa_EPSG <- NULL
})

# Fire ignition points
ignition_file <- file.path(data_paths$ignition, "Seleccion_2012-2024.shp")
if (file.exists(ignition_file)) {
  ig_6372 <- terra::vect(ignition_file)
  cat("Loaded", nrow(ig_6372), "ignition points\n")
} else {
  ig_6372 <- NULL
  warning("Ignition points file not found: ", ignition_file)
}

# PART 1: LAND COVER ANALYSIS --------------------------------------------
cat("\n=== PART 1: LAND COVER ANALYSIS ===\n")

#Land cover legend:
#4: Tropical or sub-tropical broadleef deciduous forest (Bosque tropical caducifolio o selva baja caducifolia)
#9: Tropical or sub-tropical grassland (Pastizal)
#15: Cropland (Agricultura)
#1: Temperate or sub-polar needleleaf forest (Bosque de coníferas)
#8: Temperate or sub-polar shrubland (Bosque espinoso)
#5: Temperate or sub-polar broadleaf deciduous forest (Bosque mesófilo de montaña o Bosque deciduo templado)
#6: Mixed forest (Bosque encinar-pinar)
#3: Tropical or sub-tropical broadleaf evergreen forest (selva perennifolia)
#17: Urban (Urbano)
#18: Water (Agua)

# Function to load and process land cover data
load_landcover_data <- function(file_path, boundary, target_crs) {
  if (!file.exists(file_path)) {
    stop("Land cover file not found: ", file_path)
  }
  landcover <- terra::rast(file_path)
  landcover_4326 <- terra::project(landcover, "EPSG:4326")
  landcover_cropped <- terra::crop(landcover_4326, boundary, mask = TRUE)
  landcover_projected <- terra::project(landcover_cropped, target_crs)
  return(landcover_projected)
}

# Load 2010 land cover
landcover_2010_file <- file.path(data_paths$vegetation, "MEX_NALCMS_landcover_2010v2_30m.tif")
cover_10_6372 <- load_landcover_data(landcover_2010_file, villa, config$target_crs)

# Export processed land cover
writeRaster(cover_10_6372, file.path(data_paths$output_veg, "cover_10_villa_EPSG6372.tif"), overwrite = TRUE)

# Analyze land cover composition
analyze_landcover <- function(landcover_raster) {
  freq_table <- terra::freq(landcover_raster, digits = 0)
  freq_table <- freq_table[!is.na(freq_table$value), ]
  total_pixels <- sum(freq_table$count)
  freq_table$percentage <- round((freq_table$count / total_pixels) * 100, 1)
  # Add descriptions
  freq_table$description <- sapply(freq_table$value, function(x) {
    switch(as.character(x),
           "1" = "Needleleaf Forest", "3" = "Evergreen Forest",
           "4" = "Deciduous Forest", "5" = "Temperate Deciduous",
           "6" = "Mixed Forest", "8" = "Shrubland", "9" = "Grassland",
           "15" = "Cropland", "17" = "Urban", "18" = "Water",
           paste("Unknown (", x, ")"))
  })
  return(freq_table[order(-freq_table$percentage), ])
}

composition_2010 <- analyze_landcover(cover_10_6372)
cat("Land cover composition:\n")
print(composition_2010)

# Forest classification
extract_forest_types <- function(landcover_raster, forest_codes, output_name) {
  forest_raster <- terra::classify(landcover_raster, cbind(forest_codes, forest_codes))
  forest_raster[!terra::values(forest_raster) %in% forest_codes] <- 0
  names(forest_raster) <- output_name
  return(forest_raster)
}

all_forest_codes <- unlist(config$forest_types)
forest_6372 <- extract_forest_types(cover_10_6372, all_forest_codes, "forest")
fire_prone_6372 <- extract_forest_types(cover_10_6372, config$fire_prone_types, "fire_prone")

# Analyze fire ignition points in forest areas
if (!is.null(ig_6372)) {
  ig_forest <- terra::extract(forest_6372, ig_6372) %>%
    summarise(N = n(), .by = forest)
  cat("Fire ignition points by forest type:\n")
  print(ig_forest)
}

#135 ignitions in land cover 1 = Temperate or sub-polar needleleaf forest = bosque de coniferas.
#115 ignitions in land cover 4 = Tropical or sub-tropical broadleef deciduous forest = selva baja caducifolia.
#41 ignitions in land cover 5 = Temperate or sub-polar broadleaf deciduous forest = bosque mesófilo de montaña.

#Land cover 3 = Tropical or sub-tropical broadleaf evergreen forest and Land cover 6 = Mixed forest have no ignitions.

# Export forest classifications
writeRaster(forest_6372, file.path(data_paths$output_cover, "forest.tif"), overwrite = TRUE)
writeRaster(fire_prone_6372, file.path(data_paths$output_cover, "des.tif"), overwrite = TRUE)

# PART 2: FOREST-AGRICULTURE INTERFACE AND DENSITY ANALYSIS ---------------
cat("\n=== PART 2: FOREST-AGRICULTURE INTERFACE AND DENSITY ANALYSIS ===\n")

# Create forest-agriculture interface
create_interface <- function(landcover_raster, forest_codes, agriculture_code, buffer_distance, target_crs) {
  # Create binary rasters
  forest_binary <- terra::ifel(landcover_raster %in% forest_codes, 1, NA)
  agriculture_binary <- terra::ifel(landcover_raster == agriculture_code, agriculture_code, NA)
  
  # Create buffers
  forest_buffer <- terra::buffer(forest_binary, width = buffer_distance)
  agriculture_buffer <- terra::buffer(agriculture_binary, width = buffer_distance)
  
  # Identify interface zones
  forest_near_ag <- terra::ifel(forest_binary == 1 & !is.na(agriculture_buffer), 115, NA)
  ag_near_forest <- terra::ifel(agriculture_binary == agriculture_code & !is.na(forest_buffer), 151, NA)
  
  return(list(forest_near_ag = forest_near_ag, ag_near_forest = ag_near_forest))
}

interface_zones <- create_interface(cover_10_6372, config$fire_prone_types, config$cropland, config$buffer_distance, config$target_crs)

# Export interface rasters
if (!is.null(interface_zones)) {
  writeRaster(interface_zones$forest_near_ag, file.path(data_paths$output_interface, "interfaz1con15_des_6372.tif"), overwrite = TRUE)
  writeRaster(interface_zones$ag_near_forest, file.path(data_paths$output_interface, "interfaz15con1_des_6372.tif"), overwrite = TRUE)
}

# Calculate deciduous forest-agriculture interface density
calculate_interface_density <- function(interface_zones, boundary, resolution, target_crs) {
  if (is.null(interface_zones)) {
    warning("No interface zones provided for density calculation")
    return(NULL)
  }
  
  # Convert interface rasters to polygons
  int115_p <- terra::as.polygons(interface_zones$forest_near_ag)
  int151_p <- terra::as.polygons(interface_zones$ag_near_forest)
  
  # Union and dissolve
  int_des_p <- terra::union(int115_p, int151_p)
  int_des_d <- terra::aggregate(int_des_p)
  
  # Convert to lines
  int_des_l <- terra::as.lines(int_des_d)
  
  # Create raster template
  template <- rast(ext(boundary), resolution = resolution, crs = target_crs)
  
  # Rasterize interface lines with length calculation
  int_des_r <- terra::rasterizeGeom(int_des_l, template, fun = "length", unit = "km")
  
  # Calculate density (km per km2)
  cell_area_km2 <- (resolution/1000) * (resolution/1000)
  int_des_dens <- int_des_r / cell_area_km2
  
  return(list(density_raster = int_des_dens, interface_lines = int_des_l))
}

interface_density <- calculate_interface_density(interface_zones, lim_6372, config$analysis_resolution, config$target_crs)

# Export interface density results
if (!is.null(interface_density)) {
  writeRaster(interface_density$density_raster, file.path(data_paths$output_interface, "int_des_dens_1.tif"), overwrite = TRUE)
  writeVector(interface_density$interface_lines, file.path(data_paths$output_interface, "int_des_6372.shp"), overwrite = TRUE)
  cat("Interface density calculated and exported\n")
}

# PART 3: SATELLITE DATA PROCESSING --------------------------------------
cat("\n=== PART 3: SATELLITE DATA PROCESSING ===\n")

# Function to load satellite bands
load_satellite_bands <- function(year, sensor_type, month = "Abril") {
  cat("Loading", year, "satellite data...\n")
  if (sensor_type == "L7") {
    bands <- config$landsat7_bands
    sensor_folder <- "Landsat_7"
  } else if (sensor_type %in% c("L8", "L9")) {
    bands <- config$landsat8_bands
    sensor_folder <- ifelse(sensor_type == "L8", "Landsat_8", "Landsat_9")
  }
  base_path <- file.path(ifelse(month == "Abril", data_paths$satellite_april, data_paths$satellite_march), sensor_folder, year)
  band_rasters <- list()
  for (band in bands) {
    band_file <- file.path(base_path, paste0(band, ".tif"))
    if (file.exists(band_file)) {
      band_rasters[[band]] <- terra::rast(band_file)
    } else {
      warning("Band file not found: ", band_file)
    }
  }
  if (length(band_rasters) > 0) {
    return(do.call(c, band_rasters))
  }
  return(NULL)
}

# Function to convert DN to TOA reflectance
convert_to_toa <- function(raster_stack, sun_elevation_deg) {
  multiplier <- config$toa_params$multiplier
  additive <- config$toa_params$additive
  sun_elevation_rad <- sun_elevation_deg * pi / 180
  toa_raster <- (multiplier * raster_stack + additive) / sin(sun_elevation_rad)
  return(toa_raster)
}

# Function to process satellite data
process_satellite_year <- function(year, reference_raster = NULL) {
  metadata <- satellite_metadata[[year]]
  if (is.null(metadata)) {
    if (year %in% c("2012", "2013", "2017", "2018", "2022")) {
      sensor_type <- "L7"
      month <- ifelse(year %in% c("2017", "2018", "2022"), "Marzo", "Abril")
    } else {
      return(NULL)
    }
  } else {
    sensor_type <- metadata$sensor
    month <- "Abril"
  }
  satellite_stack <- load_satellite_bands(year, sensor_type, month)
  if (is.null(satellite_stack)) return(NULL)
  if (!is.null(metadata)) {
    satellite_stack <- convert_to_toa(satellite_stack, metadata$sun_elevation)
  }
  satellite_6372 <- terra::project(satellite_stack, config$target_crs)
  satellite_cropped <- terra::crop(satellite_6372, lim_6372)
  if (!is.null(reference_raster)) {
    satellite_cropped <- terra::resample(satellite_cropped, reference_raster, method = "bilinear")
  }
  return(satellite_cropped)
}

# Process all years
available_years <- c("2012", "2013", "2014", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")
reference_2024 <- process_satellite_year("2024")
if (is.null(reference_2024)) {
  stop("Could not load 2024 reference data")
}
satellite_data <- list()
for (year in available_years) {
  cat("Processing year:", year, "\n")
  processed_data <- process_satellite_year(year, reference_2024)
  if (!is.null(processed_data)) {
    satellite_data[[year]] <- processed_data
    output_file <- file.path(data_paths$output_stacks, paste0(year, ifelse(year %in% names(satellite_metadata), "_r.tif", ".tif")))
    terra::writeRaster(processed_data, output_file, overwrite = TRUE)
    cat("Exported stack for", year, "\n")
  }
}

# PART 4: SPECTRAL INDICES CALCULATION -----------------------------------
cat("\n=== PART 4: SPECTRAL INDICES CALCULATION ===\n")

# Function to calculate spectral indices
calculate_spectral_indices <- function(satellite_stack, year, sensor_type = "L8") {
  cat("Calculating spectral indices for", year, "\n")
  if (sensor_type == "L7") {
    band_mapping <- list(blue = "B1", green = "B2", red = "B3", nir = "B4", swir2 = "B5", swir3 = "B7")
    scale_factor <- 10000
  } else {
    band_mapping <- list(blue = "B2", green = "B3", red = "B4", nir = "B5", swir2 = "B6", swir3 = "B7")
    scale_factor <- 1
  }
  tryCatch({
    indices <- RStoolbox::spectralIndices(
      satellite_stack,
      blue = band_mapping$blue,
      green = band_mapping$green,
      red = band_mapping$red,
      nir = band_mapping$nir,
      swir2 = band_mapping$swir2,
      swir3 = band_mapping$swir3,
      scaleFactor = scale_factor,
      indices = c("NDVI", "EVI", "WDVI", "MSAVI")
    )
    new_names <- paste0(c("NDVI", "EVI", "WDVI", "MSAVI"), year)
    names(indices) <- new_names
    return(indices)
  }, error = function(e) {
    cat("Error calculating indices for", year, ":", e$message, "\n")
    return(NULL)
  })
}

# Calculate indices
indices_data <- list()
for (year in names(satellite_data)) {
  sensor_type <- ifelse(year %in% c("2012", "2013", "2017", "2018", "2022"), "L7", "L8")
  indices <- calculate_spectral_indices(satellite_data[[year]], year, sensor_type)
  if (!is.null(indices)) {
    indices_data[[year]] <- indices
    output_file <- file.path(data_paths$output_indices, paste0(year, ".tif"))
    terra::writeRaster(indices, output_file, overwrite = TRUE)
  }
}

# PART 5: MULTI-TEMPORAL ANALYSIS ----------------------------------------
cat("\n=== PART 5: MULTI-TEMPORAL ANALYSIS ===\n")

# Function to extract index time series
extract_index_timeseries <- function(indices_data, index_name) {
  index_layers <- list()
  for (year in names(indices_data)) {
    layer_name <- paste0(index_name, year)
    if (layer_name %in% names(indices_data[[year]])) {
      index_layers[[year]] <- indices_data[[year]][[layer_name]]
    }
  }
  if (length(index_layers) > 0) {
    combined <- do.call(c, index_layers)
    return(combined)
  }
  return(NULL)
}

# Extract time series
ndvi_timeseries <- extract_index_timeseries(indices_data, "NDVI")
evi_timeseries <- extract_index_timeseries(indices_data, "EVI")
wdvi_timeseries <- extract_index_timeseries(indices_data, "WDVI")
msavi_timeseries <- extract_index_timeseries(indices_data, "MSAVI")

# Clamp indices
if (!is.null(ndvi_timeseries)) ndvi_clamped <- terra::clamp(ndvi_timeseries, lower = -1, upper = 1)
if (!is.null(evi_timeseries)) evi_clamped <- terra::clamp(evi_timeseries, lower = -1, upper = 1)
if (!is.null(msavi_timeseries)) msavi_clamped <- terra::clamp(msavi_timeseries, lower = -1, upper = 1)

# PART 6: VEGETATION DENSITY ANALYSIS ------------------------------------
cat("\n=== PART 6: VEGETATION DENSITY ANALYSIS ===\n")

# Dense vegetation analysis (WDVI)
if (!is.null(wdvi_timeseries)) {
  cat("Analyzing dense vegetation using WDVI...\n")
  wdvi_mean <- terra::app(wdvi_timeseries, mean, na.rm = TRUE)
  wdvi_values <- terra::values(wdvi_mean, na.rm = TRUE)
  q75_threshold <- quantile(wdvi_values, 0.75, na.rm = TRUE)
  wdvi_q3 <- wdvi_mean > q75_threshold
  wdvi_q3_clipped <- terra::crop(wdvi_q3, lim_6372, mask = TRUE)
  agg_factor <- round(100 / terra::res(wdvi_q3_clipped)[1])
  wdvi_aggregated <- terra::aggregate(wdvi_q3_clipped, fact = agg_factor, na.rm = TRUE)
  cell_res <- terra::res(wdvi_aggregated)[1]
  cell_area_ha <- (cell_res / 100)^2
  wdvi_density <- (wdvi_aggregated * cell_area_ha) * 100
  writeRaster(wdvi_q3_clipped, file.path(data_paths$output_wdvi, "WDVI_Q3_30.tif"), overwrite = TRUE)
  writeRaster(wdvi_density, file.path(data_paths$output_wdvi, "WDVI_Q3_dens_ha.tif"), overwrite = TRUE)
  cat("Dense vegetation analysis complete\n")
}

# Canopy continuity analysis (MSAVI)
if (!is.null(msavi_clamped)) {
  cat("Analyzing canopy continuity using MSAVI...\n")
  msavi_mean <- terra::app(msavi_clamped, mean, na.rm = TRUE)
  msavi_values <- terra::values(msavi_mean, na.rm = TRUE)
  q25_threshold <- quantile(msavi_values, 0.25, na.rm = TRUE)
  msavi_q1 <- msavi_mean <= q25_threshold
  msavi_q1_clipped <- terra::crop(msavi_q1, lim_6372, mask = TRUE)
  agg_factor <- round(100 / terra::res(msavi_q1_clipped)[1])
  msavi_aggregated <- terra::aggregate(msavi_q1_clipped, fact = agg_factor, na.rm = TRUE)
  cell_res <- terra::res(msavi_aggregated)[1]
  cell_area_ha <- (cell_res / 100)^2
  msavi_density <- (msavi_aggregated * cell_area_ha) * 100
  writeRaster(msavi_q1_clipped, file.path(data_paths$output_msavi, "MSAVI_Q1_30.tif"), overwrite = TRUE)
  writeRaster(msavi_density, file.path(data_paths$output_msavi, "MSAVI_Q1_dens_ha.tif"), overwrite = TRUE)
  cat("Canopy continuity analysis complete\n")
}

# PART 7: EXPORT FINAL RESULTS AND SUMMARY -------------------------------
cat("\n=== PART 7: EXPORTING FINAL RESULTS AND SUMMARY ===\n")

# Export land cover results
land_cover_exports <- list(
  "cover_10_6372.tif" = cover_10_6372,
  "forest_all_types.tif" = forest_6372,
  "fire_prone_forest.tif" = fire_prone_6372
)
for (filename in names(land_cover_exports)) {
  output_path <- file.path(data_paths$output_cover, filename)
  tryCatch({
    terra::writeRaster(land_cover_exports[[filename]], output_path, overwrite = TRUE)
    cat("Exported:", filename, "\n")
  }, error = function(e) {
    cat("Failed to export", filename, ":", e$message, "\n")
  })
}

# Create summary report
create_summary_report <- function() {
  cat("Creating summary report...\n")
  processing_log <- data.frame(
    Component = c("Land Cover Analysis", "Interface Density", "Satellite Processing", "Spectral Indices",
                  "Dense Vegetation", "Canopy Continuity"),
    Status = c("Complete", ifelse(!is.null(interface_density), "Complete", "Skipped"),
               "Complete", "Complete",
               ifelse(!is.null(wdvi_timeseries), "Complete", "Skipped"),
               ifelse(!is.null(msavi_clamped), "Complete", "Skipped")),
    Details = c(
      paste(nrow(composition_2010), "land cover classes"),
      ifelse(!is.null(interface_density), "Deciduous-agriculture interface density", "No interface data"),
      paste(length(satellite_data), "years processed"),
      paste(length(indices_data), "years with indices"),
      ifelse(!is.null(wdvi_timeseries), "WDVI Q3 analysis", "No WDVI data"),
      ifelse(!is.null(msavi_clamped), "MSAVI Q1 analysis", "No MSAVI data")
    )
  )
  write.csv(processing_log, file.path(data_paths$output_veg, "processing_log.csv"), row.names = FALSE)
  return(processing_log)
}

# Generate summary
final_summary <- create_summary_report()

# Print completion message
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Processing date:", as.character(Sys.Date()), "\n")
cat("Years processed:", length(satellite_data), "/", length(available_years), "\n")
cat("Spectral indices calculated for:", length(indices_data), "years\n")
cat("Land cover classes identified:", nrow(composition_2010), "\n")
cat("Interface density analysis:", ifelse(!is.null(interface_density), "Complete", "Skipped"), "\n")
cat("Results exported to output directories.\n")
cat("Check processing_log.csv for detailed status.\n")