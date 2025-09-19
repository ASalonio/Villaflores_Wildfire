# Villaflores Climate Analysis
# Author: Augusto Salonio
# Date: 14th Jan 2025
# Description: Downloads and processes WorldClim climate data including wind speed,
#              solar radiation, maximum temperature, and precipitation for the 
#              Villaflores region in Chiapas, Mexico

# Load required libraries ------------------------------------------------
library(terra)
library(geodata)
library(sf)
library(dplyr)

# Configuration -----------------------------------------------------------
# Set base directory (modify as needed)
base_dir <- "~/Villaflores/R"

# Define data paths
data_paths <- list(
  limits = file.path(base_dir, "data/Límites"),
  climate_output = file.path(base_dir, "output/Climatico")
)

# Create output directory if it doesn't exist
if (!dir.exists(data_paths$climate_output)) {
  dir.create(data_paths$climate_output, recursive = TRUE)
}

# WorldClim parameters
worldclim_config <- list(
  country = "MEX",
  version = "2.1",
  resolution = 0.5,  # 30 arcsec (~1km at equator)
  variables = c("wind", "srad", "tmax", "prec")
)

# Seasonal selection (month numbers)
seasonal_config <- list(
  dry_season_month = 4,   # April (end of dry season)
  wet_season_month = 10   # October (end of wet season)
)

# Load Boundary Data ------------------------------------------------------
cat("Loading boundary data...\n")

# Load Villaflores boundary
villaflores_file <- file.path(data_paths$limits, "Villaflores.shp")
if (!file.exists(villaflores_file)) {
  # Try alternative location
  villaflores_file <- file.path(data_paths$limits, "Viejos", "Villaflores.shp")
}

if (!file.exists(villaflores_file)) {
  stop("Villaflores boundary file not found. Please check the file path.")
}

lim_villa <- terra::vect(villaflores_file)

# Check and display boundary info
cat("Boundary CRS:", terra::crs(lim_villa), "\n")
cat("Boundary extent:", as.vector(terra::ext(lim_villa)), "\n")

# Download Climate Data ---------------------------------------------------
cat("Downloading WorldClim climate data for Mexico...\n")

# Function to safely download WorldClim data
download_worldclim <- function(var, country, path, version, res) {
  cat("Downloading", var, "data...\n")
  
  tryCatch({
    data <- worldclim_country(
      country = country,
      path = path,
      version = version,
      res = res,
      var = var
    )
    
    cat("Successfully downloaded", var, "data\n")
    return(data)
    
  }, error = function(e) {
    cat("Error downloading", var, ":", e$message, "\n")
    return(NULL)
  })
}

# Download all climate variables
climate_data <- list()

for (var in worldclim_config$variables) {
  climate_data[[var]] <- download_worldclim(
    var = var,
    country = worldclim_config$country,
    path = data_paths$climate_output,
    version = worldclim_config$version,
    res = worldclim_config$resolution
  )
}

# Check if all downloads were successful
successful_downloads <- sapply(climate_data, function(x) !is.null(x))
if (!all(successful_downloads)) {
  failed_vars <- names(climate_data)[!successful_downloads]
  warning("Failed to download: ", paste(failed_vars, collapse = ", "))
}

# Process Climate Data ---------------------------------------------------
cat("Processing climate data...\n")

# Function to process each climate variable
process_climate_var <- function(climate_raster, var_name, boundary) {
  if (is.null(climate_raster)) {
    cat("Skipping", var_name, "- data not available\n")
    return(NULL)
  }
  
  cat("Processing", var_name, "...\n")
  
  # Crop to study area
  cropped <- terra::crop(climate_raster, boundary, mask = TRUE)
  
  # Reproject to EPSG:6372
  reprojected <- terra::project(cropped, "EPSG:6372")
  
  return(reprojected)
}

# Process all climate variables
processed_climate <- list()

for (var in names(climate_data)) {
  processed_climate[[var]] <- process_climate_var(
    climate_data[[var]], 
    var, 
    lim_villa
  )
}

# Extract Seasonal Data --------------------------------------------------
cat("Extracting seasonal data...\n")

# Function to extract specific months
extract_seasonal_data <- function(climate_stack, month, var_name) {
  if (is.null(climate_stack)) {
    cat("Cannot extract", var_name, "- data not available\n")
    return(NULL)
  }
  
  if (terra::nlyr(climate_stack) < month) {
    warning("Month ", month, " not available for ", var_name)
    return(NULL)
  }
  
  # Extract specific month
  monthly_data <- climate_stack[[month]]
  
  # Create informative name
  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  names(monthly_data) <- paste0(var_name, "_", month_names[month])
  
  return(monthly_data)
}

# Extract seasonal climate data
seasonal_climate <- list()

# Dry season variables (April)
for (var in c("wind", "srad", "tmax")) {
  if (!is.null(processed_climate[[var]])) {
    seasonal_climate[[paste0(var, "_apr")]] <- extract_seasonal_data(
      processed_climate[[var]], 
      seasonal_config$dry_season_month, 
      var
    )
  }
}

# Wet season variable (October precipitation)
if (!is.null(processed_climate[["prec"]])) {
  seasonal_climate[["prec_oct"]] <- extract_seasonal_data(
    processed_climate[["prec"]], 
    seasonal_config$wet_season_month, 
    "prec"
  )
}

# Data Quality Checks ----------------------------------------------------
cat("Performing data quality checks...\n")

# Function to check data quality
check_data_quality <- function(raster_data, var_name) {
  if (is.null(raster_data)) {
    return(data.frame(
      Variable = var_name,
      Status = "Missing",
      Min = NA,
      Max = NA,
      Mean = NA,
      NA_cells = NA,
      Total_cells = NA
    ))
  }
  
  values_vec <- terra::values(raster_data, na.rm = FALSE)
  
  return(data.frame(
    Variable = var_name,
    Status = "Available",
    Min = round(min(values_vec, na.rm = TRUE), 2),
    Max = round(max(values_vec, na.rm = TRUE), 2),
    Mean = round(mean(values_vec, na.rm = TRUE), 2),
    NA_cells = sum(is.na(values_vec)),
    Total_cells = length(values_vec)
  ))
}

# Create quality report
quality_report <- do.call(rbind, lapply(names(seasonal_climate), function(var) {
  check_data_quality(seasonal_climate[[var]], var)
}))

print(quality_report)

# Export Results ----------------------------------------------------------
cat("Exporting climate data...\n")

# Define output filenames with proper extensions
output_files <- list(
  wind_apr = "wind_apr_6372.tif",
  srad_apr = "srad_apr_6372.tif", 
  tmax_apr = "tmax_apr_6372.tif",
  prec_oct = "prec_oct_6372.tif"
)

# Export seasonal climate data
export_summary <- data.frame(
  Variable = character(),
  Filename = character(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (var in names(output_files)) {
  filename <- output_files[[var]]
  filepath <- file.path(data_paths$climate_output, filename)
  
  if (!is.null(seasonal_climate[[var]])) {
    tryCatch({
      terra::writeRaster(seasonal_climate[[var]], filepath, overwrite = TRUE)
      export_summary <- rbind(export_summary, data.frame(
        Variable = var,
        Filename = filename,
        Status = "Success"
      ))
      cat("Exported:", filename, "\n")
    }, error = function(e) {
      export_summary <- rbind(export_summary, data.frame(
        Variable = var,
        Filename = filename,
        Status = paste("Error:", e$message)
      ))
      cat("Failed to export:", filename, "-", e$message, "\n")
    })
  } else {
    export_summary <- rbind(export_summary, data.frame(
      Variable = var,
      Filename = filename,
      Status = "No data available"
    ))
    cat("Skipped:", filename, "- no data available\n")
  }
}

# Create Metadata File ---------------------------------------------------
cat("Creating metadata file...\n")

metadata <- list(
  project_info = list(
    title = "Climate Data Processing for Villaflores Region",
    date_processed = Sys.Date(),
    description = "WorldClim climate data processed for fire risk analysis"
  ),
  data_source = list(
    source = "WorldClim",
    version = worldclim_config$version,
    resolution = paste(worldclim_config$resolution, "degrees"),
    country = worldclim_config$country,
    url = "https://www.worldclim.org/"
  ),
  processing_info = list(
    target_crs = "EPSG:6372",
    dry_season_month = seasonal_config$dry_season_month,
    wet_season_month = seasonal_config$wet_season_month,
    boundary_file = basename(villaflores_file)
  ),
  variables = list(
    wind = "Wind speed (m/s) - April",
    srad = "Solar radiation (kJ/m²/day) - April", 
    tmax = "Maximum temperature (°C) - April",
    prec = "Precipitation (mm) - October"
  )
)

# Save metadata as JSON
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(metadata, 
                       file.path(data_paths$climate_output, "metadata.json"), 
                       pretty = TRUE)
}

# Save reports as CSV
write.csv(quality_report, 
          file.path(data_paths$climate_output, "data_quality_report.csv"), 
          row.names = FALSE)

write.csv(export_summary, 
          file.path(data_paths$climate_output, "export_summary.csv"), 
          row.names = FALSE)

# Final Summary -----------------------------------------------------------
cat("\n=== CLIMATE ANALYSIS COMPLETE ===\n")
cat("Output directory:", data_paths$climate_output, "\n")
cat("Files exported:", sum(export_summary$Status == "Success"), "/", nrow(export_summary), "\n")
cat("Check the quality report and export summary for details.\n")

# Print final summary
if (nrow(export_summary) > 0) {
  cat("\nExport Summary:\n")
  print(export_summary)
}