# Villaflores Grid-Based Spatial Analysis for Fire Risk
# Author: Augusto Salonio
# Date: August 1, 2025
# Description: This script creates a 1km² grid and a 30m mold for spatial analysis in Villaflores,
#              processes topographic, climatic, demographic, and vegetation data,
#              performs statistical analysis (normality, correlation, VIF),
#              and incorporates fire occurrence data for fire risk assessment.
# Dependencies: sp, terra, sf, readr, tidyr, dplyr, stringr, tibble, openxlsx,
#               lubridate, parzer, leaflet, corrplot, usdm, tidyterra
# Usage: Place input data in the 'data' directory and run the script from the base directory.
#        Outputs are saved in the 'output' directory.

# Load Required Libraries ------------------------------------------------
library(sp)
library(terra)
library(sf)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(openxlsx)
library(lubridate)
library(parzer)
library(leaflet)
library(corrplot)
library(usdm)
library(tidyterra) # For terra-tidyverse integration

# Note: Use terra:: explicitly with dplyr to avoid function conflicts

# Configuration ----------------------------------------------------------
# Set base directory (modify as needed)
base_dir <- getwd()

# Define data paths
data_paths <- list(
  limits = file.path(base_dir, "data/Límites"),
  topographic = file.path(base_dir, "data/Topografico"),
  climatic = file.path(base_dir, "data/Climatico"),
  human = file.path(base_dir, "data/Humano"),
  vegetation = file.path(base_dir, "data/Vegetación"),
  ignition = file.path(base_dir, "data/Puntos de Ignición - Forestales"),
  output_malla = file.path(base_dir, "output/Malla"),
  output_plots = file.path(base_dir, "output/Plots"),
  output_ignition = file.path(base_dir, "output/Puntos de Ignición - Forestales"),
  output_regressors = file.path(base_dir, "output/Malla/Ocurrencia y Densidad")
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
  target_crs = "EPSG:6372",
  grid_resolution = 1000, # 1km² grid
  mold_resolution = 30,   # 30m mold
  k_neighbors = 8         # For spatial imputation
)

# Load Boundary and Input Data -------------------------------------------
cat("Loading boundary and input data...\n")

# Load boundary shapefile
tryCatch({
  lim_6372 <- terra::vect(file.path(data_paths$limits, "lim_6372.shp"))
}, error = function(e) {
  stop("Failed to load boundary data: ", e$message)
})

# Load DEM
tryCatch({
  dem_6372 <- terra::rast(file.path(data_paths$topographic, "dem_6372.tif"))
}, error = function(e) {
  stop("DEM file not found: ", e$message)
})

# Load fire ignition points
ignition_file <- file.path(data_paths$ignition, "Seleccion_2012-2024.shp")
if (file.exists(ignition_file)) {
  ig_6372 <- terra::vect(ignition_file)
  cat("Loaded", nrow(ig_6372), "ignition points\n")
} else {
  ig_6372 <- NULL
  warning("Ignition points file not found: ", ignition_file)
}

# PART 1: CREATE 1KM² GRID -----------------------------------------------
cat("\n=== PART 1: CREATE 1KM² GRID ===\n")

create_1km_grid <- function(dem, lim, config, output_path) {
  # Calculate aggregation factor for 1km² grid
  agg_factor <- config$grid_resolution / terra::res(dem)[1]
  
  # Aggregate to 1km²
  malla_1 <- terra::aggregate(dem, fact = agg_factor, na.rm = TRUE)
  values(malla_1) <- NA
  
  # Assign cell IDs
  cell_ids <- 1:terra::ncell(malla_1)
  values(malla_1) <- cell_ids
  names(malla_1) <- "ID"
  
  # Clip to boundary
  clip_malla_1 <- terra::crop(malla_1, lim, mask = TRUE)
  
  # Convert to vector
  clip_malla_1_v <- terra::as.polygons(clip_malla_1)
  
  # Export results
  terra::writeRaster(clip_malla_1, file.path(output_path, "malla.tif"), overwrite = TRUE)
  terra::writeVector(clip_malla_1_v, file.path(output_path, "malla.shp"), overwrite = TRUE)
  
  cat("1km² grid created and exported\n")
  return(list(raster = clip_malla_1, vector = clip_malla_1_v))
}

# Create grid
malla_result <- create_1km_grid(dem_6372, lim_6372, config, data_paths$output_malla)

# PART 2: CREATE 30M MOLD ------------------------------------------------
cat("\n=== PART 2: CREATE 30M MOLD ===\n")

create_30m_mold <- function(dem, lim, config, output_path) {
  # Calculate aggregation factor for 30m resolution
  agg_factor <- round(config$mold_resolution / terra::res(dem)[1])
  
  # Aggregate to 30m
  molde_30 <- terra::aggregate(dem, fact = agg_factor, na.rm = TRUE)
  values(molde_30) <- NA
  
  # Assign cell IDs
  cell_ids <- 1:terra::ncell(molde_30)
  values(molde_30) <- cell_ids
  names(molde_30) <- "ID"
  
  # Clip to boundary
  clip_molde_30 <- terra::crop(molde_30, lim, mask = TRUE)
  
  # Export result
  terra::writeRaster(clip_molde_30, file.path(output_path, "molde_30.tif"), overwrite = TRUE)
  
  cat("30m mold created and exported\n")
  return(clip_molde_30)
}

# Create mold
molde_30 <- create_30m_mold(dem_6372, lim_6372, config, data_paths$output_malla)

# PART 3: PROCESS EXPLANATORY VARIABLES -----------------------------------
cat("\n=== PART 3: PROCESS EXPLANATORY VARIABLES ===\n")

process_explanatory_vars <- function(molde, config, data_paths, output_path) {
  # Define input rasters and their names
  raster_configs <- list(
    list(file = file.path(data_paths$topographic, "dem_6372_30.tif"), name = "Altitud"),
    list(file = file.path(data_paths$topographic, "perc_sw_dens_ha.tif"), name = "Perc_SW"),
    list(file = file.path(data_paths$topographic, "perc_ridge_dens_ha.tif"), name = "Perc_Ridge"),
    list(file = file.path(data_paths$climatic, "wind_abr_6372.tif"), name = "Wind"),
    list(file = file.path(data_paths$climatic, "tmax_abr_6372.tif"), name = "Tmax"),
    list(file = file.path(data_paths$climatic, "prec_oct_6372.tif"), name = "Prec"),
    list(file = file.path(data_paths$human, "Censos/OK/OK_dens_1km.tif"), name = "Pob_Dens"),
    list(file = file.path(data_paths$human, "Censos/Distancia/euc_pob_30.tif"), name = "Euc_Pob"),
    list(file = file.path(data_paths$human, "Caminos/car_dens_1.tif"), name = "Len_Car"),
    list(file = file.path(data_paths$human, "Caminos/hab_dens_1.tif"), name = "Len_Hab"),
    list(file = file.path(data_paths$vegetation, "coberturas/percentage_des_dens_ha.tif"), name = "Perc_Des"),
    list(file = file.path(data_paths$vegetation, "Interfaz Deciduo - Agricola/int_des_dens_1.tif"), name = "Len_Int_Des"),
    list(file = file.path(data_paths$vegetation, "WDVI sobre Percentil 75/WDVI_Q3_dens_ha.tif"), name = "WDVI_Q3"),
    list(file = file.path(data_paths$vegetation, "MSAVI Percentil 10/MSAVI_Q1_dens_ha.tif"), name = "MSAVI_Q1")
  )
  
  # Load and process rasters
  raster_list <- lapply(raster_configs, function(rc) {
    if (!file.exists(rc$file)) {
      warning("File not found: ", rc$file)
      return(NULL)
    }
    rast <- terra::rast(rc$file)
    rast <- terra::resample(rast, molde, method = "bilinear")
    names(rast) <- rc$name
    return(rast)
  })
  
  # Remove NULL entries
  raster_list <- Filter(Negate(is.null), raster_list)
  
  # Combine into a stack
  s_media <- do.call(c, raster_list)
  
  # Export stack
  terra::writeRaster(s_media, file.path(output_path, "s_media_30_QII.tif"), overwrite = TRUE)
  
  cat("Explanatory variables processed and stack exported\n")
  return(s_media)
}

# Process variables
s_media <- process_explanatory_vars(molde_30, config, data_paths, data_paths$output_malla)

# PART 4: EXTRACT AND IMPUTE DATA ----------------------------------------
cat("\n=== PART 4: EXTRACT AND IMPUTE DATA ===\n")

spatial_neighbor_imputation <- function(raster_stack, malla_vector, k_neighbors) {
  # Convert vector to sf
  malla_sf <- if (inherits(malla_vector, "SpatVector")) sf::st_as_sf(malla_vector) else malla_vector
  
  # Extract mean values
  data_df <- terra::zonal(raster_stack, malla_vector, fun = "mean", na.rm = TRUE)
  
  # Impute missing values
  for (col in names(data_df)) {
    if (is.numeric(data_df[[col]]) && any(is.na(data_df[[col]]))) {
      na_indices <- which(is.na(data_df[[col]]))
      for (i in na_indices) {
        current_geom <- sf::st_centroid(malla_sf[i, ])
        distances <- sf::st_distance(current_geom, sf::st_centroid(malla_sf))
        valid_neighbors <- which(!is.na(data_df[[col]]) & 1:nrow(data_df) != i)
        
        if (length(valid_neighbors) >= k_neighbors) {
          neighbor_distances <- distances[valid_neighbors]
          k_nearest <- valid_neighbors[order(neighbor_distances)[1:k_neighbors]]
          neighbor_values <- data_df[[col]][k_nearest]
          weights <- 1 / (as.numeric(distances[k_nearest]) + 1)
          data_df[[col]][i] <- sum(neighbor_values * weights) / sum(weights)
        } else {
          data_df[[col]][i] <- median(data_df[[col]], na.rm = TRUE)
        }
      }
    }
  }
  
  return(data_df)
}

# Extract and impute
malla_1_media <- terra::zonal(s_media, malla_result$vector, fun = "mean", na.rm = TRUE)
malla_1_media_fixed <- spatial_neighbor_imputation(s_media, malla_result$vector, config$k_neighbors)

# Check for remaining NAs
cat("Remaining NAs in imputed data:", sum(is.na(malla_1_media_fixed)), "\n")

# PART 5: ANALYZE PREDICTORS ---------------------------------------------
cat("\n=== PART 5: ANALYZE PREDICTORS ===\n")

analyze_predictors <- function(data_df, output_path) {
  # Enhanced Normality Test
  cat("Performing enhanced normality testing...\n")
  
  # Initialize results dataframe
  normality_results <- data.frame(
    Variable = character(), 
    p_value = numeric(), 
    Is_Normal = logical(),
    stringsAsFactors = FALSE
  )
  
  # Perform Shapiro-Wilk test for each numeric variable
  for (col in colnames(data_df)) {
    if (is.numeric(data_df[[col]])) {  # Check if column is numeric
      # Handle sample size limitations for Shapiro-Wilk test
      sample_data <- data_df[[col]][!is.na(data_df[[col]])]
      sample_size <- length(sample_data)
      
      if (sample_size >= 3 && sample_size <= 5000) {
        test <- shapiro.test(sample_data)  # Perform Shapiro-Wilk test
        normality_results <- rbind(normality_results, 
                                   data.frame(Variable = col, 
                                              p_value = test$p.value, 
                                              Is_Normal = test$p.value > 0.05))  # p > 0.05 -> Normal
        cat(sprintf("Variable: %s, p-value: %.6f, Normal: %s\n", 
                    col, test$p.value, ifelse(test$p.value > 0.05, "Yes", "No")))
      } else {
        # For samples outside Shapiro-Wilk range
        normality_results <- rbind(normality_results, 
                                   data.frame(Variable = col, 
                                              p_value = NA, 
                                              Is_Normal = NA))
        cat(sprintf("Variable: %s - Sample size (%d) outside Shapiro-Wilk test range\n", col, sample_size))
      }
    }
  }
  
  # Export enhanced normality results
  openxlsx::write.xlsx(normality_results, file.path(output_path, "enhanced_normality_test_F11.xlsx"), rowNames = FALSE)
  write.csv(normality_results, file.path(output_path, "enhanced_normality_test_F11.csv"), row.names = FALSE)
  
  # Plot histograms and Q-Q plots
  cat("Creating histograms and Q-Q plots...\n")
  for (col in colnames(data_df)) {
    if (is.numeric(data_df[[col]])) {
      # Histogram
      png(file.path(output_path, paste0("hist_", col, ".png")), width = 800, height = 600)
      hist(data_df[[col]], main = paste("Histogram of", col), xlab = col, 
           col = "lightblue", border = "black", breaks = 30)
      dev.off()
      
      # Q-Q plot
      png(file.path(output_path, paste0("qq_", col, ".png")), width = 800, height = 600)
      qqnorm(data_df[[col]], main = paste("Q-Q Plot of", col))
      qqline(data_df[[col]], col = "red", lwd = 2)
      dev.off()
    }
  }
  
  # Correlation analysis
  cat("Performing correlation analysis...\n")
  numeric_data <- data_df[sapply(data_df, is.numeric)]
  spearman_cor <- cor(numeric_data, method = "spearman", use = "pairwise.complete.obs")
  
  png(file.path(output_path, "corrplot.png"), width = 1000, height = 1000)
  corrplot(spearman_cor, type = "lower", method = "color", 
           tl.cex = 0.8, tl.col = "black", tl.srt = 45)
  dev.off()
  
  # Export correlation matrix
  write.csv(spearman_cor, file.path(output_path, "correlation_matrix.csv"))
  
  # VIF analysis
  cat("Performing VIF analysis...\n")
  tryCatch({
    vif_results <- usdm::vifcor(numeric_data, th = 0.7, method = "spearman")
    
    # Export VIF results if available
    if (!is.null(vif_results@results)) {
      write.csv(vif_results@results, file.path(output_path, "vif_results.csv"))
    }
  }, error = function(e) {
    cat("VIF analysis failed:", e$message, "\n")
    vif_results <- NULL
  })
  
  cat("Enhanced predictor analysis completed\n")
  cat("Results saved:\n")
  cat("  - Enhanced normality test: enhanced_normality_test_F11.xlsx and .csv\n")
  cat("  - Histograms and Q-Q plots: hist_*.png and qq_*.png\n")
  cat("  - Correlation plot: corrplot.png\n")
  cat("  - Correlation matrix: correlation_matrix.csv\n")
  
  return(list(normality = normality_results, correlation = spearman_cor, vif = vif_results))
}

# Analyze predictors
analysis_results <- analyze_predictors(malla_1_media_fixed, data_paths$output_plots)

# Incorporate approved predictors
approved_vars <- c("Altitud", "Perc_Ridge", "Perc_SW", "Wind", "Prec", "Perc_Des",
                   "Len_Int_Des", "WDVI_Q3", "MSAVI_Q1", "Pob_Dens", "Euc_Pob",
                   "Len_Car", "Len_Hab")
for (var in approved_vars) {
  if (var %in% colnames(malla_1_media_fixed)) {
    malla_result$vector[[var]] <- malla_1_media_fixed[[var]]
  }
}

# PART 6: INCORPORATE FIRE OCCURRENCE ------------------------------------
cat("\n=== PART 6: INCORPORATE FIRE OCCURRENCE ===\n")

incorporate_occurrence <- function(ignition, malla_raster, malla_vector, output_path) {
  if (is.null(ignition)) {
    warning("No ignition data provided")
    return(NULL)
  }
  
  # Rasterize ignition points
  ig_rast <- terra::rasterize(ignition, malla_raster, fun = "first", background = 0)
  names(ig_rast) <- "Ocurrencia"
  
  # Extract occurrence
  ig_df <- terra::zonal(ig_rast, malla_raster, fun = "first", na.rm = TRUE)
  
  # Export occurrence
  write_csv2(ig_df, file.path(output_path, "Ocurrencia_1km.csv"))
  
  # Add to vector
  malla_vector$Ocurrencia <- ig_df$Ocurrencia
  
  cat("Fire occurrence incorporated\n")
  return(list(vector = malla_vector, df = ig_df))
}

# Incorporate occurrence
occurrence_result <- incorporate_occurrence(ig_6372, malla_result$raster, malla_result$vector, data_paths$output_ignition)

# PART 7: EXPORT FINAL RESULTS AND SUMMARY -------------------------------
cat("\n=== PART 7: EXPORT FINAL RESULTS AND SUMMARY ===\n")

# Export final results
malla_1_df <- malla_1_media_fixed
if (!is.null(occurrence_result)) {
  malla_1_df$Ocurrencia <- occurrence_result$df$Ocurrencia
  terra::writeVector(occurrence_result$vector, 
                     file.path(data_paths$output_regressors, "Regresores_1km_A3.shp"), 
                     overwrite = TRUE)
}
write_csv2(malla_1_df, file.path(data_paths$output_regressors, "Regresores_1km_A3.csv"))

# Create summary report
create_summary_report <- function(malla_df, analysis_results, output_path) {
  # Count normal variables
  normal_vars <- sum(analysis_results$normality$Is_Normal, na.rm = TRUE)
  total_vars <- nrow(analysis_results$normality)
  
  processing_log <- data.frame(
    Component = c("Grid Creation", "Mold Creation", "Explanatory Variables", 
                  "Data Extraction", "Predictor Analysis", "Fire Occurrence"),
    Status = c("Complete", "Complete", "Complete", "Complete", "Complete", 
               ifelse(!is.null(occurrence_result), "Complete", "Skipped")),
    Details = c(
      "1km² grid created and clipped",
      "30m mold created and clipped",
      paste(length(names(malla_df)) - ifelse("Ocurrencia" %in% names(malla_df), 1, 0), "explanatory variables processed"),
      paste(nrow(malla_df), "grid cells extracted"),
      paste("Enhanced normality testing:", normal_vars, "of", total_vars, "variables are normal (p > 0.05)"),
      ifelse(!is.null(occurrence_result), "Fire occurrence data incorporated", "No ignition data")
    )
  )
  
  write.csv(processing_log, file.path(output_path, "processing_log.csv"), row.names = FALSE)
  return(processing_log)
}

# Generate summary
final_summary <- create_summary_report(malla_1_df, analysis_results, data_paths$output_malla)

# Print completion message
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Processing date:", as.character(Sys.Date()), "\n")
cat("Grid cells processed:", nrow(malla_1_df), "\n")
cat("Variables analyzed:", length(approved_vars), "\n")
if (!is.null(analysis_results$normality)) {
  normal_count <- sum(analysis_results$normality$Is_Normal, na.rm = TRUE)
  total_count <- nrow(analysis_results$normality)
  cat("Normality test results:", normal_count, "of", total_count, "variables are normally distributed\n")
}
cat("Results exported to:", data_paths$output_regressors, "\n")
cat("Check processing_log.csv for detailed status.\n")