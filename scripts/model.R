# Villaflores Fire Occurrence Models
# Author: Augusto Salonio
# Date: August 1, 2025
# Description: This script performs logistic regression modeling to predict fire occurrence
#              based on topographic, climatic, demographic, and vegetation variables.
#              Includes univariate analysis, multivariate modeling, and model selection.
# Dependencies: Multiple packages for spatial analysis, modeling, and visualization
# Usage: Run after completing grid-based spatial analysis
#        Input: Regresores_seleccionados_1km_F.csv from previous analysis
#        Output: Model objects, visualizations, and processed datasets

# Load Required Libraries ------------------------------------------------
library(sp)
library(terra)
library(sf)
library(readr)
library(tidyr)
library(dplyr) # Note: Use terra:: explicitly when using with dplyr
library(stringr)
library(tibble)
library(openxlsx)
library(lubridate)
library(parzer)
library(leaflet)
library(ggplot2)
library(car)
library(splines)
library(mgcv)
library(lme4)
library(lmtest)
library(glmmTMB)
library(performance)
library(pscl)
library(DescTools)
library(pROC)
library(fmsb)
library(cvAUC)
library(MuMIn)

# Configuration ----------------------------------------------------------
# Set working directory and paths
base_dir <- getwd()
data_paths <- list(
  input = file.path(base_dir, "output/Malla/Ocurrencia y Densidad"),
  output_models = file.path(base_dir, "output/Models"),
  output_plots = file.path(base_dir, "output/Models/Plots")
)

# Create output directories
for (dir in unlist(data_paths[grepl("output", names(data_paths))])) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Configuration parameters
config <- list(
  alpha_level = 0.05,
  delta_aic = 4,           # Threshold for model selection
  spline_df = 3,           # Degrees of freedom for splines
  poly_degree = 2          # Degree for polynomial terms
)

# Load and Prepare Data --------------------------------------------------
cat("Loading and preparing data...\n")

# Load preprocessed data
datafile <- file.path(data_paths$input, "Regresores_seleccionados_1km_F.csv")
if (!file.exists(datafile)) {
  stop("Input data file not found: ", datafile)
}

reg_raw <- read_csv2(datafile)
cat("Loaded", nrow(reg_raw), "observations with", ncol(reg_raw), "variables\n")

# Display data structure
glimpse(reg_raw)

# Data Preprocessing and Standardization ---------------------------------
cat("Preprocessing and standardizing data...\n")

prepare_modeling_data <- function(data) {
  # Separate occurrence variable
  ocurrencia <- data %>%
    dplyr::select(Ocurrencia)
  
  # Identify and standardize numeric predictors
  id_cols <- c("ID", "Ocurrencia", "Densidad")
  numeric_predictors <- data %>%
    dplyr::select(!any_of(id_cols)) %>%
    dplyr::select(where(is.numeric))
  
  # Standardize numeric variables
  standardized_predictors <- scale(numeric_predictors) %>%
    as.data.frame()
  
  # Combine standardized predictors with occurrence
  modeling_data <- bind_cols(standardized_predictors, ocurrencia)
  
  # Add ID back
  if ("ID" %in% names(data)) {
    modeling_data$ID <- as.factor(data$ID)
  }
  
  # Convert occurrence to factor
  modeling_data$Ocurrencia <- as.factor(modeling_data$Ocurrencia)
  
  return(modeling_data)
}

# Prepare modeling dataset
reg_scale_fix <- prepare_modeling_data(reg_raw)

# Check for missing values
missing_summary <- colSums(is.na(reg_scale_fix))
cat("Missing values summary:\n")
print(missing_summary[missing_summary > 0])

# Display final dataset structure
cat("Final modeling dataset structure:\n")
glimpse(reg_scale_fix)

# Null Model --------------------------------------------------------------
cat("\n=== NULL MODEL ===\n")

# Fit null model for comparison
mod_null <- glm(Ocurrencia ~ 1, data = reg_scale_fix, family = binomial)
cat("Null model summary:\n")
print(summary(mod_null))
cat("Null model AIC:", AIC(mod_null), "\n")
cat("Null model deviance:", deviance(mod_null), "on", df.residual(mod_null), "degrees of freedom\n")

# Univariate Analysis -----------------------------------------------------
cat("\n=== UNIVARIATE ANALYSIS ===\n")

perform_univariate_analysis <- function(data, predictors, output_path) {
  univariate_results <- list()
  
  for (predictor in predictors) {
    if (predictor %in% names(data)) {
      cat("Analyzing:", predictor, "\n")
      
      # Fit univariate model
      formula_str <- paste("Ocurrencia ~", predictor)
      model <- glm(as.formula(formula_str), family = binomial, data = data)
      
      # Store results
      univariate_results[[predictor]] <- list(
        model = model,
        aic = AIC(model),
        deviance = deviance(model),
        p_value = summary(model)$coefficients[2, 4] # P-value of the predictor
      )
      
      # Create prediction data for plotting
      predictor_range <- seq(min(data[[predictor]], na.rm = TRUE), 
                             max(data[[predictor]], na.rm = TRUE), 
                             length.out = 100)
      prediction_data <- data.frame(x = predictor_range)
      names(prediction_data)[1] <- predictor
      prediction_data$predicted_prob <- predict(model, newdata = prediction_data, type = "response")
      
      # Create and save plot
      p <- ggplot() +
        geom_point(data = data, aes_string(x = predictor, y = "as.numeric(Ocurrencia) - 1"), 
                   alpha = 0.5, color = "gray60") +
        geom_line(data = prediction_data, aes_string(x = predictor, y = "predicted_prob"), 
                  color = "blue", size = 1.2) +
        labs(title = paste("Logistic Regression:", predictor),
             x = predictor,
             y = "Fire Occurrence Probability") +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      # Save plot
      plot_filename <- file.path(output_path, paste0("univariate_", predictor, ".png"))
      ggsave(plot_filename, p, width = 10, height = 6, dpi = 300)
    }
  }
  
  return(univariate_results)
}

# Define predictors for analysis
predictors <- c("Len_Hab", "Perc_Des", "Perc_Ridge", "Altitud", "Euc_Pob")

# Perform univariate analysis
univariate_results <- perform_univariate_analysis(reg_scale_fix, predictors, data_paths$output_plots)

# Summary of univariate results
cat("Univariate analysis summary:\n")
univariate_summary <- data.frame(
  Variable = names(univariate_results),
  AIC = sapply(univariate_results, function(x) x$aic),
  Deviance = sapply(univariate_results, function(x) x$deviance),
  P_Value = sapply(univariate_results, function(x) x$p_value)
)
print(univariate_summary)

# Export univariate results
write_csv(univariate_summary, file.path(data_paths$output_models, "univariate_results.csv"))

# Interaction Effects Analysis -------------------------------------------
cat("\n=== INTERACTION EFFECTS ANALYSIS ===\n")

create_interaction_plots <- function(data, interactions, output_path) {
  for (i in seq_along(interactions)) {
    interaction <- interactions[[i]]
    cat("Creating interaction plot for:", interaction$formula, "\n")
    
    # Fit model
    model <- glm(as.formula(paste("Ocurrencia ~", interaction$formula)), 
                 family = binomial, data = data)
    
    # Create prediction grid
    vars <- interaction$vars
    grid_data <- expand.grid(
      x = seq(min(data[[vars[1]]], na.rm = TRUE), max(data[[vars[1]]], na.rm = TRUE), length.out = 50),
      y = seq(min(data[[vars[2]]], na.rm = TRUE), max(data[[vars[2]]], na.rm = TRUE), length.out = 50)
    )
    names(grid_data) <- vars
    
    # Add additional variables if needed for complex interactions
    if (length(vars) > 2) {
      for (var in vars[3:length(vars)]) {
        grid_data[[var]] <- mean(data[[var]], na.rm = TRUE)
      }
    }
    
    # Predict probabilities
    grid_data$predicted_prob <- predict(model, newdata = grid_data, type = "response")
    
    # Create heatmap
    p <- ggplot(grid_data, aes_string(x = vars[1], y = vars[2], fill = "predicted_prob")) +
      geom_tile() +
      scale_fill_gradient(low = "blue", high = "red", name = "Probability") +
      labs(title = interaction$title,
           x = interaction$x_label,
           y = interaction$y_label) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    # Save plot
    plot_filename <- file.path(output_path, paste0("interaction_", names(interactions)[i], ".png"))
    ggsave(plot_filename, p, width = 10, height = 8, dpi = 300)
  }
}

# Define interaction effects to analyze
interactions <- list(
  car_wdvi = list(
    formula = "Len_Car:poly(WDVI_Q3, 2)[,2]",
    vars = c("Len_Car", "WDVI_Q3"),
    title = "Interaction: Road Density × Forest Greenness",
    x_label = "Road Density",
    y_label = "Forest Greenness (WDVI_Q3)"
  ),
  int_msavi = list(
    formula = "Len_Int_Des:poly(MSAVI_Q1, 2)",
    vars = c("Len_Int_Des", "MSAVI_Q1"),
    title = "Interaction: Forest-Agriculture Interface × Herbaceous Layer",
    x_label = "Forest-Agriculture Interface Density",
    y_label = "Herbaceous Layer Density (MSAVI_Q1)"
  ),
  sw_prec = list(
    formula = "Perc_SW:ns(Prec, df = 3)[,1]",
    vars = c("Perc_SW", "Prec"),
    title = "Interaction: Southwest Slope × October Precipitation",
    x_label = "Southwest Slope Percentage",
    y_label = "October Precipitation"
  ),
  pop_wind = list(
    formula = "Pob_Dens:ns(Wind, df = 3)",
    vars = c("Pob_Dens", "Wind"),
    title = "Interaction: Population Density × Wind Speed",
    x_label = "Population Density",
    y_label = "Wind Speed"
  )
)

# Create interaction plots
create_interaction_plots(reg_scale_fix, interactions, data_paths$output_plots)

# Multivariate Model Development -----------------------------------------
cat("\n=== MULTIVARIATE MODEL DEVELOPMENT ===\n")

# Define full model with all interactions
full_model_formula <- "Ocurrencia ~ 
  Altitud +
  Euc_Pob +
  Perc_Des +
  Perc_Ridge +
  Perc_SW:ns(Prec, df = 3)[, 1] +
  Pob_Dens*ns(Wind, df = 3) +
  Len_Car:poly(WDVI_Q3, 2)[, 2] +
  Len_Hab +
  Len_Int_Des*poly(MSAVI_Q1, 2) +
  ns(Wind, df = 3):poly(WDVI_Q3, 2)*poly(MSAVI_Q1, 2)"

# Fit full model
cat("Fitting full multivariate model...\n")
full_model <- glm(as.formula(full_model_formula), 
                  data = reg_scale_fix, 
                  family = binomial)

# Model diagnostics
cat("Full model summary:\n")
print(summary(full_model))

# Model performance metrics
cat("Model performance:\n")
print(model_performance(full_model))

# Pseudo R-squared values
cat("Pseudo R-squared values:\n")
print(PseudoR2(full_model, which = "all"))

# Check multicollinearity
cat("Collinearity check:\n")
print(check_collinearity(full_model))

# Model Selection Using Information Criteria -----------------------------
cat("\n=== MODEL SELECTION ===\n")

# Set up for model dredging
options(na.action = "na.fail")

cat("Performing automated model selection (this may take a while)...\n")
# Model dredging to find best combinations
model_dredge <- dredge(full_model, rank = "AIC")

# Select models within delta AIC threshold
model_subset <- subset(model_dredge, delta < config$delta_aic)

cat("Number of models within delta AIC <", config$delta_aic, ":", nrow(model_subset), "\n")

# Show top models
cat("Top models:\n")
print(model_subset)

# Model averaging
cat("Performing model averaging...\n")
model_averaged <- model.avg(model_subset)

cat("Model averaging summary:\n")
print(summary(model_averaged))

# Final Selected Model ---------------------------------------------------
cat("\n=== FINAL MODEL ===\n")

# Define final model based on model selection results
final_model_formula <- "Ocurrencia ~ 
  Len_Hab +
  Perc_Des +
  Perc_Ridge +
  Len_Car:poly(WDVI_Q3, 2)[, 2] +
  Len_Int_Des*poly(MSAVI_Q1, 2) +
  Perc_SW:ns(Prec, df = 3)[, 1] +
  Pob_Dens*ns(Wind, df = 3) +
  ns(Wind, df = 3):poly(WDVI_Q3, 2)*poly(MSAVI_Q1, 2)"

# Fit final model
final_model <- glm(as.formula(final_model_formula),
                   data = reg_scale_fix,
                   family = binomial(link = "logit"))

# Final model diagnostics
cat("Final model summary:\n")
print(summary(final_model))

cat("Final model performance:\n")
print(model_performance(final_model))

cat("Final model pseudo R-squared:\n")
print(PseudoR2(final_model, which = "all"))

cat("Final model collinearity check:\n")
print(check_collinearity(final_model))

# Model validation plots
cat("Creating model diagnostic plots...\n")
png(file.path(data_paths$output_plots, "model_diagnostics.png"), width = 1200, height = 800)
check_model(final_model)
dev.off()

# Export Results ----------------------------------------------------------
cat("\n=== EXPORTING RESULTS ===\n")

# Save model objects
save(final_model, full_model, model_averaged, 
     file = file.path(data_paths$output_models, "fire_occurrence_models.RData"))

# Save processed dataset
write_csv2(reg_scale_fix, 
           file.path(data_paths$input, "Reg_Scale_Fix_F.csv"))

# Create model comparison summary
model_comparison <- data.frame(
  Model = c("Null", "Full", "Final"),
  AIC = c(AIC(mod_null), AIC(full_model), AIC(final_model)),
  Deviance = c(deviance(mod_null), deviance(full_model), deviance(final_model)),
  DF = c(df.residual(mod_null), df.residual(full_model), df.residual(final_model)),
  Pseudo_R2_McFadden = c(0, PseudoR2(full_model, "McFadden"), PseudoR2(final_model, "McFadden"))
)

write_csv(model_comparison, file.path(data_paths$output_models, "model_comparison.csv"))

# Create summary report
create_model_report <- function(final_model, model_comparison, output_path) {
  report <- list(
    processing_date = as.character(Sys.Date()),
    n_observations = nrow(reg_scale_fix),
    n_fire_events = sum(as.numeric(reg_scale_fix$Ocurrencia) - 1),
    fire_rate = mean(as.numeric(reg_scale_fix$Ocurrencia) - 1),
    final_model_aic = AIC(final_model),
    final_model_pseudo_r2 = PseudoR2(final_model, "McFadden"),
    n_predictors = length(coefficients(final_model)) - 1,
    significant_predictors = sum(summary(final_model)$coefficients[-1, 4] < 0.05)
  )
  
  # Convert to dataframe for easy export
  report_df <- data.frame(
    Metric = names(report),
    Value = unlist(report)
  )
  
  write_csv(report_df, file.path(output_path, "modeling_report.csv"))
  return(report_df)
}

# Generate final report
final_report <- create_model_report(final_model, model_comparison, data_paths$output_models)

# Print completion summary
cat("\n=== MODELING COMPLETE ===\n")
cat("Processing date:", as.character(Sys.Date()), "\n")
cat("Total observations:", nrow(reg_scale_fix), "\n")
cat("Fire events:", sum(as.numeric(reg_scale_fix$Ocurrencia) - 1), "\n")
cat("Fire rate:", round(mean(as.numeric(reg_scale_fix$Ocurrencia) - 1), 3), "\n")
cat("Final model AIC:", round(AIC(final_model), 2), "\n")
cat("Final model pseudo R²:", round(PseudoR2(final_model, "McFadden"), 3), "\n")
cat("Results exported to:", data_paths$output_models, "\n")
cat("Check modeling_report.csv for detailed summary.\n")