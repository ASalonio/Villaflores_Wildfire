# Villaflores Forest Fires Model Evaluation
# Author: Augusto Salonio
# Date: August 1, 2025
# Description: Script for analyzing fire occurrence prediction using logistic regression,
#              including AUC evaluation, variable importance, model performance metrics,
#              and merging predictions with spatial data.
# Dependencies: Run after model training (model.R)
# Input: Reg_Scale_Fix_F.csv, Regresores_1km_F.shp, model_F15.RData
# Output: Model discriminatory graphs and malla_mod_F15.shp (SpatVector with predictions and regressors)

# --- Load Libraries ---
library(sp)
library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(PresenceAbsence)
library(pROC)
library(ROCR)

# Note: `terra` and `dplyr` can conflict; use `terra::` when needed.

# --- Set Working Directory and Load Data ---
setwd("~/Villaflores/R/Puntos")

# Input paths
datapath_malla <- "~/Villaflores/R/output/Malla"
datapath_ocurrencia <- file.path(datapath_malla, "Ocurrencia y Densidad")  # Corrected from "Ocurrencia"

# Load data with error handling
datafile_csv <- file.path(datapath_ocurrencia, "Reg_Scale_Fix_F.csv")
if (!file.exists(datafile_csv)) stop("File not found: ", datafile_csv)
reg_scale_fix <- read_csv2(datafile_csv)

datafile_shp <- file.path(datapath_ocurrencia, "Regresores_1km_F.shp")
if (!file.exists(datafile_shp)) stop("File not found: ", datafile_shp)
malla_mod_F <- terra::vect(datafile_shp)

# Load pre-trained model
modelfile <- "model_F15.RData"
if (!file.exists(modelfile)) stop("Model file not found: ", modelfile)
load(modelfile)

# Verify required columns in reg_scale_fix
required_cols <- c("ID", "Ocurrencia", "Len_Hab", "Perc_Des", "Perc_Ridge", "Len_Car", 
                   "WDVI_Q3", "Len_Int_Des", "MSAVI_Q1", "Perc_SW", "Prec", "Pob_Dens", "Wind")
missing_cols <- setdiff(required_cols, colnames(reg_scale_fix))
if (length(missing_cols) > 0) stop("Missing columns in reg_scale_fix: ", paste(missing_cols, collapse = ", "))

# --- Model Discrimination Accuracy: AUC of ROC ---
# Initialize AUC storage
auc_values <- numeric(100)

# Perform 100 repetitions for AUC calculation
set.seed(123)  # Ensure reproducibility
for (i in 1:100) {
  # Split data (70% train, 30% test)
  train_indices <- sample(nrow(reg_scale_fix), size = 0.7 * nrow(reg_scale_fix))
  train_set <- reg_scale_fix[train_indices, ] %>% select(-ID, -any_of("geometry"))
  test_set <- reg_scale_fix[-train_indices, ] %>% select(-ID, -any_of("geometry"))
  
  # Train logistic regression model
  logit_model <- glm(
    Ocurrencia ~ Len_Hab +
      Perc_Des +
      Perc_Ridge +
      Len_Car:poly(WDVI_Q3, 2)[, 2] +
      Len_Int_Des:poly(MSAVI_Q1, 2) +
      Perc_SW:ns(Prec, df = 3)[, 1] +
      Pob_Dens:ns(Wind, df = 3) +
      ns(Wind, df = 3):poly(WDVI_Q3, 2):poly(MSAVI_Q1, 2),
    data = train_set,
    family = binomial(link = "logit"),
    control = glm.control(maxit = 50)
  )
  
  # Predict on test set
  pred_probs <- predict(logit_model, test_set, type = "response")
  
  # Handle NA predictions
  valid_idx <- !is.na(pred_probs)
  if (sum(valid_idx) == 0) {
    warning("All predictions are NA in iteration ", i, "; skipping")
    next
  }
  pred_probs <- pred_probs[valid_idx]
  actual <- test_set$Ocurrencia[valid_idx]
  
  # Calculate AUC
  pred_obj <- prediction(pred_probs, actual)
  auc_values[i] <- performance(pred_obj, "auc")@y.values[[1]]
  
  # Plot ROC curve for first iteration
  if (i == 1) {
    perf_obj <- performance(pred_obj, "tpr", "fpr")
    plot(perf_obj, main = "ROC Curve (Iteration 1)", col = "red", lwd = 2)
    abline(a = 0, b = 1, lwd = 2, lty = 3, col = "black")
  }
}

# Summarize AUC results
cat("Mean AUC over 100 repetitions:", mean(auc_values, na.rm = TRUE), "\n")
cat("Standard deviation of AUC:", sd(auc_values, na.rm = TRUE), "\n")

# Plot AUC histogram
hist(
  auc_values,
  main = "AUC Distribution",
  sub = paste("Mean AUC:", round(mean(auc_values, na.rm = TRUE), 3),
              "/ SD AUC:", round(sd(auc_values, na.rm = TRUE), 3)),
  xlab = "AUC",
  col = "skyblue",
  border = "white"
)

# --- Variable Importance Using Permutation ---
# Predict probabilities using model_F15
original_predictions <- predict(model_F15, reg_scale_fix, type = "response")
actual <- as.numeric(as.factor(reg_scale_fix$Ocurrencia)) - 1

# Calculate original AUC
original_auc <- roc(actual, original_predictions)$auc

# Define variables and interactions
variables_to_include <- setdiff(colnames(reg_scale_fix), c("Ocurrencia", "Altitud", "Euc_Pob", "geometry"))
interaction_terms <- c(
  "Pob_Dens:Wind",
  "Len_Car:WDVI_Q3",
  "Len_Int_Des:MSAVI_Q1",
  "Perc_SW:Prec",
  "Wind:WDVI_Q3:MSAVI_Q1"
)
all_variables <- c(variables_to_include, interaction_terms)

# Initialize results
variable_contribution <- data.frame(Variable = all_variables, AUC_Drop = NA)

# Permutation importance
for (var in all_variables) {
  permuted_data <- reg_scale_fix
  if (var %in% interaction_terms) {
    vars_in_interaction <- strsplit(var, ":")[[1]]
    for (int_var in vars_in_interaction) {
      permuted_data[[int_var]] <- sample(permuted_data[[int_var]])
    }
  } else {
    permuted_data[[var]] <- sample(permuted_data[[var]])
  }
  
  # Generate predictions with permuted data
  permuted_predictions <- predict(model_F15, permuted_data, type = "response")
  
  # Calculate AUC drop
  permuted_auc <- roc(actual, permuted_predictions)$auc
  variable_contribution$AUC_Drop[variable_contribution$Variable == var] <- original_auc - permuted_auc
}

# Sort and print results
variable_contribution <- variable_contribution[order(-variable_contribution$AUC_Drop), ]
print(variable_contribution)

# Plot variable importance
ggplot(variable_contribution, aes(x = reorder(Variable, AUC_Drop), y = AUC_Drop)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(x = "Variable", y = "AUC Drop") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

# --- Model Discrimination Accuracy: Sensitivity and Specificity ---
# Prepare evaluation data
pred_model_F15 <- predict(model_F15, reg_scale_fix, type = "response")
EvalData <- reg_scale_fix %>%
  mutate(
    ID = as.integer(ID),
    Ocurrencia = as.numeric(Ocurrencia),
    Pred_F15 = round(pred_model_F15, 3)
  ) %>%
  select(ID, Ocurrencia, Pred_F15)

# Calculate accuracy metrics
accuracy_F15 <- presence.absence.accuracy(EvalData, threshold = 0.1, find.auc = TRUE)
accuracy_F15[, -c(1, 2)] <- signif(accuracy_F15[, -c(1, 2)], digits = 2)
print(accuracy_F15[c("threshold", "sensitivity", "specificity")])

# Metrics for a specific threshold
meva_F15 <- ecospat.meva.table(EvalData$Pred_F15, EvalData$Ocurrencia, 0.2)
print(meva_F15)

# --- Merge Predictions with Spatial Data ---
# Validate lengths
if (nrow(malla_mod_F) != length(pred_model_F15)) {
  stop("Length mismatch between malla_mod_F and pred_model_F15")
}

# Verify required columns in malla_mod_F
required_spatial_cols <- c("ID", "Perc_Ridge", "Perc_SW", "Wind", "Prec", "Len_Int_De", 
                           "Perc_Des", "MSAVI_Q1", "WDVI_Q3", "Len_Hab", "Len_Car", "Pob_Dens", "Ocurrencia")
missing_spatial_cols <- setdiff(required_spatial_cols, names(malla_mod_F))
if (length(missing_spatial_cols) > 0) stop("Missing columns in malla_mod_F: ", paste(missing_spatial_cols, collapse = ", "))

# Add predictions to spatial data
malla_mod_F$Pred_F15 <- pred_model_F15

# Select and rename variables in spatial data
malla_mod_F15 <- malla_mod_F %>%
  tidyterra::select(
    ID,
    Perc_Ridge,
    Perc_SW,
    Wind,
    Prec,
    Len_Int_Des = Len_Int_De,  # Rename Len_Int_De to Len_Int_Des
    Perc_Des,
    MSAVI_Q1,
    WDVI_Q3,
    Len_Hab,
    Len_Car,
    Pob_Dens,
    Ocurrencia,
    Pred_F15
  )

# Save the output spatial data
output_shp <- file.path(datapath_ocurrencia, "malla_mod_F15.shp")
terra::writeVector(malla_mod_F15, output_shp, overwrite = TRUE)
cat("Output saved to:", output_shp, "\n")