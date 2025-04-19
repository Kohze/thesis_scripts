# --- thesis_scripts/Chapter3/04_xgboost_classification.R ---
# Purpose: Train and evaluate an XGBoost model to classify sequences
#          (nORF vs ORF) based on their PCM features. Includes feature importance
#          analysis and hyperparameter tuning using caret.
# Chapter Relevance: 3.4 High-Dimensional Statistical Analysis

# Load necessary packages
library(data.table)
library(xgboost)
library(caret)
library(pROC)
library(ggplot2)

# Define Input/Output files (Ensure these match output from 02_...R)
norf_features_file <- "nORFs_pcm_features.csv"
orf_features_file <- "ORFs_pcm_features.csv"
# Note: Random features are not typically used for classification training/testing

output_dir <- "."
model_output_file <- "xgboost_norf_orf_classifier.model"
importance_plot_file <- "xgboost_feature_importance.png"
roc_plot_file <- "xgboost_roc_curve.png"
tuning_plot_file <- "xgboost_caret_tuning.png"

# --- Load Processed PCM Features --- 
cat("Loading PCM feature files...\n")
if (!file.exists(norf_features_file)) stop("nORF features file not found: ", norf_features_file)
if (!file.exists(orf_features_file)) stop("ORF features file not found: ", orf_features_file)

nORF <- fread(norf_features_file)
ORF <- fread(orf_features_file)
cat("Features loaded.\n")

# Clean up potential index columns
nORF[, V1 := NULL]
ORF[, V11 := NULL]
if("X" %in% names(nORF)) nORF[, X := NULL]
if("X.1" %in% names(nORF)) nORF[, X.1 := NULL]
if("X" %in% names(ORF)) ORF[, X := NULL]
if("X.1" %in% names(ORF)) ORF[, X.1 := NULL]

# --- Prepare Data for XGBoost --- 
cat("Preparing data for XGBoost classification...\n")

# Combine nORF and ORF data
dataset <- rbindlist(list(nORF, ORF), fill = TRUE)

# Remove rows with NAs introduced during feature calculation or joining
cat("Initial rows:", nrow(dataset), "\n")
dataset <- na.omit(dataset)
cat("Rows after na.omit:", nrow(dataset), "\n")

# Convert type to numeric label (nORF=1, ORF=0)
dataset[, label := ifelse(type == "nORF", 1, 0)]

# Select feature columns (exclude id, type, original label)
feature_cols <- setdiff(names(dataset), c("id", "type", "label"))

# Optional: Check for and remove zero-variance columns again after combining/cleaning
feature_data_matrix <- as.matrix(dataset[, ..feature_cols])
variances <- apply(feature_data_matrix, 2, var, na.rm = TRUE)
zero_var_cols <- names(which(variances < 1e-10))
if (length(zero_var_cols) > 0) {
    cat("Removing zero/low variance columns:", paste(zero_var_cols, collapse=", "), "\n")
    feature_cols <- setdiff(feature_cols, zero_var_cols)
}

# Prepare final data matrix and labels
data_matrix <- as.matrix(dataset[, ..feature_cols])
labels <- as.numeric(dataset$label)

# Check if data exists
if (nrow(data_matrix) == 0 || length(labels) == 0) {
    stop("No data remaining after preparation for XGBoost.")
}

cat("Data prepared. Features:", ncol(data_matrix), "Samples:", nrow(data_matrix), "\n")

# --- Split Data (Train/Test/Validation) --- 
set.seed(123) # for reproducibility

# First split: 70% train, 30% temp (test + validation)
train_ind <- createDataPartition(labels, p = 0.70, list = FALSE)
train_data <- data_matrix[train_ind, ]
train_labels <- labels[train_ind]
temp_data <- data_matrix[-train_ind, ]
temp_labels <- labels[-train_ind]

# Second split: 50% test, 50% validation from temp
test_ind <- createDataPartition(temp_labels, p = 0.50, list = FALSE)
test_data <- temp_data[test_ind, ]
test_labels <- temp_labels[test_ind]
validation_data <- temp_data[-test_ind, ]
validation_labels <- temp_labels[-test_ind]

cat("Data split complete:\n")
cat(" Train:", nrow(train_data), "samples\n")
cat(" Test:", nrow(test_data), "samples\n")
cat(" Validation:", nrow(validation_data), "samples\n")

# Convert to DMatrix format for XGBoost
dtrain <- xgb.DMatrix(data = train_data, label = train_labels)
dtest <- xgb.DMatrix(data = test_data, label = test_labels)
dvalid <- xgb.DMatrix(data = validation_data, label = validation_labels)

# --- Initial XGBoost Model Training --- 
cat("\nTraining initial XGBoost model...\n")

# Define parameters (taken from original script, adjust as needed)
# Using 'binary:logistic' for binary classification, 'auc' for evaluation
params_initial <- list(
    objective = "binary:logistic",
    eval_metric = "auc", # Area Under Curve - good for binary classification
    eta = 0.1,         # Learning rate
    max_depth = 10,      # Max tree depth
    gamma = 7.5,         # Minimum loss reduction for split
    alpha = 2,           # L1 regularization
    lambda = 0.4         # L2 regularization
    # subsample = 0.8,   # Subsample ratio of the training instance
    # colsample_bytree = 0.8 # Subsample ratio of columns when constructing each tree
)

nrounds_initial <- 250
watchlist <- list(train = dtrain, validation = dvalid)

bst_initial <- xgb.train(
    params = params_initial,
    data = dtrain,
    nrounds = nrounds_initial,
    watchlist = watchlist,
    early_stopping_rounds = 20, # Stop if validation AUC doesn't improve for 20 rounds
    print_every_n = 25,
    verbose = 1
)

cat("Initial model training complete. Best iteration:", bst_initial$best_iteration, "\n")

# --- Evaluate Initial Model --- 
cat("Evaluating initial model on test set...\n")

preds_test_prob <- predict(bst_initial, dtest)
preds_test_class <- as.integer(preds_test_prob > 0.5) # Threshold at 0.5

# Confusion Matrix
conf_mat_test <- confusionMatrix(data = factor(preds_test_class, levels=c(0,1)),
                                 reference = factor(test_labels, levels=c(0,1)))
cat("Test Set Confusion Matrix:\n")
print(conf_mat_test)

# Accuracy
accuracy_test <- conf_mat_test$overall['Accuracy']
cat("Test Set Accuracy:", sprintf("%.4f", accuracy_test), "\n")

# ROC Curve
roc_obj <- roc(response = test_labels,
               predictor = preds_test_prob,
               levels=c(0, 1))
auc_test <- auc(roc_obj)
cat("Test Set AUC:", sprintf("%.4f", auc_test), "\n")

# Plot ROC Curve
png(file.path(output_dir, roc_plot_file), width=600, height=600)
plot(roc_obj, main = paste("XGBoost ROC Curve (Test Set) - AUC:", sprintf("%.3f", auc_test)),
     print.auc = TRUE, legacy.axes = TRUE, lwd = 2)
dev.off()
cat("Saved ROC plot to:", roc_plot_file, "\n")


# --- Feature Importance --- 
cat("Calculating and plotting feature importance...\n")
importance_matrix <- xgb.importance(model = bst_initial)
print(head(importance_matrix, 20)) # Print top 20 features

# Plot importance
png(file.path(output_dir, importance_plot_file), width=800, height=1000)
tryCatch(xgb.plot.importance(importance_matrix = importance_matrix, top_n = 50),
         error=function(e) { plot.new(); title("Error plotting importance") })
dev.off()
cat("Saved feature importance plot to:", importance_plot_file, "\n")

# --- Optional: Hyperparameter Tuning with Caret --- 
run_tuning <- FALSE # Set to TRUE to run the potentially time-consuming tuning

if (run_tuning) {
    cat("\nStarting hyperparameter tuning with caret (this may take a while)...\n")

    # Define grid (smaller example than original for speed)
    tune_grid <- expand.grid(
      nrounds = c(100, 250, 500),
      eta = c(0.05, 0.1, 0.2),
      max_depth = c(5, 8, 10),
      gamma = c(1, 5, 10),
      min_child_weight = 1,
      colsample_bytree = 0.8,
      subsample = 0.8
    )

    # Define control parameters
    tune_control <- caret::trainControl(
      method = "cv",             # Cross-validation
      number = 3,                # Number of folds
      verboseIter = TRUE,       # Show progress
      allowParallel = TRUE,
      classProbs = TRUE,        # Needed for AUC summary
      summaryFunction = twoClassSummary # Use AUC for selection
    )

    # Run tuning (Ensure labels are factors for caret classification)
    train_labels_factor <- factor(ifelse(train_labels == 1, "nORF", "ORF"), levels=c("ORF", "nORF"))

    xgb_tune <- tryCatch(caret::train(
      x = train_data,
      y = train_labels_factor,
      trControl = tune_control,
      tuneGrid = tune_grid,
      method = "xgbTree",
      metric = "ROC", # Select based on AUC
      verbose = FALSE
    ), error = function(e) { cat("Error during caret tuning:", conditionMessage(e), "\n"); NULL })

    if (!is.null(xgb_tune)) {
        cat("Caret tuning complete. Best parameters:\n")
        print(xgb_tune$bestTune)

        # Plot tuning results
        png(file.path(output_dir, tuning_plot_file), width=800, height=600)
        print(ggplot(xgb_tune) + theme_bw() + labs(title="XGBoost Hyperparameter Tuning Results (caret)"))
        dev.off()
        cat("Saved tuning plot to:", tuning_plot_file, "\n")

        # Optionally, train final model with best parameters found by caret
        # final_params <- c(xgb_tune$bestTune, objective = "binary:logistic", eval_metric = "auc")
        # final_nrounds <- xgb_tune$bestTune$nrounds
        # bst_final <- xgb.train(params=final_params, data=dtrain, nrounds=final_nrounds, watchlist=watchlist, early_stopping_rounds=20)
        # ... re-evaluate bst_final ...
    }
} else {
    cat("\nSkipping hyperparameter tuning.\n")
}

# --- Save the Initial Model --- 
cat("Saving initial trained XGBoost model to:", model_output_file, "\n")
xgb.save(bst_initial, file.path(output_dir, model_output_file))

cat("--- 04_xgboost_classification.R finished ---\n") 