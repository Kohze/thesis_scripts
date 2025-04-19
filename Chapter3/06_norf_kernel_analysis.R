# --- thesis_scripts/Chapter3/06_norf_kernel_analysis.R ---
# Purpose: Develop and analyze the "nORF Kernel", a score derived from a reduced set
#          of top PCM features. Includes PCA visualization of the kernel score
#          and classification using the reduced feature set.
# Chapter Relevance: 3.5 The nORF-Kernel

# Load necessary packages
library(data.table)
library(ggplot2)
library(xgboost)
library(caret)
library(pROC)
library(xtable) # For LaTeX table output

# Define Input/Output files
norf_features_file <- "nORFs_pcm_features.csv"
orf_features_file <- "ORFs_pcm_features.csv"
model_output_file <- "xgboost_norf_kernel_classifier.model"
importance_plot_file <- "xgboost_kernel_feature_importance.png"
roc_plot_file <- "xgboost_kernel_roc_curve.png"
kernel_plot_file <- "norf_kernel_vs_length.png"
pca_plot_file <- "norf_kernel_pca.png"

# Define the features used for the kernel (based on original script)
kernel_features <- c("prop3.G1.residue100", "VS577", "DAYM780201.lag30",
                     "BHAR880101.lag29", "CIDH920105.lag10", "BIGC670101.lag30")

# --- Load PCM Features --- 
cat("Loading PCM feature files...\n")
if (!file.exists(norf_features_file)) stop("nORF features file not found: ", norf_features_file)
if (!file.exists(orf_features_file)) stop("ORF features file not found: ", orf_features_file)

nORF <- fread(norf_features_file)
ORF <- fread(orf_features_file)
cat("Features loaded.\n")

# Clean up potential index columns
nORF[, V1 := NULL]; ORF[, V11 := NULL]
if("X" %in% names(nORF)) nORF[, X := NULL]; if("X.1" %in% names(nORF)) nORF[, X.1 := NULL]
if("X" %in% names(ORF)) ORF[, X := NULL]; if("X.1" %in% names(ORF)) ORF[, X.1 := NULL]

# Combine nORF and ORF data
corrFrame <- rbindlist(list(nORF, ORF), fill = TRUE)

# --- Calculate Sequence Lengths (if not already present) ---
# Requires loading the sequences again, or having saved lengths
# Assuming lengths need to be recalculated from loaded sequences:
loaded_data_file <- "chapter3_loaded_data.RData"
if (file.exists(loaded_data_file)) {
    load(loaded_data_file) # Loads nORFAAseq, nORFGeneIds, cORFSequencesAA, cORFSeqNames
    norf_lengths <- data.table(id = nORFGeneIds, length = nchar(nORFAAseq))
    orf_lengths <- data.table(id = cORFSeqNames, length = nchar(cORFSequencesAA))
    all_lengths <- rbindlist(list(norf_lengths, orf_lengths))
    corrFrame <- merge(corrFrame, all_lengths, by="id", all.x=TRUE)
    rm(nORFAAseq, nORFGeneIds, cORFSequencesAA, cORFSeqNames, norf_lengths, orf_lengths, all_lengths) # Clean up memory
} else {
    warning("Sequence data file not found. Cannot calculate lengths.")
    # Add placeholder length if needed, or stop if length is critical
    # corrFrame[, length := NA_integer_]
}

# --- Prepare Data for Kernel --- 
cat("Preparing data for nORF Kernel analysis...\n")

# Check if kernel features exist
missing_kernel_features <- setdiff(kernel_features, names(corrFrame))
if (length(missing_kernel_features) > 0) {
    stop("Missing required kernel features in the data: ", paste(missing_kernel_features, collapse=", "))
}

# Select relevant columns and handle NAs
kernel_data <- corrFrame[, .SD, .SDcols = c("id", "type", "length", kernel_features)]
kernel_data <- na.omit(kernel_data) # Remove rows with NAs in kernel features or length

# --- Calculate nORF Kernel Score --- 
cat("Calculating nORF Kernel score...\n")
# Simple sum of the selected features (as in original script)
kernel_data[, kernel_score := rowSums(.SD), .SDcols = kernel_features]
# Original script also calculated log(sum), which might be sensitive to negative sums
# kernel_data[, kernel_score_log := log(rowSums(.SD) + 1e-6), .SDcols = kernel_features] # Add small offset for log

# --- Visualize Kernel Score --- 
cat("Visualizing nORF Kernel score...\n")

p_kernel_vs_length <- ggplot(kernel_data, aes(x = length, y = kernel_score, color = type)) +
    geom_point(alpha = 0.3) +
    ggtitle("nORF Kernel Score vs Sequence Length") +
    xlab("Amino Acid Length") +
    ylab("nORF-Kernel Score (Sum of Selected Features)") +
    scale_color_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800"), name = "Type") +
    theme_minimal() +
    theme(legend.position = "bottom")

print(p_kernel_vs_length)
ggsave(file.path(".", kernel_plot_file), plot = p_kernel_vs_length, width = 8, height = 6)

# --- PCA on Kernel Features --- 
cat("Performing PCA on kernel features...\n")

pca_data <- kernel_data[, .SD, .SDcols = kernel_features]
pca_result <- prcomp(pca_data, scale = TRUE, center = TRUE)

# Extract first two principal components
pcaDF <- data.table(
    PCA1 = pca_result$x[, 1],
    PCA2 = pca_result$x[, 2],
    type = kernel_data$type
)

# Plot PCA results
p_pca <- ggplot(pcaDF, aes(x = PCA1, y = PCA2, color = type)) +
    geom_point(alpha = 0.3) +
    ggtitle("PCA of nORF Kernel Features") +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    theme_minimal() +
    scale_color_manual(values = c("nORF" = "#00AFBB", "ORF" = "#E7B800"), name = "Type") +
    theme(legend.position = "bottom")

print(p_pca)
ggsave(file.path(".", pca_plot_file), plot = p_pca, width = 7, height = 6)


# --- XGBoost Classification using Kernel Features --- 
cat("\nPreparing data for XGBoost classification using kernel features...\n")

# Prepare data matrix and labels
data_matrix_kernel <- as.matrix(kernel_data[, .SD, .SDcols = kernel_features])
labels_kernel <- as.numeric(ifelse(kernel_data$type == "nORF", 1, 0))

# Check data
if (nrow(data_matrix_kernel) == 0) stop("No data remaining for kernel classification.")

# Split Data (using same seed and proportions as in script 04)
set.seed(123)
train_ind_k <- createDataPartition(labels_kernel, p = 0.70, list = FALSE)
train_data_k <- data_matrix_kernel[train_ind_k, ]
train_labels_k <- labels_kernel[train_ind_k]
temp_data_k <- data_matrix_kernel[-train_ind_k, ]
temp_labels_k <- labels_kernel[-train_ind_k]

test_ind_k <- createDataPartition(temp_labels_k, p = 0.50, list = FALSE)
test_data_k <- temp_data_k[test_ind_k, ]
test_labels_k <- temp_labels_k[test_ind_k]
validation_data_k <- temp_data_k[-test_ind_k, ]
validation_labels_k <- temp_labels_k[-test_ind_k]

# Convert to DMatrix
dtrain_k <- xgb.DMatrix(data = train_data_k, label = train_labels_k)
dtest_k <- xgb.DMatrix(data = test_data_k, label = test_labels_k)
dvalid_k <- xgb.DMatrix(data = validation_data_k, label = validation_labels_k)

cat("Training XGBoost model using kernel features...\n")

# Use similar parameters as before, or retune if necessary
params_kernel <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = 0.1,
    max_depth = 5, # Might need less depth with fewer features
    gamma = 1,
    alpha = 1,
    lambda = 1
)
nrounds_kernel <- 100 # May converge faster
watchlist_k <- list(train = dtrain_k, validation = dvalid_k)

bst_kernel <- xgb.train(
    params = params_kernel,
    data = dtrain_k,
    nrounds = nrounds_kernel,
    watchlist = watchlist_k,
    early_stopping_rounds = 15,
    print_every_n = 10,
    verbose = 1
)

cat("Kernel feature model training complete. Best iteration:", bst_kernel$best_iteration, "\n")

# --- Evaluate Kernel Feature Model --- 
cat("Evaluating kernel feature model on test set...\n")

preds_test_prob_k <- predict(bst_kernel, dtest_k)
preds_test_class_k <- as.integer(preds_test_prob_k > 0.5)

# Confusion Matrix
conf_mat_test_k <- confusionMatrix(data = factor(preds_test_class_k, levels=c(0,1)),
                                   reference = factor(test_labels_k, levels=c(0,1)))
cat("Kernel Features - Test Set Confusion Matrix:\n")
print(conf_mat_test_k)

# Accuracy
accuracy_test_k <- conf_mat_test_k$overall['Accuracy']
cat("Kernel Features - Test Set Accuracy:", sprintf("%.4f", accuracy_test_k), "\n")

# ROC Curve
roc_obj_k <- roc(response = test_labels_k, predictor = preds_test_prob_k, levels=c(0, 1))
auc_test_k <- auc(roc_obj_k)
cat("Kernel Features - Test Set AUC:", sprintf("%.4f", auc_test_k), "\n")

# Plot ROC Curve
png(file.path(".", roc_plot_file), width=600, height=600)
plot(roc_obj_k, main = paste("XGBoost ROC Curve (Kernel Features) - AUC:", sprintf("%.3f", auc_test_k)),
     print.auc = TRUE, legacy.axes = TRUE, lwd = 2)
dev.off()
cat("Saved kernel ROC plot to:", roc_plot_file, "\n")

# --- Feature Importance (Kernel Model) --- 
cat("Calculating kernel feature importance...\n")
importance_matrix_k <- xgb.importance(model = bst_kernel)
print(importance_matrix_k)

# Save importance as LaTeX table (as in original script)
latex_table <- xtable(importance_matrix_k,
                      caption = "XGBoost Feature Importance (Kernel Features Model)",
                      align = c("c","c","c","c","c","c"))
print(latex_table, type = "latex", file = file.path(".", "xgboost_kernel_importance.tex"),
      include.rownames = FALSE, floating = FALSE, hline.after = c(-1, 0))
cat("Saved kernel feature importance LaTeX table.\n")

# Plot importance
png(file.path(".", importance_plot_file), width=600, height=400)
tryCatch(xgb.plot.importance(importance_matrix = importance_matrix_k, top_n = length(kernel_features)),
         error=function(e) { plot.new(); title("Error plotting kernel importance") })
dev.off()
cat("Saved kernel feature importance plot to:", importance_plot_file, "\n")

# --- Save the Kernel Model --- 
cat("Saving trained kernel XGBoost model to:", model_output_file, "\n")
xgb.save(bst_kernel, file.path(".", model_output_file))

cat("--- 06_norf_kernel_analysis.R finished ---\n") 