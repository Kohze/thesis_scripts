# --- thesis_scripts/Chapter3/02_pcm_feature_extraction.R ---
# Purpose: Define functions for calculating proteochemometric (PCM) features
#          and apply them to the nORF, canonical ORF, and random sequences loaded
#          in the previous step. Saves the calculated features.
# Chapter Relevance: 3.3 Proteochemometric Analysis

# Load necessary packages
library(protr)
library(stringi)
library(data.table) # Or library(dplyr) if preferred for data frame manipulation

# Define Input/Output files (Ensure these match output from 01_...R)
loaded_data_file <- "chapter3_loaded_data.RData"
output_dir <- "."
norf_features_file <- "nORFs_pcm_features.csv"
orf_features_file <- "ORFs_pcm_features.csv"
random_features_file <- "Random_pcm_features.csv"

# --- Load Data --- 
cat("Loading sequence data from:", loaded_data_file, "\n")
if (!file.exists(loaded_data_file)) stop("Loaded data file not found: ", loaded_data_file)
load(loaded_data_file)
cat("Data loaded successfully.\n")

# --- PCM Feature Extraction Functions --- 

# Wrapper functions with error handling for protr::extract* functions
brotoFunc <- function (x) {
  return(tryCatch(protr::extractMoreauBroto(x), error = function(e) NULL))
}

triadsFunc <- function (x) {
  return(tryCatch(protr::extractCTriad(x), error = function(e) NULL))
}

fAACFunc <- function (x) {
  return(tryCatch(protr::extractAAC(x), error = function(e) NULL))
}

fDCFunc <- function (x) {
  return(tryCatch(protr::extractDC(x), error = function(e) NULL))
}

fCTDCFunc <- function (x) {
  return(tryCatch(protr::extractCTDC(x), error = function(e) NULL))
}

fCTDTFunc <- function (x) {
  return(tryCatch(protr::extractCTDT(x), error = function(e) NULL))
}

fCTDDFunc <- function (x) {
  return(tryCatch(protr::extractCTDD(x), error = function(e) NULL))
}

fSOCNFunc <- function (x) {
  return(tryCatch(protr::extractSOCN(x), error = function(e) NULL))
}

fQSOFunc <- function (x) {
  return(tryCatch(protr::extractQSO(x), error = function(e) NULL))
}

fPAACFunc <- function (x) {
  return(tryCatch(protr::extractPAAC(x), error = function(e) NULL))
}

fAPAACFunc <- function (x) {
  return(tryCatch(protr::extractAPAAC(x), error = function(e) NULL))
}

# Function to convert list of feature vectors to a data frame
# Handles cases where some feature extractions might return NULL
conversionFrame <- function(feature_list) {
  # Find the first non-NULL element to get column names
  first_valid_index <- Position(Negate(is.null), feature_list)
  if (is.na(first_valid_index)) {
      warning("All feature extractions resulted in NULL.")
      return(data.frame()) # Return empty data frame
  }
  col_names <- names(feature_list[[first_valid_index]])
  if (is.null(col_names)) {
      warning("Feature names could not be determined.")
      # Attempt to create generic names if unnamed list/vector
      num_cols <- length(feature_list[[first_valid_index]])
      col_names <- paste0("V", seq_len(num_cols))
  }

  # Convert list to matrix, handling NULLs by filling with NA
  # Ensure all elements processed are vectors of the same expected length
  expected_length <- length(col_names)
  feature_matrix <- t(sapply(feature_list, function(item) {
    if (is.null(item) || length(item) != expected_length) {
      rep(NA, expected_length)
    } else {
      unlist(item)
    }
  }))

  df <- as.data.frame(feature_matrix)
  colnames(df) <- col_names
  return(df)
}

# --- Main PCM Analysis Function --- 
pchemAnalysis <- function(sequences, sequence_ids, type_label) {

  cat("\nStarting PCM analysis for:", type_label, "(", length(sequences), "sequences )\n")

  # Apply each feature extraction function
  # Consider using BiocParallel::bplapply for parallel processing if beneficial
  cat(" Extracting MoreauBroto...\n")
  mBroto <- lapply(sequences, brotoFunc)
  cat(" Extracting CTriad...\n")
  mTriads <- lapply(sequences, triadsFunc)
  cat(" Extracting AAC...\n")
  mAAC <- lapply(sequences, fAACFunc)
  cat(" Extracting DC...\n")
  mDC <- lapply(sequences, fDCFunc)
  cat(" Extracting CTDC...\n")
  mCTDC <- lapply(sequences, fCTDCFunc)
  cat(" Extracting CTDT...\n")
  mCTDT <- lapply(sequences, fCTDTFunc)
  cat(" Extracting CTDD...\n")
  mCTDD <- lapply(sequences, fCTDDFunc)
  cat(" Extracting SOCN...\n")
  mSOCN <- lapply(sequences, fSOCNFunc)
  cat(" Extracting QSO...\n")
  mQSO <- lapply(sequences, fQSOFunc)
  cat(" Extracting PAAC...\n")
  mPAAC <- lapply(sequences, fPAACFunc)
  cat(" Extracting APAAC...\n")
  mAPAAC <- lapply(sequences, fAPAACFunc)

  # Convert results to data frames
  cat("Converting results to data frames...\n")
  Fmbroto <- conversionFrame(mBroto)
  Fmtriads <- conversionFrame(mTriads)
  FmAAC <- conversionFrame(mAAC)
  FmDC <- conversionFrame(mDC)
  FmCTDC <- conversionFrame(mCTDC)
  FmCTDT <- conversionFrame(mCTDT)
  FmCTDD <- conversionFrame(mCTDD)
  FmSOCN <- conversionFrame(mSOCN)
  FmQSO <- conversionFrame(mQSO)
  FmPAAC <- conversionFrame(mPAAC)
  FmAPAAC <- conversionFrame(mAPAAC)

  # Aggregate features into a single data frame
  cat("Aggregating features...\n")
  # Use cbind only if data frames are non-empty and have correct number of rows
  # Create a base data frame with ID and type
  base_df <- data.frame(id = sequence_ids, type = type_label, stringsAsFactors = FALSE)

  # List of feature data frames to bind
  feature_dfs <- list(Fmbroto, Fmtriads, FmAAC, FmDC, FmCTDC, FmCTDT, FmCTDD, FmSOCN, FmQSO, FmPAAC, FmAPAAC)

  # Check dimensions and combine safely
  combined_df <- base_df
  for (df in feature_dfs) {
      if (nrow(df) == nrow(base_df) && ncol(df) > 0) {
          combined_df <- cbind(combined_df, df)
      } else {
          warning("Skipping a feature set due to incorrect dimensions or empty data frame.")
          # Optional: Add NA columns as placeholders if needed
      }
  }

  cat("Finished PCM analysis for:", type_label, "\n")
  return(combined_df)
}

# --- Calculate Features for Each Sequence Type ---

# nORFs
nORF_features <- pchemAnalysis(nORFAAseq, nORFGeneIds, "nORF")
cat("Writing nORF features to:", norf_features_file, "\n")
fwrite(nORF_features, file.path(output_dir, norf_features_file))

# Canonical ORFs
# Note: Canonical ORFs might have issues with certain extractors if sequences are too short/long or contain unusual characters
# Filter out potentially problematic sequences first if needed
valid_corf_indices <- which(nchar(cORFSequencesAA) > 10) # Example: minimum length filter
cORF_features <- pchemAnalysis(cORFSequencesAA[valid_corf_indices],
                               cORFSeqNames[valid_corf_indices],
                               "ORF")
# Remove rows with many NAs resulting from failed extractions
cORF_features <- na.omit(cORF_features) # Or use a threshold for NA count
cat("Writing canonical ORF features to:", orf_features_file, "\n")
fwrite(cORF_features, file.path(output_dir, orf_features_file))

# Random Sequences
random_features <- pchemAnalysis(randomSequencesAA, paste0("Random_", 1:length(randomSequencesAA)), "Random")
cat("Writing random sequence features to:", random_features_file, "\n")
fwrite(random_features, file.path(output_dir, random_features_file))

# --- Save Session Info (Optional) ---
# save(nORF_features, cORF_features, random_features, file="chapter3_pcm_features.RData")

cat("--- 02_pcm_feature_extraction.R finished ---\n") 