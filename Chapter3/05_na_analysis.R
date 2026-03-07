# --- 0. Load Necessary Libraries ---
# Ensure these are installed: install.packages(c("dplyr", "stringr"))
library(dplyr)
library(stringr)
# library(tidyr) # Not strictly needed if not pivoting for LaTeX prep

# --- 1. Data Loading and Initial Setup ---
# IMPORTANT: This script assumes the following R objects are ALREADY LOADED
#            and correctly structured in your R environment:
#
#   - merged_features_only_df: Data frame.
#     Rows are observations (e.g., sequences), columns are calculated features.
#     Expected to have 233,020 rows if it's the full nORF set before NA filtering for this analysis.
#     Column names MUST be in a format like "METHOD.FEATURE_X" (e.g., "mAAC.AAC_1")
#     to allow extraction of the feature calculation method.
#
#   - b: Data frame.
#     Must have the same number of rows as `merged_features_only_df` and be perfectly row-aligned.
#     Must contain a column named 'group' (e.g., `b$group`) which holds the
#     biological category for each observation (e.g., "norfs", "proteincoding").

# --- 2. Pre-flight Checks for Actual Data ---
message("--- Running Pre-flight Checks ---")
if (!exists("merged_features_only_df") || !is.data.frame(merged_features_only_df) || nrow(merged_features_only_df) == 0) {
  stop("ERROR: merged_features_only_df is not loaded, not a data frame, or is empty. Please load your data.")
}
if (!exists("b") || !is.data.frame(b) || nrow(b) == 0 || !("group" %in% colnames(b))) {
  stop("ERROR: b is not loaded, not a data frame, is empty, or does not contain a 'group' column. Please load your data.")
}
if (nrow(merged_features_only_df) != nrow(b)) {
  stop(paste0("ERROR: Row counts of merged_features_only_df (", nrow(merged_features_only_df),
              ") and b (", nrow(b), ") do not match. Ensure they are aligned."))
}
if (is.null(colnames(merged_features_only_df)) || length(colnames(merged_features_only_df)) == 0) {
    stop("ERROR: merged_features_only_df has no column names.")
}
message(paste("Dimensions of merged_features_only_df:",
            nrow(merged_features_only_df), "rows,",
            ncol(merged_features_only_df), "cols"))
message(paste("Dimensions of b:",
            nrow(b), "rows,",
            ncol(b), "cols"))
# Ensure b$group is treated as a factor for consistent counting
if (!is.factor(b$group)) {
  b$group <- factor(b$group)
  message("Converted b$group to a factor.")
}
message(paste("Unique levels in b$group:", paste(levels(b$group), collapse=", ")))


# --- 3. NA Row Drop Analysis & Biological Group Impact (for Table 1 summary) ---
message("\n--- Analyzing NA Row Drop Statistics & Impact on Biological Groups ---")

original_row_count_analysis <- nrow(merged_features_only_df)
complete_rows_logical_analysis <- complete.cases(merged_features_only_df)
rows_without_na_count_analysis <- sum(complete_rows_logical_analysis)
rows_dropped_count_analysis <- original_row_count_analysis - rows_without_na_count_analysis
percentage_rows_dropped_analysis <- (rows_dropped_count_analysis / original_row_count_analysis) * 100

message(paste("Original number of rows for analysis:", original_row_count_analysis))
message(paste("Number of rows with NO NAs (complete cases):", rows_without_na_count_analysis))
message(paste("Number of rows that would be dropped:", rows_dropped_count_analysis))
message(paste("Percentage of rows that would be dropped:", round(percentage_rows_dropped_analysis, 2), "%"))

# Distribution of groups in the ORIGINAL dataset
original_distribution_analysis <- b %>%
  dplyr::count(group, name = "Count_Original", .drop = FALSE) %>%
  mutate(Percentage_Original = (Count_Original / sum(Count_Original)) * 100)

# Distribution of groups AFTER removing rows with NAs
after_na_removal_distribution_analysis <- b[complete_rows_logical_analysis, , drop = FALSE] %>%
  dplyr::count(group, name = "Count_After_NA_Removal", .drop = FALSE) %>%
  mutate(Percentage_After_NA_Removal = (Count_After_NA_Removal / sum(Count_After_NA_Removal)) * 100)

# Combine distributions
summary_table_bio_groups <- full_join(original_distribution_analysis,
                                      after_na_removal_distribution_analysis,
                                      by = "group")
# Replace NAs with 0 if a group disappears entirely (e.g., if a group had 100% NAs and was small)
summary_table_bio_groups <- summary_table_bio_groups %>%
  mutate_at(vars(starts_with("Count_"), starts_with("Percentage_")), ~ifelse(is.na(.), 0, .))


# Calculate additional metrics for Table 1
summary_table_bio_groups <- summary_table_bio_groups %>%
  mutate(
    Proportion_Retained_within_Group = ifelse(Count_Original > 0, Count_After_NA_Removal / Count_Original, 0),
    Percentage_Point_Change = Percentage_After_NA_Removal - Percentage_Original
  ) %>%
  # Select and order columns for the final table structure
  select(group,
         Count_Original, Percentage_Original,
         Count_After_NA_Removal, Percentage_After_NA_Removal,
         Proportion_Retained_within_Group,
         Percentage_Point_Change) %>%
  arrange(desc(Count_Original))

# Add TOTAL row for Table 1
total_row_bio_groups <- data.frame(
  group = "TOTAL",
  Count_Original = sum(summary_table_bio_groups$Count_Original),
  Percentage_Original = 100.00,
  Count_After_NA_Removal = sum(summary_table_bio_groups$Count_After_NA_Removal),
  Percentage_After_NA_Removal = if (sum(summary_table_bio_groups$Count_After_NA_Removal) > 0) 100.00 else 0, # Handle case of 0 rows after removal
  Proportion_Retained_within_Group = if (sum(summary_table_bio_groups$Count_Original) > 0) {
                                       sum(summary_table_bio_groups$Count_After_NA_Removal) / sum(summary_table_bio_groups$Count_Original)
                                     } else { 0 },
  Percentage_Point_Change = NA_real_,
  stringsAsFactors = FALSE
)
summary_table_bio_groups_with_total <- bind_rows(summary_table_bio_groups, total_row_bio_groups)

message("\nSummary Table for Biological Group Distribution (Console Output):")
print(summary_table_bio_groups_with_total, row.names = FALSE)
message("Note: This table provides data for LaTeX Table 1.")


# --- 4. NA Metrics by Feature Calculation Method (for Table 2 summary) ---
message("\n--- Calculating NA Metrics by Feature Calculation Method ---")

na_counts_per_column_analysis <- colSums(is.na(merged_features_only_df))
feature_methods_from_colnames_analysis <- str_extract(names(na_counts_per_column_analysis), "^[^\\.]+")

if(all(is.na(feature_methods_from_colnames_analysis))) {
    stop("ERROR: Could not extract any feature methods from column names of merged_features_only_df. Ensure names are like 'METHOD.FEATURE_X'.")
}

na_summary_df_methods_analysis <- data.frame(
  full_col_name = names(na_counts_per_column_analysis),
  feature_method = feature_methods_from_colnames_analysis,
  na_count = na_counts_per_column_analysis,
  stringsAsFactors = FALSE
)

total_na_by_method_analysis <- na_summary_df_methods_analysis %>%
  filter(!is.na(feature_method)) %>% # Only process rows where a method was successfully extracted
  group_by(feature_method) %>%
  summarise(
    Total_NA_Cells = sum(na_count, na.rm = TRUE),
    Num_Feature_Columns = n(), # Number of columns for this method
    Avg_NAs_Per_Column = mean(na_count, na.rm = TRUE),
    .groups = 'drop' # Good practice for summarise
  ) %>%
  # Add the new percentage of NA cells for this method
  mutate(
    Perc_NA_Cells_for_Method = (Total_NA_Cells / (Num_Feature_Columns * original_row_count_analysis)) * 100
  ) %>%
  # Arrange by average NAs per column (descending) as in your example LaTeX output
  arrange(desc(Avg_NAs_Per_Column))

message("\nSummary Table for NA Metrics by Feature Method (Console Output):")
print(as.data.frame(total_na_by_method_analysis))
message("Note: This table provides data for LaTeX Table 2.")

message("\n--- Script Finished ---") 