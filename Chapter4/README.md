# Chapter 4: Gene Isoform Diversity and Feature Analysis

This directory contains scripts for the analyses presented in Chapter 4 of the thesis, focusing on gene isoform diversity and its relationship with various gene features.

## Scripts Overview

The scripts are designed to be run sequentially:

1.  **`01_isoform_stats.R`**: 
    *   **Purpose**: Calculates statistics related to isoform usage per gene (e.g., number of isoforms, primary isoform usage). Merges these stats with gene feature data.
    *   **Inputs**: 
        *   `data/uniprot_isoforms_filtered.RDS` (Filtered isoform data)
        *   `data/processed/gene_features_filtered.RDS` (Filtered gene features)
    *   **Outputs**: 
        *   `data/processed/isoform_stats_with_features.RDS` (Combined data frame)
    *   **Packages**: `dplyr`, `tidyr`

2.  **`02_pearson_corr.R`**: 
    *   **Purpose**: Performs Pearson correlation analysis between isoform statistics (from script 01) and continuous gene features. Generates correlation plots and saves correlation results.
    *   **Inputs**: 
        *   `data/processed/isoform_stats_with_features.RDS` (Output from script 01)
    *   **Outputs**:
        *   `results/chapter4/isoform_feature_correlations.pdf` (Correlation plots)
        *   `results/chapter4/isoform_feature_correlation_table.csv` (Correlation coefficients and p-values)
    *   **Packages**: `dplyr`, `ggplot2`, `ggpubr`, `Hmisc`, `corrplot`, `broom` 

3.  **`03_isoform_classification.R`**:
    *   **Purpose**: Classifies genes based on the number of isoforms (e.g., single vs. multi-isoform). Performs statistical tests (Wilcoxon rank-sum test) to compare gene features between these classes and generates comparison boxplots.
    *   **Inputs**:
        *   `data/processed/isoform_stats_with_features.RDS` (Output from script 01, assumed to be loaded as `gene_features`)
    *   **Outputs**:
        *   `results/chapter4/gene_feature_by_isoform_class_boxplots.pdf` (Boxplots comparing features)
        *   `results/chapter4/gene_feature_by_isoform_class_stats.csv` (Statistical test results)
    *   **Packages**: `ggplot2`, `dplyr`, `ggpubr`, `tidyr`, `broom`

## Workflow

1.  Ensure the required input files are present in the specified `data/` subdirectories.
2.  Run the scripts in the numbered order (01, 02, 03).
3.  Check the `results/chapter4/` directory for output plots and tables.

## Required Packages

Make sure the following R packages are installed:
- `dplyr`
- `tidyr`
- `ggplot2`
- `ggpubr`
- `Hmisc`
- `corrplot`
- `broom` 