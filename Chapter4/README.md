# Chapter 4: Gene Isoform Diversity and Feature Analysis

This directory contains scripts for the analyses presented in Chapter 4 of the thesis, focusing on gene isoform diversity and its relationship with various gene features, as well as the primary Differential Transcript Usage (DTU) analysis.

## Scripts Overview

The scripts are designed to be run sequentially:

1.  **`01_run_isoformSwitchAnalyzeR.R`**:
    *   **Purpose**: Performs the core DTU analysis using the `isoformSwitchAnalyzeR` package to identify significant isoform switching events across timepoints for Forebrain, Midbrain, and Hindbrain.
    *   **Inputs**: Kallisto quantification results, annotation file (GTF/GFF3), genome sequence (FASTA).
    *   **Outputs**: `switchAnalyzeRlist` RDS objects (e.g., `Forebrain_switchAnalyzeRlist_final.rds`) containing detailed analysis results for each tissue. Also optional CSV tables of significant switches and consequences.
    *   **Packages**: `isoformSwitchAnalyzeR`, `dplyr`, `Biostrings`, `ggplot2`, `DEXSeq`.

2.  **`02_isoform_stats.R`**: 
    *   **Purpose**: Calculates statistics related to isoform usage per gene (e.g., number of isoforms, primary isoform usage). Merges these stats with gene feature data.
    *   **Inputs**: 
        *   `data/uniprot_isoforms_filtered.RDS` (Filtered isoform data)
        *   `data/processed/gene_features_filtered.RDS` (Filtered gene features)
    *   **Outputs**: 
        *   `data/processed/isoform_stats_with_features.RDS` (Combined data frame)
    *   **Packages**: `dplyr`, `tidyr`

3.  **`03_pearson_corr.R`**: 
    *   **Purpose**: Performs Pearson correlation analysis between isoform statistics (from script 02) and continuous gene features. Generates correlation plots and saves correlation results.
    *   **Inputs**: 
        *   `data/processed/isoform_stats_with_features.RDS` (Output from script 02)
    *   **Outputs**:
        *   `results/chapter4/isoform_feature_correlations.pdf` (Correlation plots)
        *   `results/chapter4/isoform_feature_correlation_table.csv` (Correlation coefficients and p-values)
    *   **Packages**: `dplyr`, `ggplot2`, `ggpubr`, `Hmisc`, `corrplot`, `broom` 

4.  **`04_isoform_classification.R`**:
    *   **Purpose**: Classifies genes based on the number of isoforms (e.g., single vs. multi-isoform). Performs statistical tests (Wilcoxon rank-sum test) to compare gene features between these classes and generates comparison boxplots.
    *   **Inputs**:
        *   `data/processed/isoform_stats_with_features.RDS` (Output from script 02, assumed to be loaded as `gene_features`)
    *   **Outputs**:
        *   `results/chapter4/gene_feature_by_isoform_class_boxplots.pdf` (Boxplots comparing features)
        *   `results/chapter4/gene_feature_by_isoform_class_stats.csv` (Statistical test results)
    *   **Packages**: `ggplot2`, `dplyr`, `ggpubr`, `tidyr`, `broom`

## Workflow

1.  Ensure the required input files are present in the specified `data/` subdirectories.
2.  Run the scripts in the numbered order (01, 02, 03, 04).
3.  Check the `results/chapter4/` directory for output plots and tables.

## Required Packages

Make sure the following R packages are installed:
- `isoformSwitchAnalyzeR`
- `DEXSeq`
- `Biostrings`
- `dplyr`
- `tidyr`
- `ggplot2`
- `ggpubr`
- `Hmisc`
- `corrplot`
- `broom` 

# Chapter 4: `isoformSwitchAnalyzeR` Analysis README

## Purpose

This directory contains the R script (`01_run_isoformSwitchAnalyzeR.R`) responsible for performing the core Differential Transcript Usage (DTU) analysis for the RNA-seq data presented in Chapter 4 of the thesis, "The Elusive Dynamic of the Genome."

The script uses the `isoformSwitchAnalyzeR` package to identify genes exhibiting significant changes in isoform usage across the developmental time course for Forebrain, Midbrain, and Hindbrain tissues separately.

## Prerequisites

1.  **R Environment:** A working R installation (version 4.0 or later recommended) with the `isoformSwitchAnalyzeR` package and its dependencies (including `dplyr`, `Biostrings`, `ggplot2`, `DEXSeq`) installed.
2.  **External Tools (Optional):** If sequence-based consequence analyses (e.g., protein domains, coding potential, signal peptides) are desired, the corresponding external tools must be installed and configured, and their paths correctly specified in the R script:
    *   [CPAT](https://cpat.readthedocs.io/en/latest/)
    *   [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-5.0) (Requires academic license)
    *   [PfamScan](https://www.ebi.ac.uk/Tools/pfa/pfamscan/) (Requires Pfam database download)
    *   [NetSurfP-2](https://services.healthtech.dtu.dk/service.php?NetSurfP-2.0)
3.  **Input Data:** The necessary input data files must be available (see below).

## Input Data

The script requires the following input data:

1.  **Transcript Quantification Results:** Output files from a transcript quantification tool, typically [Kallisto](https://pachterlab.github.io/kallisto/) or [Salmon](https://salmon.readthedocs.io/en/latest/). The script expects one directory per sample, containing the `abundance.tsv` (or `quant.sf`) file. The paths to these directories need to be correctly configured in the `kallisto_path_template` and `sample_data` variables within the script.
2.  **Annotation File:** A gene annotation file in GTF or GFF3 format, corresponding to the reference transcriptome used for quantification. The path should be specified in the `annotation_path` variable.
3.  **Genome Sequence File:** A FASTA file containing the genome sequence corresponding to the annotation file. The path should be specified in the `fasta_path` variable.

## Script Usage (`01_run_isoformSwitchAnalyzeR.R`)

1.  **Configuration:**
    *   **CRITICAL:** Before running, you **must** edit the `01_run_isoformSwitchAnalyzeR.R` script to provide the correct paths to your input data (quantification directories, annotation file, genome FASTA) and any optional external tools.
    *   Verify and customize the `sample_data` data frame to accurately reflect your sample IDs, their corresponding conditions (e.g., timepoints), tissue types, and Kallisto output directory names.
    *   Adjust analysis parameters (e.g., `alpha_value`, `dIF_cutoff`, `analysis_threads`, filtering cutoffs) as needed for your specific analysis.
    *   Set the boolean flags (`run_orf_analysis`, `run_cpat_analysis`, etc.) to `TRUE` only for the sequence analyses you wish to perform *and* have the necessary external tools installed and configured.
2.  **Execution:**
    *   Navigate to the `thesis_scripts/chapter4` directory in your terminal.
    *   Run the script using R: `Rscript 01_run_isoformSwitchAnalyzeR.R`

## Output

The primary outputs of the script are saved in the specified `output_dir` (default: `./thesis_scripts/chapter4/isoformSwitchAnalyzeR_output`):

1.  **`switchAnalyzeRlist` Objects (RDS):** For each tissue processed (Forebrain, Midbrain, Hindbrain), a final `switchAnalyzeRlist` object is saved as an `.rds` file (e.g., `Forebrain_switchAnalyzeRlist_final.rds`). This object contains all imported data and the results of the filtering, statistical testing, and consequence analyses. These RDS files are the primary inputs for downstream analysis scripts (e.g., those used in `fullscript.txt` to generate `combinedForebrain`, etc.).
2.  **Result Tables (CSV - Optional):**
    *   `*_significant_switches.csv`: A table listing the statistically significant isoform switches identified by `isoformSwitchTestDEXSeq`.
    *   `*_switch_consequences.csv`: A table detailing the predicted functional consequences of the identified switches.
3.  **External Tool Output (Optional):** If external tools were run, their raw output files will be saved in subdirectories within the `output_dir` (e.g., `Forebrain_external_analsysis`).

## Structure

The script processes each tissue (Forebrain, Midbrain, Hindbrain) sequentially within a loop. For each tissue, it performs the full workflow: data import, filtering, statistical testing, optional sequence analysis, consequence analysis, and saving results. 