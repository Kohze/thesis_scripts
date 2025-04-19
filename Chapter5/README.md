# Chapter 5: Epigenetic Regulation and DNA Shape

This directory contains R scripts for the analyses presented in Chapter 5 of the thesis, focusing on the interplay between DNA methylation, gene expression, chromatin states, motif occurrences, and DNA shape.

## Workflow Overview

The scripts generally follow this analysis flow:

1.  `01_methylation_expression.R`: Analyze the basic relationship between methylation and expression.
2.  `02_methylation_expression_divergence.R`: Analyze the variability (divergence) of methylation and expression across conditions.
3.  `03_hmm_enrichment.R`: Perform enrichment analysis of features within chromatin HMM states.
4.  `04_monaLisa.R`: Identify enriched sequence motifs within specific genomic regions (e.g., HMM states, DMRs).
5.  `05_shape_analysis_post_analysis.R`: Analyze DNA shape parameters, possibly in relation to motifs or regions identified earlier.
6.  `06_gviz_plotting.R`: Visualize results by combining different data tracks for specific genomic loci.

## Input Data Requirements

These scripts assume various data objects are loaded into the R environment or read from files. Key potential inputs include:

*   **Methylation Data**: Gene/region-level methylation values (e.g., `methylation_data` df or `methylation_matrix`).
*   **Expression Data**: Gene expression values (e.g., `expression_data` df or `expression_matrix`).
*   **Gene Annotations**: Data frame mapping IDs to features like biotype (`gene_annotations`).
*   **HMM States**: GRanges object with chromatin state calls (`hmm_states_gr`).
*   **Genomic Features**: Named list of GRanges objects for enrichment analysis (`features_list`).
*   **Background Regions**: GRanges defining the background for enrichment (`background_regions`).
*   **Target Regions (Motif)**: GRanges list for monaLisa (`target_regions`).
*   **Sequence Data**: BSgenome object or DNAStringSet for genome sequence (`sequence_data`).
*   **Motif Database**: List or PWMatrixList of motifs (`motif_database`).
*   **Shape Predictions**: Output from DNA shape prediction tools (`shape_predictions`).
*   **Regions for Gviz**: Data frame/GRanges specifying regions to plot (`regions_to_plot`).
*   **Genome Assembly ID**: String like "hg38" or "hg19" (`genome_assembly`).

Refer to individual script headers for more specific input expectations.

## Scripts

1.  **`01_methylation_expression.R`**:
    *   **Purpose**: Analyzes correlation and group differences between DNA methylation and gene expression.
    *   **Actions**: Merges data, calculates Pearson correlation, performs t-tests/Wilcoxon tests between high/low methylation groups, generates scatter and box plots.
    *   **Key Output**: Correlation plots, box plots, summary statistics CSV.

2.  **`02_methylation_expression_divergence.R`**:
    *   **Purpose**: Analyzes the variability (e.g., standard deviation) of methylation and expression across samples/conditions.
    *   **Actions**: Calculates SD for methylation and expression per gene, correlates these divergence metrics, visualizes divergence vs. mean levels or other features.
    *   **Key Output**: Divergence statistics CSV, correlation plots, divergence-vs-feature plots.

3.  **`03_hmm_enrichment.R`**:
    *   **Purpose**: Calculates enrichment of genomic features within chromatin HMM states.
    *   **Actions**: Uses base-pair overlap and Fisher's exact test to assess enrichment for each feature type in each HMM state compared to a background.
    *   **Key Output**: Enrichment statistics CSV (including fold enrichment, p-value, FDR), optional heatmap visualization.

4.  **`04_monaLisa.R`**:
    *   **Purpose**: Performs motif enrichment analysis in specified target regions using `monaLisa`.
    *   **Actions**: Runs `calcMotifEnrichment` considering background regions, generates summary tables and plots (empirical p-values, scores, heatmaps) for enriched motifs.
    *   **Key Output**: Full monaLisa results RDS object, summary plots PDF, heatmap PDF.

5.  **`05_shape_analysis_post_analysis.R`**:
    *   **Purpose**: Post-processes DNA shape predictions (e.g., from `DNAshapeR`).
    *   **Actions**: Summarizes shape parameters (MGW, ProT, Roll, HelT) around centers of regions of interest (e.g., motifs), generates average shape profile plots. (*Note: Core summarizing function is a placeholder requiring implementation based on specific shape data format*).
    *   **Key Output**: Shape summary CSV, average profile plots PDF.

6.  **`06_gviz_plotting.R`**:
    *   **Purpose**: Creates multi-track genomic visualizations for specific regions using `Gviz`.
    *   **Actions**: Defines and combines Gviz tracks (genes, methylation, HMM states, motifs, etc.), plots specified genomic regions with these tracks.
    *   **Key Output**: PDF plots for each specified genomic region.

## Required Packages

*   **Core**: `dplyr`, `tidyr`, `stats`
*   **Genomic Ranges**: `GenomicRanges`, `rtracklayer` (potentially)
*   **Bioinformatics**: `SummarizedExperiment`, `BiocParallel`, `BSgenome.*` (specific assembly), `monaLisa`, `DNAshapeR` (if used for input), `TxDb.*` or `EnsDb.*` (potentially for gene models)
*   **Statistics**: `matrixStats`, `broom`, `Hmisc` (potentially)
*   **Plotting**: `ggplot2`, `ggpubr`, `Gviz`, `RColorBrewer`, `patchwork` (potentially)
*   **Tables**: `xtable` (potentially), `kableExtra` (potentially)

Please ensure all necessary packages are installed. 