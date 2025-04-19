# Thesis Scripts: The Elusive Dynamic of the Genome

## Overview

This repository contains the R scripts supporting the analyses presented in the thesis "The Elusive Dynamic of the Genome: Novel Sequence Patterns, Splicing and Regulation" by Robin Kohze (University of Cambridge, 2025).

The scripts cover investigations into:

*   **Novel Open Reading Frame (nORF) characterization:** Including proteochemometric analysis and machine learning classification (Chapter 3).
*   **Gene Isoform Diversity:** Analyzing isoform usage, correlations with gene features, and comparing single vs. multi-isoform genes (Chapter 4).
*   **Epigenetic Regulation of Splicing:** Examining the roles of DNA methylation, chromatin states, DNA shape, and sequence motifs in alternative splicing during mouse brain development (Chapter 5).

## Repository Structure

The scripts are organized into subdirectories corresponding to thesis chapters or general setup tasks. Each chapter directory contains scripts numbered (`0X_script_name.R`) according to their typical execution order within that chapter's workflow.

*   **[General_Setup/](./General_Setup/)**: Scripts for metadata processing, WGBS key management, and initial methylation data handling. See the [General Setup README](./General_Setup/README.md).
*   **[Chapter2/](./Chapter2/)**: Contains a reference to the external code repository for the nORFs.org platform. See the [Chapter 2 README](./Chapter2/README.md).
*   **[Chapter3/](./Chapter3/)**: Scripts for the high-dimensional analysis of nORFs, including PCM feature extraction, correlation analysis, XGBoost classification, topology analysis, and nORF kernel development. See the [Chapter 3 README](./Chapter3/README.md).
*   **[Chapter4/](./Chapter4/)**: Scripts analyzing gene isoform diversity, its correlation with gene features, and comparing single- vs. multi-isoform gene characteristics. See the [Chapter 4 README](./Chapter4/README.md).
*   **[Chapter5/](./Chapter5/)**: Scripts investigating the interplay between epigenetics and splicing, including methylation/expression correlations, divergence analysis, HMM state enrichment, motif enrichment, DNA shape analysis, and Gviz plotting. See the [Chapter 5 README](./Chapter5/README.md).

Please refer to the `README.md` file within each subdirectory for detailed information on the specific inputs, outputs, workflow, and package requirements for the scripts therein.

## Dependencies

These scripts rely on numerous R packages from CRAN and Bioconductor. Key packages include `dplyr`, `ggplot2`, `GenomicRanges`, `SummarizedExperiment`, `Gviz`, `monaLisa`, `xgboost`, `stats`, and others specific to bioinformatics and data analysis.

**Recommendation:** Use `renv` for managing package dependencies to ensure reproducibility. If a `renv.lock` file is present, run `renv::restore()` in R to install the exact package versions used.

## Data

Input data (e.g., sequence files, annotation files, methylation calls, expression matrices, HMM states) are generally assumed to be available locally and are referenced within the scripts. Check individual script headers and chapter READMEs for specific input requirements and expected data formats.

## Usage

Consult the README file within each chapter directory for the recommended workflow and dependencies specific to that chapter's analyses.

## License

MIT

## Contact

Robin Kohze - [GitHub Profile](https://github.com/Kohze)
