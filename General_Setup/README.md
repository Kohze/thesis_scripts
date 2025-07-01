# General Setup Scripts

This directory contains general setup and processing scripts for the thesis analysis.

- `wgbs_keys.txt`: Contains keys or identifiers used in WGBS data processing.
- `rnaseq_keys.txt`: Contains ENCODE accession keys for RNA-seq data from mouse developmental tissue samples, covering embryonic day 10.5 through postnatal day 0 across multiple tissues (forebrain, midbrain, hindbrain, heart, liver, limb, and lung).

## Mouse ENCODE Matrix Dataset

The **Mouse ENCODE Matrix** is a comprehensive developmental time series dataset from the Barbara Wold lab at Caltech, part of the ENCODE Project Consortium. This dataset provides RNA-seq data across multiple mouse tissues during embryonic development.

### Dataset Overview
- **Source**: Barbara Wold lab, Caltech
- **Reference**: ENCODE Project Consortium
- **Data Type**: RNA-seq
- **Species**: Mouse (Mus musculus)
- **Developmental Range**: Embryonic day 10.5 (E10.5) through Postnatal day 0 (P0)

### Tissue Types Included
1. **Forebrain** - Complete time series (E10.5 to P0)
2. **Midbrain** - Complete time series (E10.5 to P0)
3. **Hindbrain** - Complete time series (E10.5 to P0)
4. **Heart** - Complete time series (E10.5 to P0)
5. **Liver** - Partial time series (E11.5 to P0)
6. **Limb** - Partial time series (E10.5 to E15.5)
7. **Lung** - Partial time series (E14.5 to P0)

### Access
- **ENCODE Portal**: [Mouse ENCODE Matrix](https://www.encodeproject.org/mouse-development-matrix/?type=Experiment&status=released&related_series.@type=OrganismDevelopmentSeries&replicates.library.biosample.organism.scientific_name=Mus+musculus)
- **Data Repository**: ENCODE Project Consortium
- **Citation**: ENCODE Project Consortium. (2004). The ENCODE (ENCyclopedia Of DNA Elements) Project. Science, 306(5696), 636-640.

### Usage in Thesis
This dataset is used for developmental time series analysis, providing insights into gene expression patterns across multiple tissues during mouse embryonic development. The RNA-seq data complements the WGBS data for integrated multi-omics analysis. 