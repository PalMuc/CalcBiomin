# WGCNA on DESeq2-normalized Sycon RNA-seq

This script runs variance-stabilized normalization with DESeq2 (v1.42.1) and uses WGCNA module detection, plotting, and gene list export. The published data were analyzed with DESeq2 v1.42.1 and WGCNA v1.72.5.

## Requirements
- R  
- R packages: `DESeq2`, `genefilter`, `WGCNA`, `topGO`, `ggplot2`, `reshape2`, `svglite`

## Input files
All required input files are included in this repository:
- `inputfiles/count_data/gene_counts_combined/Sci_gene_counts_combined.tsv`  
  → Gene count matrix (tab-separated)
- `inputfiles/count_data/counts_info/info_Sci_body_parts+regenI+II.csv`  
  → Sample info with `condition` column

## How to run
From the **repository root**:
```r
source("R_scripts/wgcna/wgcna.R")```