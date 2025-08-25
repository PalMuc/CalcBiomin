# CalcBiomin

This repository contains scripts and input data for the analyses described in the paper:  
**“Genetic parallels in biomineralization of the calcareous sponge *Sycon ciliatum* and stony corals”**
by Voigt O, Wilde MV, Fröhlich T, Fradusco B, Vargas S, Wörheide G.

A DOI will be added as soon as it is available.

The repository is organized to replicate RNA-seq analyses (DESeq2 DGE, WGCNA), annotation, and GO enrichment.

---

## Repository structure

~~~
inputfiles/                     # Input data used in the analyses
    peptides_spicules/          # Peptide annotation data for spicule proteins
    annotation_files/           # Gene/transcript annotation data
        petide_blast_and_annotation/
    count_data/                 # # Read count tables and metadata (see details below)
        counts_regeneration/
        gene_counts_combined/
        counts_body_parts/
        counts_info/
    GOterms/                    # GO terms obtained with the script perl/Get_GO_Annotations_v2.pl
        gene_GO_terms/
        transcript_GOterms/

R_scripts/                      # R scripts for statistical analyses
    DESeq2/                     # Differential gene expression with DESeq2
        plots_DESeq_body-parts/
        plots_DESeq_regeneration/
    GOs/                        # GO term enrichment analysis
    wcgna/                      # WGCNA co-expression analysis
        plots/

perl/                           # Perl scripts for preprocessing and annotation
transcriptome/                  # Transcriptome reference data
~~~
## Count Matrices for DGE and WGCNA Analyses

Five specimens of *Sycon ciliatum* were dissected into three body parts:

1. Oscular region  
2. Inner sponge wall  
3. Outer sponge wall  

RNA-seq was conducted for each body part, and the raw reads were mapped against a *Sycon* transcriptome. Mapped reads were further filtered to remove sequences from commensal organisms. Gene and transcript counts for each filtered set were obtained with **SALMON** and combined into count matrices for the body parts experiment and for the regeneration experiment.

- Count matrices for the **body parts** experiment are in:  
  `inputfiles/count_data/counts_body_parts/`

- For the **regeneration** experiment, raw reads (PRJNA628727) were processed in the same way. Count matrices are in:  
  `inputfiles/count_data/counts_regeneration/`

- A matrix combining the **gene counts** of both experiments is in:  
  `inputfiles/count_data/gene_counts_combined/`

- Sample information for these matrices is in:  
  `inputfiles/count_data/counts_info/`


---

## Requirements

- **R** (≥ 4.0)
- R packages commonly used across scripts:  
  `DESeq2`, `genefilter`, `WGCNA`, `topGO`, `ggplot2`, `reshape2`, `svglite`  
  *(Some scripts may require additional packages; see comments within each script.)*
- **Perl** (for scripts in the `perl/` directory)

---

## Input files

All required input files are included under `inputfiles/`, including:
- Count data tables (gene/transcript level) and sample metadata
- Gene/transcript annotation and BLAST results
- GO term mapping files

---

## How to run

1. Open the relevant script from `R_scripts/` or `perl/`.
2. Follow the comments in each script regarding working directories and paths  
   (several scripts set their own working directory and use repo-relative paths).
3. Run from R or the command line, e.g.:

   ~~~bash
   Rscript R_scripts/<subfolder>/<script_name>.R
   ~~~

---

## Outputs

Depending on the script:
- Differential expression result tables (CSV)
- Plots (SVG/PDF) saved in subfolders such as `R_scripts/DESeq2/plots_*` or `R_scripts/wcgna/plots`
- GO enrichment result tables
- WGCNA module gene lists

---

