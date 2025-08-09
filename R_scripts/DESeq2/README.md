## R Scripts for Differential Gene Expression (DGE) Analysis with DESeq2

This folder contains R scripts used to analyze differential gene expression in *Sycon* sponge samples using **DESeq2**.

### Scripts Included
- `DESeqSycon_body-parts.R`
-   → Analyze gene expression between **osculum** and other body parts on gene level.  
- `DESeqSycon_body-parts_transcripts.R`  
  → Analyze gene expression between **osculum** and other body parts on transcript level, and compare to the translated transcripts that were idenfied in the spicule matrix (to extract those that were also over-expressed in the **osculum** region
  
- `DESeq2-Sycon_regen-I+II.R`  
  → Compare **non-spicule-forming** regeneration stages vs. **spicule-forming** stages.

### Input Data
- Count data is located in:  
  `/input files/count_data/`

- `biomin_genes.tsv`:  
  Maps biomineralization gene names to gene IDs. Used as input for plotting.

### Output
- `Sci_body-parts_over_DEGs_osculum_vs_bodywall_p001-L2FC2.csv`:  
  Contains gene IDs significantly **overexpressed in the osculum**  
  (Log2 fold change ≥ 2, adjusted p-value < 0.01).  
  → Output from `DESeqSycon_body-parts.R`.
