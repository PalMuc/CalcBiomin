# **Count Matrices for DGE and WCGNA Analyses**

Five specimens of *Sycon ciliatum* were dissected into three body parts:

1. Oscular region
2. Inner sponge wall
3. Outer sponge wall

RNAseq was conducted for each body part, and we mapped the raw reads against a *Sycon* transcriptome. Only mapping reads were further processed to exclude reads from commensals. Gene and transcript counts for each filtered set were obtained with SALMON and combined into count matrices for the body parts experiment and the regeneration data. The count matrices are provided in the folder **counts\_body\_parts/**

We downloaded raw reads from a regeneration experiment (PRJNA628727) and processed them in the same way. The count matrices are provided in the folder **counts\_regeneration**/

A matrix combining the gene counts of both experiments is provided in **gene\_counts\_combined/**

Sample information of these matrices is included in the **folder count\_info/**
