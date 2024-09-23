##DEQ-Analysis to identify transcripts of proteins from spicules that are also over-expressed in spicule-fomring osculm region
##Compare with all info-table to obtain gene IDs for these transcripts 
library("DESeq2")
library("vegan")
warnings()
#setwd #complete your path to repository/R_scripts_DeSeq2 in order to repeat analyses 
#without needing to change paths to all inputfiles. Replace "##write" with "write" in order to 
#produce the output files. 
setwd=("/Users/ovoigt/Dropbox/_Calcarins/submission\ files/repository/R_scripts/DESeq2")

input_counts <- "../../inputfiles/count_data/counts_body_parts/Sci_body-parts_transcript_counts_matrix.tsv"
input_sample_information <- "../../inputfiles/count_data/counts_info/info_Sci_body_parts.csv"

# Read data
DE_Matrix <- read.csv(input_counts, header = TRUE, sep = "\t")
DE_Matrix_INFO <- read.csv(input_sample_information, header = TRUE)

# Filter rows based on counts, only keep those in which at least four samples have higher than 10
DE_Matrix <- DE_Matrix[rowSums(DE_Matrix > 10) >= 4, ]

# Create DESeqDataSet object
DE_MatrixSeq <- DESeqDataSetFromMatrix(countData = DE_Matrix, colData = DE_Matrix_INFO, design = ~condition)
DE_MatrixSeq$condition <- relevel(DE_MatrixSeq$condition, ref = "body_wall")

# Preprocess and transform data
#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
DE_MatrixSeq <- estimateSizeFactors(DE_MatrixSeq)

#get rows with all non-zero counts
non_zero_rows<-apply(counts(DE_MatrixSeq), 1, function(x){all(x>0)})

#get all rows
all_rows<-apply(counts(DE_MatrixSeq), 1, function(x){all(x>=0)})

#number of non-zero rows:
sum(non_zero_rows)

#number of rows:
sum(all_rows)
#cummulative distribution of normalized counts for non-zero rows
multiecdf(counts(DE_MatrixSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))

#density of normalized counts
multidensity(counts(DE_MatrixSeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))
#######

DE_MatrixSeq <- estimateDispersions(DE_MatrixSeq)
DE_MatrixSeq <- nbinomWaldTest(DE_MatrixSeq)

#plot dispersion vs. mean of normalized counts
plotDispEsts(DE_MatrixSeq)
##plotDispEsts(dds)

#PCA

# Apply rlog transformation to the DESeqDataSet DE_MatrixSeq for variance stabilization, ignoring experimental
DE_MatrixSeq_rlog <- rlogTransformation(DE_MatrixSeq, blind = TRUE)


#PCA of samples
# Plot PCA of condition, sample, lib
lapply(c("condition", "sample", "lib"), function(x) plotPCA(DE_MatrixSeq_rlog, intgroup=c(x)))


#wald test
DE_MatrixSeq_Results<-results(DE_MatrixSeq, pAdjustMethod = "BH")

#number of DE genes at 0.01 significance level
table(DE_MatrixSeq_Results$padj < 0.01)
#number of DE genes at 0.01 significance level and log fold change >= 2
table(DE_MatrixSeq_Results$padj < 0.01 & abs(DE_MatrixSeq_Results$log2FoldChange) >= 2)
#number of overexpressed genes at 0.01 significance level and log fold change >= 2
table(DE_MatrixSeq_Results$padj < 0.01 & DE_MatrixSeq_Results$log2FoldChange >= 2)
#number of underexpressed genes at 0.01 significance level and log fold change <= -2
table(DE_MatrixSeq_Results$padj < 0.01 & DE_MatrixSeq_Results$log2FoldChange <= -2)
hist(DE_MatrixSeq_Results$pvalue, main = "Treatment vs. Control", xlab="p-values")
#plotMA(DE_MatrixSeq_Results, alpha=0.01)

####
#test for global differences in expression using the distance matrix and the function adonis2.
#####

#produce a matrix with the transformed distances
osculum_bodywall <- data.frame(Condition=c("body_wall" , "osculum" , "body_wall" , "body_wall" , "osculum" , "body_wall" , "body_wall" , "osculum" , "body_wall" , "body_wall" , "osculum" , "body_wall" , "body_wall" , "osculum" , "body_wall"))
#now use the function adonis to test for differences between treatments. By default adonis uses the bray-curtis distance.
adonis2(t(assay(DE_MatrixSeq_rlog))~Condition, data=osculum_bodywall , method="euclidean")
#######


# Get differential expressed genes #combine with next section to use the same variable.
overexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange >= 2)
over_deg_names<-rownames(overexpressed)

# Get all_info for the peptide_IDs found in the spicules
spic_prot_all_info <- read.csv("../../inputfiles/peptides_spicules/spicule_peptides1FDR-min_pep_IDs_from_Sci_HBWS01_all_info.tsv", header=TRUE, sep="\t")

#Keep only those transcripts in spic_prot_all_info that are over-expressed in osculum region (increased spicule formation)
spic_prot_all_info_overexp <-spic_prot_all_info[spic_prot_all_info$Transcript_accession %in% over_deg_names, ]
write.table(spic_prot_all_info_overexp, "Spic_prots+overexp_osc_all_info.tsv", sep="\t",  col.names=T, row.names = F, quote = F)
#####END HERE

