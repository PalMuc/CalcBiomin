library("DESeq2")
library("RColorBrewer")
library("vegan")
library("geneplotter")
library("ComplexHeatmap")
library('RColorBrewer')
library('circlize')

#setwd #complete your path to repository/R_scripts_DeSeq2 in order to repeat analyses 
#without needing to change paths to all inputfiles. Replace "##write" with "write" in order to 
#produce the output files. 

input_counts <- "../../inputfiles/counts_regeneration/Sci_regen_setI+II_gene_counts.tsv"
input_sample_information <- "../../inputfiles/counts_info/info_Sci_regenI+II.csv"

# Read data
DE_Matrix <- read.csv(input_counts, header = TRUE, sep = "\t")
DE_Matrix_INFO <- read.csv(input_sample_information, header = TRUE)

# Filter rows based on counts
DE_Matrix <- DE_Matrix[rowSums(DE_Matrix > 10) >= 4, ]

# Create DESeqDataSet object
DE_MatrixSeq <- DESeqDataSetFromMatrix(countData = DE_Matrix, colData = DE_Matrix_INFO, design = ~condition)
DE_MatrixSeq$condition <- relevel(DE_MatrixSeq$condition, ref = "no_spic")

# Preprocess and transform data
#estimate size factors
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

# Rlog, generate heatmap and PCA

# Apply rlog transformation to the DESeqDataSet DE_MatrixSeq for variance stabilization, ignoring experimental
DE_MatrixSeq_rlog <- rlogTransformation(DE_MatrixSeq, blind = TRUE)

#PCA of samples
# Plot PCA of condition, sample, lib
lapply(c("condition", "lib","day"), function(x) plotPCA(DE_MatrixSeq_rlog, intgroup=c(x)))

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

# Get differential expressed genes #combine with next section to use the same variable.
overexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange >= 2)
over_deg_names<-rownames(overexpressed)

#overexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange >= 2)
underexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange <= -2)
top_overexpressed <- head(overexpressed[order(overexpressed$log2FoldChange, decreasing = TRUE), ], 50)
top_underexpressed <- head(underexpressed[order(underexpressed$log2FoldChange, decreasing = TRUE), ], 50)
top_genes <- rbind(top_overexpressed, top_underexpressed)
top_genes <- top_genes[order(top_genes$log2FoldChange, decreasing = TRUE), ]

#Creating heatmap
#Plot heatmap of 50 over and underexpressed genes
heatmap_data <- assay(DE_MatrixSeq_rlog)[rownames(top_genes), , drop = F]
#rename rows of libs with experiment info:
day_labels <- DE_Matrix_INFO$day[match(colnames(heatmap_data), rownames(DE_Matrix_INFO))]
# Rename the columns of heatmap_data
colnames(heatmap_data) <- day_labels

pheatmap(heatmap_data, col = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),
         scale = "row", cluster_rows = TRUE, border_color = NA)

#get significant genes names; uncomment the write.csv lines and change the path if write to a file is required/wanted
deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 & abs(log2FoldChange) >= 2)
##write.csv(deg_all_info, "Sci_regen-I+II_DEGs_no-spics_vs_spics_p001-L2FC2_all_info.csv")
degs_names<-rownames(subset(DE_MatrixSeq_Results, padj<0.01 & abs(log2FoldChange) >= 2))
#write.csv(degs_names, "Sci_regen-I+II_DEGs_no-spics_vs_spics_p001-LFC2.csv")###CLEAN
##write.table(degs_names, "Sci_regen-I+II_DEGs_no-spics_vs_spics_p001-LFC2.csv", sep=",",  col.names=F, row.names = F)

#get significantly overexpressed in treatment ##(here spic):
over_deg_names<-rownames(overexpressed)
##write.csv(overexpressed, "Sci_regen-I+II_DEGs_spic_over_vs_no-spics_p001-L2FC2_all_info.csv")
##write.table (over_deg_names, "Sci_regen-I+II_DEGs_spic_over_vs_no-spics_p001-L2FC2.csv", sep=",",  col.names=F, row.names = F)

#get significantly underexpressed in treatment: (##here:spic)
under_deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 & log2FoldChange <= -2)
under_deg_names<-rownames(under_deg_all_info)
##write.csv(under_deg_all_info, "Sci_regen-I+II_DEGs_spic_under_vs_no-spic_p001-L2FC2_all_info.csv")
##write.table ( over_deg_names, "Sci_regen-I+II_DEGs_spic_under_vs_no-spic_p001-L2FC2.csv", sep=",",  col.names=F, row.names = F)

####
#test for global differences in expression using the distance matrix and the function adonis2.
#####

#produce a matrix with the transformed distances
spic_nospic <- data.frame(Condition = DE_Matrix_INFO$condition)

#now use the function adonis to test for differences between treatments. By default adonis uses the bray-curtis distance.
adonis2(t(assay(DE_MatrixSeq_rlog))~Condition, data=spic_nospic, method="euclidean")

#plot complex-heatmap of Z-scores for comparison of changes in expression of genes per sample
## for biomin genes

#condition <-factor(c("osculum","osculum","osculum","osculum","osculum","body_wall","body_wall","body_wall","body_wall","body_wall","body_wall","body_wall","body_wall","body_wall","body_wall"))
#coldata <- data.frame(row.names = colnames(DE_MatrixSeq), condition)#
coldata <- data.frame(row.names = colnames(DE_MatrixSeq), DE_Matrix_INFO$condition)#
res <- results(DE_MatrixSeq, contrast= c ("condition", "no_spic", "spic"))
df<- as.data.frame(res)
biomingenename_map <- read.csv('biomin_genes.tsv', header=FALSE, sep="\t")

keys <- biomingenename_map$V2
values <- biomingenename_map$V1

l <-list()
for (i in 1:length(keys)){
  l[keys[i]] <- values[i]
} 
#for non-mapped genes
no_values <- setdiff(rownames(df), keys)
for (i in 1:length(no_values)){l[no_values[i]] <- no_values[i]
}
df$symbol <- unlist(l[rownames(df)], use.names= FALSE)

# Subset 'df' using the specific gene identifiers
biomingenes <- df[biomingenename_map$V2, ]
# Use the rownames of 'biomingenes' to subset the assay data
mat <- assay(DE_MatrixSeq_rlog)[rownames(biomingenes), rownames(coldata)]

base_mean <- rowMeans(mat)
# Apply z-score normalization to each row (gene) in 'mat' and transpose the result
mat.scaled <- t(apply(mat, 1, scale)) 
colnames (mat.scaled)<- colnames(mat)

#getting log2fold-change for the genes of interest
l2_val <- as.matrix(biomingenes$log2FoldChange) 
colnames(l2_val) <- "logFC"
##get rows from df of for the genes of interest
print(rownames(biomingenes) %in% rownames(df))#test if all genes of interest occur in df
df_biomin <- df[rownames(biomingenes), ]

# Convert baseMean values into a matrix and rename the column
mean <- as.matrix(df_biomin$baseMean) 
colnames(mean) <- "AveExpr"

# Define color schemes for the log fold-change values and average expression
col_logFC <-colorRamp2(c(min(l2_val), 0, max (l2_val)), c("white", "salmon", "red"))
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "orange"))

# Define the annotation for the heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp=gpar(fill=2), height= unit(2, "cm")))

# Generate three heatmaps: one for Z-scores, one for logFC and one for AveExpr
# Match the column names of mat.scaled with the rownames of DE_Matrix_INFO
# to get the corresponding "day" values
day_labels <- DE_Matrix_INFO$day[match(colnames(mat.scaled), rownames(DE_Matrix_INFO))]
h1 <- Heatmap(mat.scaled, cluster_rows=F,column_labels=day_labels, name="Z-score",cluster_columns=T)

#h1 <- Heatmap(mat.scaled, cluster_rows=F,column_labels=colnames(mat.scaled), name="Z-score",cluster_columns=T)
h2 <- Heatmap(l2_val, row_labels = biomingenes$symbol, cluster_rows =F, name="logFC", top_annotation= ha, col =  col_logFC, cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i,j],2),x,y)})

h3 <- Heatmap(mean, row_labels =biomingenes$symbol, cluster_rows =F, name = "AveExpr", col= col_AveExpr, cell_fun =function(j,i,x,y,w,h,col){grid.text(round(mean[i,j],2),x,y)})

# Combine the three heatmaps into one
h <- h1+h2+h3

# Print the heatmap
h



