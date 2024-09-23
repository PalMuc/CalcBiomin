library("DESeq2")
library("RColorBrewer")
library("vegan")
library("geneplotter")
library("ComplexHeatmap")
library("circlize")

warnings()
#setwd #complete your path to repository/R_scripts_DeSeq2 in order to repeat analyses 
#without needing to change paths to all inputfiles. Replace "##write" with "write" in order to 
#produce the output files. 
#setwd("/repository/R_scripts/DESeq2")
input_counts <- "../../inputfiles/count_data/counts_body_parts/Sci_body-parts_gene_counts.tsv"
input_sample_information <- "../../inputfiles/count_data/counts_info/info_Sci_body_parts.csv"

# Read data
DE_Matrix <- read.csv(input_counts, header = TRUE, sep = "\t")
DE_Matrix_INFO <- read.csv(input_sample_information, header = TRUE)

# Filter rows based on counts
DE_Matrix <- DE_Matrix[rowSums(DE_Matrix) > 5, ]

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
DE_MatrixSeq <- estimateDispersions(DE_MatrixSeq)
DE_MatrixSeq <- nbinomWaldTest(DE_MatrixSeq)

#plot dispersion vs. mean of normalized counts
plotDispEsts(DE_MatrixSeq)


# Apply rlog transformation to the DESeqDataSet DE_MatrixSeq for variance stabilization, ignoring experimental
DE_MatrixSeq_rlog <- rlogTransformation(DE_MatrixSeq, blind = TRUE)

#export the rlog count for downstream analysis
rlog_out <- assay (DE_MatrixSeq_rlog)#[rownames(DEMatrixSeq), rownames(coldata)]
##write.csv(rlog_out, file = "Sci_body-parts_rlog_counts.csv")

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

#overexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange >= 2)
underexpressed <- subset(DE_MatrixSeq_Results, padj < 0.01 & log2FoldChange <= -2)
top_overexpressed <- head(overexpressed[order(overexpressed$log2FoldChange, decreasing = TRUE), ], 50)
top_underexpressed <- head(underexpressed[order(underexpressed$log2FoldChange, decreasing = TRUE), ], 50)
top_genes <- rbind(top_overexpressed, top_underexpressed)
top_genes <- top_genes[order(top_genes$log2FoldChange, decreasing = TRUE), ]


# Generate heatmap
###
#Plot heatmap of 50 over and underexpressed genes
heatmap_data <- assay(DE_MatrixSeq_rlog)[rownames(top_genes), , drop = F]

# replace geneIDs with biomingenes for heatmap
#Get map Biomingene symbols to GeneID
biomingenename_map <- read.csv('biomin_genes.tsv', header=FALSE, sep="\t")

#Identify geneIDs of biomingenes in heatmap.data
overlapping_gene_ids <- intersect(row.names(heatmap_data), biomingenename_map$V2)

# Create a lookup vector using the mapping matrix for overlapping gene IDs
lookup_vector <- biomingenename_map$V1[match(overlapping_gene_ids, biomingenename_map$V2)]

# Replace row names with gene symbols for overlapping gene IDs
rownames(heatmap_data)[rownames(heatmap_data) %in% overlapping_gene_ids] <- lookup_vector

#Creating heatmap
pheatmap(heatmap_data, col = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)),
         scale = "row", cluster_rows = TRUE)

#Get all genes names and info of significant DEG
deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 )
##write.csv(deg_all_info, "Sci_body-parts_DEGs_osculum_vs_bodywall_p001_all_info.csv");

#Get all genes names and info of significant DEG with log2FOLDchange >=2 
deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 & abs(log2FoldChange) >= 2)
##write.csv(deg_all_info, "Sci_body-parts_DEGs_osculum_vs_bodywall_p001-L2FC2_all_info.csv")
degs_names<-rownames(subset(DE_MatrixSeq_Results, padj<0.01 & abs(log2FoldChange) >= 2))
##write.table(degs_names, "Sci_body-parts_DEGs_osculum-vs_bodywall_p001-LFC2.csv", sep=",",  col.names=F, row.names = F)

#Get names and info of DEGs significantly overexpressed in "treatment" (here osculum):
over_deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 & log2FoldChange >= 2)
over_deg_names<-rownames(over_deg_all_info)
##write.csv(over_deg_all_info, "Sci_body-parts_over_DEGs_osculum_vs_bodywall_p001-L2FC2_all_info.csv")
##write.table ( over_deg_names, "Sci_body-parts_over_DEGs_osculum_vs_bodywall_p001-L2FC2.csv", sep=",",  col.names=F, row.names = F)

#Get names and info of DEGs significantly underexpressed in "treatment" (here:osculum)
under_deg_all_info<-subset(DE_MatrixSeq_Results, padj<0.01 & log2FoldChange <= -2)
under_deg_names<-rownames(under_deg_all_info)
##write.csv(under_deg_all_info, "Sci_body-parts_under_DEGs_osculum_vs_bodywall_p001-L2FC2_all_info.csv")
##write.table (under_deg_names, "Sci_body-parts_under_DEGs_osculum_vs_bodywall_p001-L2FC2.csv", sep=",",  col.names=F, row.names = F)

#plot complex-heatmap of Z-scores for comparison of changes in expression of genes per sample for biomin genes

coldata <- data.frame(row.names = colnames(DE_MatrixSeq), DE_Matrix_INFO$condition)#
res <- results(DE_MatrixSeq, contrast= c ("condition", "osculum", "body_wall"))
df<- as.data.frame(res)

#genes of interest: Rename biomin-genes (from biomingenename_map)

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

# Subset 'df' using the specific gene identifiers in biomingenename_map
biomingenes <- df[biomingenename_map$V2, ]

# Use the rownames of 'biomingenes' to subset the assay data
mat <- assay(DE_MatrixSeq_rlog)[rownames(biomingenes), rownames(coldata)]

base_mean <- rowMeans(mat)
# Apply z-score normalization to each row (gene) in 'mat' and transpose the result
mat.scaled <- t(apply(mat, 1, scale)) 
colnames (mat.scaled)<- colnames(mat)

#get log2fold-change for the genes of interest
l2_val <- as.matrix(biomingenes$log2FoldChange) 
colnames(l2_val) <- "logFC"
#get rows from df of for the genes of interest
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
h1 <- Heatmap(mat.scaled, cluster_rows=F,column_labels=colnames(mat.scaled), name="Z-score",cluster_columns=T)
h2 <- Heatmap(l2_val, row_labels = biomingenes$symbol, cluster_rows =F, name="logFC", top_annotation= ha, col =  col_logFC, cell_fun = function(j, i, x, y, w, h, col){grid.text(round(l2_val[i,j],2),x,y)})

h3 <- Heatmap(mean, row_labels =biomingenes$symbol, cluster_rows =F, name = "AveExpr", col= col_AveExpr, cell_fun =function(j,i,x,y,w,h,col){grid.text(round(mean[i,j],2),x,y)})

# Combine the three heatmaps into one and print it
h <- h1+h2+h3
h

