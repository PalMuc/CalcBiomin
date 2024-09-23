library(DESeq2)
library(genefilter)
library(WGCNA)
library(topGO)
library(ggplot2)
library(reshape2)
library (svglite)
#define the working directory. In the repository folder, the relative path should be ".R_scripts/wcgna")
setwd("/Users/ovoigt/Dropbox/_Calcarins/repository/R_scripts/wcgna/")
#get normalized_counts from DESeq2
input_counts<-"../../inputfiles/count_data/gene_counts_combined/Sci_gene_counts_combined.tsv"
raw_counts <-read.csv(input_counts, head=T, sep="\t")
# Keep only the rows in raw_counts where at least 10 values are greater than 10. 
# The threshold of 10 was chosen since there are min. 13 samples in the conditions 
# (low_spic_formation vs. high_spic_formation) to compare.
# This allows for some samples to have lower counts while still retaining the row.
raw_counts <- raw_counts[rowSums(raw_counts > 10) >= 10, ]

#get info of conditions from library-info csv file
infofile <-"../../inputfiles/count_data/counts_info/info_Sci_body_parts+regenI+II.csv"# Read the info csv file
info <- read.csv(infofile, header=T, sep="\t")
#define condition
condition <- factor(info$condition)
#creating dataframe
coldata <- data.frame(row.names = colnames(raw_counts), condition)
dds <- DESeqDataSetFromMatrix(countData= raw_counts, colData = coldata, design= ~condition)

#use vst to normalize
vst <-vst(dds,blind=T)
counts <- assay (vst) [rownames(dds), rownames(coldata)]
allowWGCNAThreads(8)
t_counts<-t(counts)
##


# modified from: https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
# 11.09.2018
gsg<-goodSamplesGenes(t_counts, verbose = 3)
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(t_counts)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(t_counts)[!gsg$goodSamples], collapse=", ")))
  t_counts <- t_counts[gsg$goodSamples, gsg$goodGenes]
}else{
  
  print("all OK")
}
powers<-seq(from =1, to=35, by=1) #choosing a set of soft-thresholding powers
sft<-pickSoftThreshold(t_counts, powerVector=powers, verbose =5, networkType="unsigned", blockSize = 10000) #call network topology analysis function

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.9, col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")

#build a adjacency "correlation" matrix
softPower<-22
adjacency<-adjacency(t_counts, power = softPower)
TOM<-TOMsimilarity(adjacency)

# Call the hierarchical clustering function
geneTree<-hclust(as.dist(1-TOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize<-30
cutHeigthTree<-0.95

# Module identification using dynamic tree cut:
dynamicMods<-cutreeDynamic(dendro = geneTree, distM = 1-TOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, cutHeight = cutHeigthTree);

# Convert numeric lables into colors
dynamicColors<-labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList<-moduleEigengenes(t_counts, colors = dynamicColors)
MEs<-MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss<-1-cor(MEs);

# Cluster module eigengenes
METree<-hclust(as.dist(MEDiss), method = "average");

# Plot the result
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "", cex=0.5)

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge<-mergeCloseModules(t_counts, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors<-merge$colors;

# Eigengenes of the new merged modules:
mergedMEs<-merge$newMEs;

# Rename to moduleColors
moduleColors <- mergedColors
table(moduleColors)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

MEs <- mergedMEs;

MEs_withClass<-MEs
MEs_withClass$Group<-condition

# permT for null t-statistic distribution via permutations (n=1000) across MEs columns, and compares with actual t-stats to assess WGCNA module significance.
permT<-function(x, n){
  tvector<-c()

  for(i in 1:n){
    tvector<-c(tvector, t.test(sample(x)~as.factor(condition))$statistic)
  }
  
  return(tvector)
  
}

##permutations
permTMatrix<-as.data.frame(apply(MEs, 2, permT, n=1000))
rownames(permTMatrix)<-NULL
##unpermuted t
unPermTs<-as.data.frame(t(apply(MEs,2, function(x) t.test(x~as.factor(condition))$statistic)))

par(mfrow=c(1,1))
#Creating histograms 
for(i in 1:ncol(permTMatrix)){
  print(i)
  print(colnames(permTMatrix)[i])
  print(colnames(unPermTs)[i])
  hist(permTMatrix[,i], main = colnames(permTMatrix)[i])
  abline(v=quantile(permTMatrix[,i], probs=c(0.025,0.975)), col="red", lwd=2, lty=2)
  abline(v=unPermTs[i], col="blue", lwd=3)
}

melted_MEs_withClass<-melt(MEs_withClass, id.vars=ncol(MEs_withClass))
colnames(melted_MEs_withClass)<-c("Group","Module","ME")

significantModules<-c("MEblack","MEmidnightblue","MEpink", "MEdarkgrey")

moduleColorsOfInterest<-c("black","midnightblue","pink", "darkgrey")

#plot all modules
ggplot(melted_MEs_withClass, aes(Group, ME)) + geom_boxplot() + facet_wrap(vars(Module)) + geom_jitter()
ggsave("plots/boxplots_all_metamodules.svg", device = "svg", width = 11, height = 8.5, units = "in")

#plot only significant modules
ggplot(subset(melted_MEs_withClass, Module %in% significantModules), aes(Group, ME)) + geom_boxplot() + facet_wrap(vars(Module)) + geom_jitter() + theme_bw()
ggsave("plots/boxplots_significant_metamodules.svg", device = "svg", width = 11, height = 8.5, units = "in")


#############################################
#20220412: export gene names per metamodule
##############################################################

for(colorOfModule in moduleColorsOfInterest){
#for(colorOfModule in significantModules){  
  gene_names<-colnames(t_counts)[mergedColors == colorOfModule]
  write.table(gene_names, paste(colorOfModule,".gene.list", sep=""), row.names = F, col.names = F,quote = F)
  
}

