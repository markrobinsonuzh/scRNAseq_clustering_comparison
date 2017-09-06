# This is demonstration of using SINCERA to analysis human IPF and normal epithelial single cell RNA-seq data (Xu et al., JCI Insight 2016)
# Author: Minzhe Guo (minzhe.guo@cchmc.org)
library(SINCERA)

# load IPF data, which contains two data frames
# - the ipf.exprmatrix contains the TPM values of 540 single cell transcriptomes
# - the ipf.cells contains the cell sample and condition information
data(IPF)

ls()

dim(ipf.exprmatrix)

head(ipf.cells)

table(ipf.cells$Diagnosis) # normal or ipf cells

table(ipf.cells$DataID) # cell sample


# The analysis starts with running the construct function to create an R S4 object, which will hold all the data and analysis results.
# The function takes expression matrix and cell sample information as input
sc <- construct(exprmatrix=ipf.exprmatrix, samplevector=ipf.cells$DataID)

# After contruction, you can use setCellMeta to add more cell information to sincera
# such as adding the condition information (CONTROL or IPF) of cells into sincera object
sc <- setCellMeta(sc, name="CONDITION", value=ipf.cells$Diagnosis)

# use getCellMeta functoin to assess a specific meta data
# In most of the SINCERA functions, cell grouping will be based on the GROUP meta data
# The GROUP meta was initialized to sample information during the construction
table(getCellMeta(sc, name="GROUP"))

# Identify and remove low quality cells.
# The key parameters of running this function include: “min.expression”,
# which specifies the minimum expression value for a gene to be considered as expressed,
# and “min.genes”, which specifies the lower bound of the number of expressed genes in a cell.
sc <- filterLowQualityCells(sc, min.expression=1, min.genes=1000, do.plot = T)

# Set the minimum value of expression to 0.01
sc <- expr.minimum(sc, value=0.01)

# Run batch.analysis function to generate plots that may help identify potential batch differences.
sc <- batch.analysis(sc, analysis=c("distribution"), min.expression=1)



# Filter out non-expressed genes
sc <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=1, min.samples=1)

# Perform z-score scaling
sc <- normalization.zscore(sc, pergroup=FALSE)

# do PCA using all genes
sc <- doPCA(sc, genes=NULL, use.fast = T)

# plot the standard deviation of PCA components
plotPCASD(sc, num.pcs = 20)

# do tSNE using PCA components
sc <- doTSNE(sc, genes=NULL, dims = 1:5, use.fast = T)

# plot cells in tSNE spaces
plotRDS(sc, feature.type="tsne")




## first iteration


#In the first iteration, we selected genes for clustering using the following criteria
# genes expressed in at least 10 cells in at least 6 samples
min.samples <- 6
obj <- prefilterGenes(sc, pergroup=TRUE, min.expression=1, min.cells=10, min.samples=min.samples)
# genes with at least 0.7 specificity in at least 6 samples
obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7)

# set the selected genes for clustering
sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))


# use gap statistics to determine the number of clusters
if (FALSE) {
  x <- getExpression(sc, scaled=T, genes=getGenesForClustering(sc))
  
  library(cluster)
  cordist <- function(y) as.dist((1-cor(t(y)))/2)
  hclustForGap <- function(y, k) list(cluster=cutree(hclust(cordist(y), method = "average"),k=k))
  gapstats <- clusGap(t(x),
                      FUN = hclustForGap,
                      K.max = 10,
                      B = 100)
}

library(cluster)
print(ipf.gapstats$iter1) # the gap statistics suggested 5 clusters

# use the default HC algorithm to identify 5 clusters
# the clustering results were saved to the "CLUSTER" metadata
# if update.cellgroup is TRUE, the GROUP meta data will also be updated
sc <- cluster.assignment(sc, k=5)


# different display of the dendrogram tree
plotHC(sc, horiz=T, show.labels=T, do.radial=FALSE)
plotHC(sc, do.radial=T)
plotHC(sc)

# The condition information of cells in different clusters
print(table(getCellMeta(sc, "CONDITION"), getCellMeta(sc, "CLUSTER")))

# Plot marker expression
plotMarkers(sc, genes=c("SFTPB","SFTPC","SLC34A2","NAPSA"))
plotMarkers(sc, genes=c("KRT5","KRT14","TP63","ITGB4"))
plotMarkers(sc, genes=c("MUC5AC","MUC5B","SPDEF","MUC16"))


# Cluster 1 and 4 are in the same major branch in the dendrogram tree
# And they are likely subclusters of control cells.
# Since in this study, we are less interested in the subcluster structure of the control cells
# we merged the two clusters
cluster.1 <- getCellMeta(sc, "CLUSTER")
cluster.1[which(cluster.1==4)]=1
sc <- setCellMeta(sc, name="CLUSTER.1", value=cluster.1)
sc <- copyCellMeta(sc, "CLUSTER.1", "GROUP")

# tSNE Plot the clustering results
plotRDS(sc, feature.type="tsne")

# predicted cluster specific signature genes based on the following criteria
# <0.05 pvalue of one-tailed welch's t test
# >=1 average expression in the cluster
# expressed in at least 50% of the cluster cells
sc <- signature.prediction(sc, use.logireg = FALSE, diff.method = "welch", use.fdr = FALSE, diff.thresh = 0.05,
                           avg.thresh=1, min.expression = 1, pct.thresh = 0.5, fc.thresh = -Inf)

# get the signature genes as a data frame
# Use those genes for functional enrichment
siggenes <- getSigGenes(sc)
head(siggenes)

# map cell clusters to cell types
sc <- setCellType(sc, do.reset = T)
sc <- setCellType(sc, groups=c("1", "2", "3","5"), types=c("AT2", "Basal", "Indeterminate","Goblet"))
getCellType(sc)

# add cell type marker information to sincera
markers <- data.frame(TYPE=c(rep("AT2", 4), rep("Basal", 3), rep("Goblet", 4)),
                      SYMBOL=c("SFTPB","ABCA3","SLC34A2","LPCAT1","KRT5","KRT14","TP63","SPDEF","MUC5AC","MUC5B","SCGB3A2"))
sc <- setCellTypeMarkers(sc, markers)

head(getCellTypeMarkers(sc))

# perform rank-aggregation based cell type validation using marker expression
sc <- celltype.validation(sc)


# second iteration

# In the second iteration, We used the signature genes of the clusters identified in the first iteration to performclustering
siggenes <- as.character(unique(getSigGenes(sc)$SYMBOL))

sc <- setGenesForClustering(sc, value=siggenes)

# use gap statistics to determine the number of clusters
print(ipf.gapstats$iter2)

# use the default algorithm to identify 4 clusters
sc <- cluster.assignment(sc, k=4)

# tSNE plot of the new clusters
plotRDS(sc, feature.type="tsne")

# plot marker patterns based on the new clustering results
plotMarkers(sc, genes=c("SFTPB","SFTPC","SLC34A2","NAPSA"))
plotMarkers(sc, genes=c("KRT5","KRT14","TP63","ITGB4"))
plotMarkers(sc, genes=c("MUC5AC","MUC5B","SPDEF","MUC16"))

# Predict signature genes based on the new cluster assignment
if (FALSE) {
  sc <- signature.prediction(sc, use.logireg = FALSE, diff.method = "welch", use.fdr = FALSE,
                             diff.thresh = 0.05, avg.thresh=1, min.expression = 1, pct.thresh = 0.5, fc.thresh = -Inf)
  siggenes <- getSigGenes(sc)
}
siggenes <- sigs$iter2

table(siggenes$GROUP)

# map the second iteration clusters to cell types
sc <- setCellType(sc, do.reset = T)
sc <- setCellType(sc, groups=c("1", "2", "3","4"), types=c("AT2", "Basal", "Indeterminate","Goblet"))
getCellType(sc)

# rank aggregation based cell type validation
sc <- celltype.validation(sc)

# save the 2nd iteration clustering results, since the 3rd round clustering will override the "CLUSTER" metadata
sc <- copyCellMeta(sc, from="CLUSTER", to="CLUSTER.2")


# third iteration

# used the signature genes from the 2nd iteration to perform clustering in the 3rd iteration
siggenes <- as.character(unique(siggenes$SYMBOL))

print(length(siggenes))

sc <- setGenesForClustering(sc, value=siggenes)

# gap statistics suggested 3 clusters
print(ipf.gapstats$iter3)

# use default algorithm to assign cells to 3 clusters
sc <- cluster.assignment(sc, k=3)

# plot marker expression
plotMarkers(sc, genes=c("SFTPB","SFTPC","SLC34A2","NAPSA"))
plotMarkers(sc, genes=c("KRT5","KRT14","TP63","ITGB4"))
plotMarkers(sc, genes=c("MUC5AC","MUC5B","SPDEF","MUC16"))

# tSNE plot of the 3rd iteration clustering
plotRDS(sc, feature.type="tsne")

# map 3rd iteration clusters to cell types
sc <- setCellType(sc, do.reset = T)
sc <- setCellType(sc, groups=c("1", "2", "3"), types=c("AT2","Goblet", "Basal"))
getCellType(sc)

# rank aggregation based cell type validation
sc <- celltype.validation(sc)

# save the 3rd iteration clustering results
sc <- copyCellMeta(sc, from="CLUSTER", to="CLUSTER.3")

# since the 3rd iteration does not improve the performance, we stopped iteration and the 2nd iteration result as the cell clustering result
sc <- copyCellMeta(sc, from="CLUSTER.2", to="GROUP")

genes <- c("SFTPC","SFTPB","ABCA3","SLC34A2","LPCAT1","HOPX", "AGER","PDPN","TP63","KRT5", "KRT14","NGFR","IL10","PLA2G4C","SLC26A8","AKT2","COL1A1","LTBP4","MUC5AC","MUC5B","SPDEF","SCGB1A1")
plotHeatmap(sc, genes=genes, show.labRow = T)


# use the signature genes of the 2nd iteration clusters to perform PCA and tSNE
print(length(siggenes))
sc <- doPCA(sc, genes=siggenes, use.fast = T)
plotPCASD(sc, num.pcs = 20)
sc <- doTSNE(sc, genes=NULL, dims = 1:5, use.fast = T)
plotRDS(sc, feature.type="tsne")



# clustering comparison

sc <- copyCellMeta(sc, "CLUSTER.2", "ihc")


sc <- cluster.assignment(sc, feature.type = "pca", rds.dims=1:5, clustering.method = "graph")
sc <- copyCellMeta(sc, "CLUSTER","pca.graph")



sc <- cluster.assignment(sc, feature.type = "tsne", rds.dims=1:2, clustering.method = "pam", k=4)
sc <- copyCellMeta(sc, "CLUSTER","tsne.pam")



m <- pData(sc@data)[, c("pca.graph","tsne.pam", "ihc")]

mt <- clusterings.compare(m)



viz <- data.frame(getTSNE(sc, name="rds"), Type=factor(getCellMeta(sc, name="ihc")))
viz$Instability=factor(round(mt[rownames(viz), "Instability"],2))
g <- ggplot(viz, aes(x=tSNE1, y=tSNE2, shape=Type, col=Instability)) + geom_point(size=2)
g <- g + scale_color_manual(values=c("grey","green","orange","blue"))
g <- g + ggtitle("tSNE plot of cell clustering assignment instability")
g <- g + theme_bw() + theme(panel.grid = element_blank())
print(g)



## save all analyses and results for reproduction

save(sc, file="ipf-analyzed.rda")