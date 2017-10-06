####################################################
# PLots for Data set TRAPNELL 2014 GSE52529-GPL16791
# PCA plots by library total count and cluster
####################################################
# load libraries
library(gridExtra)
library(ggfortify)

library(scater)
library(plyr)
library(cowplot)
# load the data
##############################
DATA_DIR <- "data"
files <- list(
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda")
)
METHOD_NAME <- as.character(c("tSNEkmeans",
                              "SNNCliq",
                              "SIMLR",
                              "Seurat",
                              "SC3",
                              "pcaReduce",
                              "dbscan", 
                              "cidr"))


# load data sets

data <- labels<- vector("list", length(files))

names(data) <-names(labels) <-  names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load the the labels:

# create filenames
RES_DIR <- "results"

fileslabels <- list(
  tSNEkmeans = file.path(RES_DIR, "tSNEkmeans/tSNEkmeans_labels_trapnell2014.txt"),
  SNNCliq = file.path(RES_DIR, "SNNCliq/SNNCliq_labels_trapnell2014.txt"),
  SIMLR = file.path(RES_DIR, "SIMLR/SIMLR_labels_trapnell2014.txt"),
  Seurat= file.path(RES_DIR, "Seurat/Seurat_labels_trapnell2014.txt"),
  SC3 = file.path(RES_DIR, "Seurat/Seurat_labels_trapnell2014.txt"),
  pcaReduce = file.path(RES_DIR, "Seurat/Seurat_labels_trapnell2014.txt"),
  dbscan = file.path(RES_DIR, "dbscan/dbscan_labels_trapnell2014.txt"),
  cidr = file.path(RES_DIR, "cidr/cidr_labels_trapnell2014.txt")
  
)

# load cell labels
labels<- vector("list", length(fileslabels))

names(labels) <-  names(fileslabels)

for(i in 1:length(labels)) {
  labels[[i]] <- read.csv(fileslabels[[i]], sep = "\t")
}
Labels <- unlist(labels[[1]])
Labels <- as.factor(Labels)
# load the clustering results from methods:
filesclusters <- list(
  tSNEkmeans = file.path(RES_DIR, "tSNEkmeans/tSNEkmeans_clus_trapnell2014.txt"),
  SNNCliq = file.path(RES_DIR, "SNNCliq/SNNCliq_clus_trapnell2014.txt"),
  SIMLR = file.path(RES_DIR, "SIMLR/SIMLR_clus_trapnell2014.txt"),
  Seurat= file.path(RES_DIR, "Seurat/Seurat_clus_trapnell2014.txt"),
  SC3 = file.path(RES_DIR, "SC3/sc3_clus_trapnell2014.txt"),
  pcaReduce = file.path(RES_DIR, "PCAreduce/PCAreduce_clus_trapnell2014.txt"),
  dbscan = file.path(RES_DIR, "dbscan/dbscan_clus_trapnell2014.txt"),
  cidr = file.path(RES_DIR, "cidr/cidr_clus_trapnell2014.txt")
  
)


# load clusters
clusters<- vector("list", length(filesclusters))

names(clusters) <-  names(filesclusters)

for(i in 1:length(clusters)) {
  clusters[[i]] <- as.factor(unlist(read.csv(filesclusters[[i]], sep = "\t")))
}

# Plot on PC1 and PC2

plot.method <- NULL

# PCA dim reduce on log2 transformed, normalized count_lstpm
data[[1]] <- plotPCA(object=data[[1]],  exprs_values="exprs" ,return_SCESet=TRUE, scale_features=TRUE)
# extract components
pc.data <-as.data.frame(data[[1]]@reducedDimension)

# plot it
# color grouping by totalcount of library
total.counts <- pData(data[[1]])$total_count
# mid point of color range defined by median
mid <- median(pData(data[[1]])$total_count)

for (i in 1:length(clusters)){
  
  plot.method[[i]] <- ggplot( data = pc.data, aes(x=PC1,y=PC2, color=total.counts) )+ 
    geom_point(aes_string(shape=clusters[[i]]))+
    labs(colour=METHOD_NAME[i])+
    scale_color_gradient2(midpoint = mid, low="blue", mid = "yellow",high = "red")+
    scale_shape_manual(values=c(1:10), guide=FALSE)
  
}


plot2by3 <- plot_grid(plotlist=plot.method, labels = "auto", nrow = 3,ncol =3)
save_plot("results/plots/plot_clusterlibsize_trapnell2014.pdf", plot2by3, base_height = 10, base_width = 15)


## appendix
mid <- median(pData(data[[1]])$total_count)
ggplot(data = pc.data, aes(x=PC1,y=PC2, color=pData(data[[1]])$total_count, shape = clusters[[1]])) + geom_point() + scale_color_gradient2(midpoint = mid, low="blue", mid = "brown",high = "red")




