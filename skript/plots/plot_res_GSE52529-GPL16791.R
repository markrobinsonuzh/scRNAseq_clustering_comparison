####################################################
# PLots for Data set TRAPNELL 2014 GSE52529-GPL16791
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
                              "dbscan"))
# define method name
method <- list(
  tSNEkmeans = NULL,
  SNNCliq = NULL,
  SIMLR = NULL,
  Seurat= NULL,
  SC3 = NULL,
  pcaReduce = NULL,
  dbscan=NULL
  
)

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
  dbscan = file.path(RES_DIR, "dbscan/dbscan_labels_trapnell2014.txt")
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
  dbscan = file.path(RES_DIR, "dbscan/dbscan_clus_trapnell2014.txt")
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

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000001","#000002","#000003","#000005")

vec <- c(1:6)

for (i in 1:length(clusters)){
  
  plot.method[[i]] <- ggplot(data = pc.data , mapping = aes(x=PC1,y=PC2, group=Labels, shape=Labels))+
    geom_point(aes_string(color=clusters[[i]]))+scale_colour_manual(values=cbbPalette)+labs(colour=METHOD_NAME[i])

}

plot2by3 <- plot_grid(plotlist=plot.method, labels = "auto", nrow = 3,ncol =3)
save_plot("results/plots/plot_cluster_trapnell2014.pdf", plot2by3, base_height = 10, base_width = 15)


## appendix



