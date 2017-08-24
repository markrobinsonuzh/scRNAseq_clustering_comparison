####################################################
# PLots for Data set KUMAR 2015 GSE60749-GPL13112
####################################################
library(gridExtra)
library(ggfortify)

library(Seurat)
library(scater)
library(plyr)
# load the data
##############################
DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds")
)
METHOD_NAME <- as.character(c("tSNEkmeans",
                 "SNNCliq",
                 "SIMLR",
                 "Seurat",
                 "SC3",
                 "pcaReduce"))
# define method name
method <- list(
  tSNEkmeans = NULL,
  SNNCliq = NULL,
  SIMLR = NULL,
  Seurat= NULL,
  SC3 = NULL,
  pcaReduce = NULL
  
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

# create dir
RES_DIR <- "~/Desktop/masterthesis/results"

fileslabels <- list(
  tSNEkmeans = file.path(RES_DIR, "tSNEkmeans/tSNEkmeans_labels_kumar2015.txt"),
  SNNCliq = file.path(RES_DIR, "SNNCliq/SNNCliq_labels_kumar2015.txt"),
  SIMLR = file.path(RES_DIR, "SIMLR/SIMLR_labels_kumar2015.txt"),
  Seurat= file.path(RES_DIR, "Seurat/Seurat_labels_kumar2015.txt"),
  SC3 = file.path(RES_DIR, "Seurat/Seurat_labels_kumar2015.txt"),
  pcaReduce = file.path(RES_DIR, "Seurat/Seurat_labels_kumar2015.txt")
)

# load cell labels
labels<- vector("list", length(fileslabels))

names(labels) <-  names(fileslabels)

for(i in 1:length(labels)) {
  labels[[i]] <- read.csv(fileslabels[[i]], sep = "\t")
}
# rename cell levels

labels <- unlist(labels[[1]])
labels <- revalue(labels, c("Dgcr8 knockout mouse embryonic stem cells.culture conditions: serum+LIF"="Dgcr8 knockout mouse serum+LIF", "v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF"="v6.5 mouse 2i+LIF","v6.5 mouse embryonic stem cells.culture conditions: serum+LIF"="v6.5 mouse serum+LIF"))
Labels <- as.factor(labels)

# load the clustering results from method:
filesclusters <- list(
  tSNEkmeans = file.path(RES_DIR, "tSNEkmeans/tSNEkmeans_clus_kumar2015.txt"),
  SNNCliq = file.path(RES_DIR, "SNNCliq/SNNCliq_clus_kumar2015.txt"),
  SIMLR = file.path(RES_DIR, "SIMLR/SIMLR_clus_kumar2015.txt"),
  Seurat= file.path(RES_DIR, "Seurat/Seurat_clus_kumar2015.txt"),
  SC3 = file.path(RES_DIR, "SC3/sc3_clus_kumar2015.txt"),
  pcaReduce = file.path(RES_DIR, "PCAreduce/PCAreduce_clus_kumar2015.txt")
)


# load clusters
clusters<- vector("list", length(filesclusters))

names(clusters) <-  names(filesclusters)

for(i in 1:length(clusters)) {
  clusters[[i]] <- as.factor(unlist(read.csv(filesclusters[[i]], sep = "\t")))
}


# Plot on PC1 and PC2

plot.method <- list(
  tSNEkmeans = NULL,
  SNNCliq = NULL,
  SIMLR = NULL,
  Seurat= NULL,
  SC3 = NULL,
  pcaReduce = NULL
)



# PCA dim reduce on log2 transformed, normalized count_lstpm
data[[1]] <- plotPCA(object=data[[1]],  exprs_values="exprs" ,return_SCESet=TRUE, scale_features=TRUE)
# extract components
pc.data <-as.data.frame(data[[1]]@reducedDimension)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

vec <- c(1:6)

for (i in 1:length(clusters)){
  
  plot.method[[i]] <- ggplot(data = pc.data , mapping = aes(x=PC1,y=PC2, group=Labels, shape=Labels))+geom_point(aes_string(color=clusters[[i]]))+scale_colour_manual(values=cbbPalette)+labs(colour=METHOD_NAME[i])

}

grid.arrange(plot.method[[1]],plot.method[[2]],plot.method[[3]],plot.method[[4]],plot.method[[5]],plot.method[[6]], ncol=2)

