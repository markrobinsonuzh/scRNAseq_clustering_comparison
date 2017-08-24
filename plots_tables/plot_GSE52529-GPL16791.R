####################################################
# PLots for Data set TRAPNELL 2014 GSE52529-GPL16791
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
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds")
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

# files labels
RES_DIR <- "~/Desktop/masterthesis/results"
# define dataset
dataset <- "trapnell2014.txt"
fileslabels <- list(
  tSNEkmeans = file.path(RES_DIR, paste0("tSNEkmeans/tSNEkmeans_labels_",dataset)),
  SNNCliq = file.path(RES_DIR, paste0("SNNCliq/SNNCliq_labels_",dataset)),
  SIMLR = file.path(RES_DIR, paste0("SIMLR/SIMLR_labels_",dataset)),
  Seurat= file.path(RES_DIR, paste0("Seurat/Seurat_labels_",dataset)),
  SC3 = file.path(RES_DIR, paste0("Seurat/Seurat_labels_",dataset)),
  pcaReduce = file.path(RES_DIR, paste0("Seurat/Seurat_labels_",dataset))
)

# load cell labels
labels<- vector("list", length(fileslabels))

names(labels) <-  names(fileslabels)

for(i in 1:length(labels)) {
  labels[[i]] <- read.csv(fileslabels[[i]], sep = "\t")
}
# rename cell levels

labels <- unlist(labels[[1]])
T0=grep("T0",labels, value=TRUE)
length(grep("T24",labels, value=FALSE))
length(grep("T48",labels, value=TRUE))
length(c(grep("T0",labels, value=TRUE)))
labels <- mapvalues(labels, from = c(grep("T0",labels, value=TRUE),grep("T24",labels, value=TRUE),grep("T48",labels, value=TRUE)) , to=c(array("T0", dim=96),array("T24", dim=96),array("T48", dim=96)) )
Labels <- as.factor(labels)
# load the clustering results from method:
filesclusters <- list(
  tSNEkmeans = file.path(RES_DIR, paste0("tSNEkmeans/tSNEkmeans_clus_",dataset)),
  SNNCliq = file.path(RES_DIR, paste0("SNNCliq/SNNCliq_clus_",dataset)),
  SIMLR = file.path(RES_DIR, paste0("SIMLR/SIMLR_clus_",dataset)),
  Seurat= file.path(RES_DIR, paste0("Seurat/Seurat_clus_",dataset)),
  SC3 = file.path(RES_DIR, paste0("SC3/sc3_clus_",dataset)),
  pcaReduce = file.path(RES_DIR, paste0("PCAreduce/PCAreduce_clus_",dataset))
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

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000001","#000002","#000003","#000005")

vec <- c(1:6)

for (i in 1:length(clusters)){
  
  plot.method[[i]] <- ggplot(data = pc.data , mapping = aes(x=PC1,y=PC2, group=Labels, shape=Labels))+geom_point(aes_string(color=clusters[[i]]))+scale_colour_manual(values=cbbPalette)+labs(colour=METHOD_NAME[i])
  
}

grid.arrange(plot.method[[1]],plot.method[[2]],plot.method[[3]],plot.method[[4]],plot.method[[5]],plot.method[[6]], ncol=2)

plot.method[[5]]
