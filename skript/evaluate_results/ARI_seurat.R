#####################################
# adjusted Rand Index for tSNEkmeans
#####################################

#load libraries
source("skript/helper_files/Helper_functions.R")

library(MLmetrics)
library(caret)
library(mclust)

# define method
METHOD <- "Seurat"

# file paths to the clustering results
DATA_DIR <- "results/Seurat/"

### files with the cell labels, "ground truth":
files_labels <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_labels_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_labels_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_labels_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_labels_koh2016.txt"))
)

# read in labels
labels <- read.labels(files_labels = files_labels)


# read in cluster results 

files_clusters <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_clus_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_clus_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_clus_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_clus_koh2016.txt"))
)

cluster <- read.cluster(files_clusters=files_clusters)


################################
## compute adjusted Rand Index
################################

ari.cluster<- vector("list", length(files_labels))
names(ari.cluster) <- names(files_labels) 



for (i in names(files_labels)) {
  
  pred_i   <- cluster[[i]]
  lab_i    <- labels[[i]]
  ari.cluster[[i]]$ARI <- adjustedRandIndex(pred_i,lab_i)
  
}

# store result file
store.ari(ari.cluster=ari.cluster, DATA_DIR=DATA_DIR, METHOD=METHOD)
