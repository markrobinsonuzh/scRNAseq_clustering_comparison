#################################
### Compute the F1 scores    ####
#################################
#load libraries
library(dplyr)
library(plyr)

#load helper file
source("/Users/angeloduo/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/helper_calc_f1_scores.R")
source("/Users/angeloduo/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/Helper_functions.R")

# define method
METHOD <- "pcaReduce"

# file paths to the clustering results
DATA_DIR <-  "results/PCAreduce/"

### files with the cell labels, "ground truth":
files_labels <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_labels_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_labels_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_labels_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_labels_koh2016.txt"))
)

# read in labels
labels <- read.labels(files_labels = files_labels)%>%sapply(as.factor)
# create lookupfile
look <- sapply(labels, function(x) levels(x) <- c(1:length(levels(x))) )
# change levels to integer
for (i in seq_len(length(labels))) {
  {
  levels(labels[[i]]) <- look[[i]] 
  }
  lbls <- sapply(labels,as.integer)
}

# read in cluster results, format as integer

files_clusters <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_clus_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_clus_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_clus_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_clus_koh2016.txt"))
)
clus <- read.cluster(files_clusters=files_clusters)%>%sapply(as.integer)

### calculate F1 scores
for (i in seq_len(length(clus))) {
         print(calc_f1_scores(lbls[[i]],clus[[i]])  ) 
}

#####


