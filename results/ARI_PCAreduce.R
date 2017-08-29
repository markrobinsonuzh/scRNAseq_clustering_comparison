#####################################
# adjusted Rand Index for PCAreduce
#####################################

library(MLmetrics)
library(caret)
library(mclust)

# define method
METHOD <- "pcaReduce"

# file paths to the clustering results
DATA_DIR <-  "~/Desktop/masterthesis/results/PCAreduce/"

### files with the cell labels, "ground truth":
files_labels <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_labels_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_labels_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_labels_xue2013.txt"))
)

# read in labels
labels <- vector("list", length(files_labels))
names(labels) <- names(files_labels) 

for (i in 1:length(files_labels)) {
  lab_i <-  read.csv(files_labels[[i]], sep="\t")
  labels[[i]] <- as.vector(unlist(lab_i))
}

# read in cluster results 

files_clusters <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_clus_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_clus_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_clus_xue2013.txt"))
)


cluster <- vector("list", length(files_clusters))
names(cluster) <- names(files_clusters) 

for (i in 1:length(files_clusters)) {
  clus_i <-  read.csv(files_clusters[[i]], sep="\t")
  cluster[[i]] <- as.vector(unlist(clus_i))
}


## compute adjusted Rand Index
stat.cluster<- vector("list", length(files_labels))
names(stat.cluster) <- names(files_labels) 



for (i in 1:length(files_labels)) {
  
  pred_i   <- cluster[[i]]
  lab_i    <- labels[[i]]
  stat.cluster[[i]]$ARI <- adjustedRandIndex(pred_i,lab_i)
  
}

# store result file

file_results <- paste0(DATA_DIR,METHOD ,"_ARI_", names(files_labels), ".txt")

for (i in 1:length(file_results)) {
  res_i <- stat.cluster[[i]]
  write.table(res_i, file = file_results[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

