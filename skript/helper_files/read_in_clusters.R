#################################
### create results file      ####
#################################

#load libraries
library(plyr)
library(dplyr)


source("/Users/angeloduo/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/Helper_functions.R")



# define method
METHOD <- c("pcaReduce","dbscan", "RtSNEkmeans", "SC3", "Seurat", "SIMLR", "SNNCliq")

# file paths to the clustering results
DATA_DIR <-  "results"
DATASET <-   "Koh2016"



files_labels <- file.path(DATA_DIR, METHOD,paste0(METHOD,"_labels_",DATASET,".txt"))%>%as.list()
names(files_labels) <- METHOD
  
labels <- vector("list", length(files_labels))
names(labels) <- names(files_labels) 
  
for (i in names(labels)) {
    lab_i <-  read.csv(files_labels[[i]], sep="")
    labels[[i]] <- as.vector(unlist(lab_i))
    
  }


# read in cluster results

files_cluster <- file.path(DATA_DIR, METHOD,paste0(METHOD,"_clus_",DATASET,".txt"))%>%as.list()
names(files_cluster) <- METHOD

clusters <- vector("list", length(files_cluster))
names(clusters) <- names(files_cluster) 

for (i in names(clusters)) {
  lab_i <-  read.csv(files_cluster[[i]], sep="")
  clusters[[i]] <- as.vector(unlist(lab_i))
  
}

#### generate object with cluster and files
clusters$labels <- labels$pcaReduce
### Save object
save(clusters,file = paste0("results/run_results/",DATASET,".rda"))
