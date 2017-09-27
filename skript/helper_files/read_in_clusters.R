#################################
### create results file      ####
#################################

# This file reads in the text files from the cluster runs for each method as well as the labels with the ground truth for the respective dataset. 
# The cluster and labels are then saved as a single list per dataset and stored as a .rdata file.


#load libraries
library(plyr)
library(dplyr)


source("/Users/angeloduo/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/Helper_functions.R")



# define method
METHOD <- c("PCAReduce","dbscan", "RtSNEkmeans", "SC3", "Seurat", "SIMLR", "SNNCliq")

# file paths to the clustering results
DATA_DIR <-  "results"
DATASET <-   "trapnell2014"



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
clusters$labels <- labels$RtSNEkmeans
### Save object
save(clusters,file = paste0("results/run_results/",DATASET,".rda"))

