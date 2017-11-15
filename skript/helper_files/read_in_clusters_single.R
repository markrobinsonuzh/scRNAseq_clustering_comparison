#################################
### create results file      ####
#################################

# This file reads in the text files from the cluster runs for each method as well as the labels with the ground truth for the respective dataset. 
# The cluster and labels are then saved as a single list per dataset and stored as a "datasetname".rdata file in the run_results directory


#load libraries
library(plyr)
library(dplyr)


source("skript/helper_files/Helper_functions.R")



# define method "pcaReduce","dbscan", "RtSNEkmeans", "SC3", "SIMLR", "SNNCliq", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm"
METHOD <- c("pcaReduce","dbscan", "RtSNEkmeans", "SC3", "SIMLR", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm")   

# file paths to the clustering results
DATA_DIR <-  "results"
DATASET <-"zhengmix2016"   # "kumar2015" ,"trapnell2014" ,"zhengmix2016" , "koh2016" , "simDataKumar"



files_labels <- file.path(DATA_DIR, METHOD,paste0(METHOD,"_labels_",DATASET,".txt"))%>%as.list()
names(files_labels) <- METHOD
  
labels <- vector("list", length(files_labels))
names(labels) <- names(files_labels) 
  
for (i in names(labels)) {
    lab_i <-  read.csv(files_labels[[i]], sep="")
    labels[[i]] <- as.vector(unlist(lab_i))
    
  }

# read in cluster results

files_cluster <- file.path(DATA_DIR, METHOD,paste0(METHOD,"_clus_",DATASET,".txt"))%>%as.list() # file path
names(files_cluster) <- METHOD # gives names

clusters <- vector("list", length(files_cluster))
names(clusters) <- names(files_cluster) 

for (i in names(clusters)) {
  lab_i <-  read.csv(files_cluster[[i]], sep="")
  clusters[[i]] <- as.vector(unlist(lab_i))
  
}

#### generate object with cluster and files
clusters$labels <- labels$pcaReduce
### Save object
save(clusters,file = paste0("results/run_results/cluster_single_",DATASET,".rda"))


