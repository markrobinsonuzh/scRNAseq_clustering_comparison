#################################
### store cluster results under single object     ####
#################################

# This file reads in the text files from the cluster runs for each method as well as the labels with the ground truth for the respective dataset. 
# The cluster and labels are then saved as a single list per dataset and stored as a "datasetname".rdata file in the run_results directory
# To do: if method not available, load labels from general file...

# load the helper function
source("skript/helper_files/Helper_functions.R")
# define which method to load
# methods: "pcaReduce","dbscan", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "SNNCliq", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm"
METHOD <- c("pcaReduce", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm" )   

# file paths to the clustering results, change the path according to the processed datasets
DATA_DIR <-  "results"
# datasets: "default", "filtered","unfiltered","optimalk"
datatype <- "optimalk"
#
DATASET <-c("kumar2015" ,"trapnell2014" ,"zhengmix2016" , "koh2016" , "simDataKumar", "simDataKumar2")   
# store .rda objects , per dataset
for (i in seq_len(length(DATASET)) ) {
save_cluster_single3(  METHOD,DATA_DIR, DATASET[i], datatype )
}

