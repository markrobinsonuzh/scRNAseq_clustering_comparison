#######################################################
### Compute the ARI for the single cluster computations
### For all Methods and datasets
######################################################

#load libraries
library(dplyr)
library(mclust)
# load data files
DATASET <- c("koh2016", "kumar2015",  "trapnell2014","zhengmix2016", "simDataKumar")
datatype<-  "unfiltered"

for (h in seq_len(length(DATASET)) ) {
  
  data_files <- as.list(paste0("results/run_results/cluster_single_",datatype,"_",DATASET[h],".rda"))
  
  load(data_files[[1]])
  
  # reformat labels as factors
  
  clusters$labels <- clusters$labels%>%as.factor
  
  # create lookupfile
  levels(clusters$labels) <- c(1:length(levels(clusters$labels)))
  clusters$labels <- as.integer(clusters$labels)%>%as.data.frame
  ### calculate ARI scores
  res.ari <- vector("list", length(clusters))
  names(res.ari) <- names(clusters)
  
  for (i in names(clusters)) {
    lab_i    <- clusters$labels[,1]
    pred_i   <- unlist( clusters[[i]])
    res.ari[[i]]$ARI <- tryCatch(adjustedRandIndex(pred_i,lab_i), error=function(e) NA )
    
  }
  save(res.ari, file = paste0("results/run_results/ari_single_",datatype,"_",DATASET[h],".rda"))
  
  print(res.ari)
  
}

