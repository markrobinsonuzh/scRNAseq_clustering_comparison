#######################################################
### Compute the ARI for a  range of clusters k 
### For all Methods and datasets
######################################################

#load libraries
library(dplyr)
library(mclust)
# load data files
DATASET <- c("koh2016", "kumar2015", "trapnell2014", "zhengmix2016")
for (h in seq_len(length(DATASET)) ) {
  
  data_files <- as.list(paste0("results/run_results/cluster_krange_",DATASET[h],".rda"))
  
  load(data_files[[1]])
  # reformat labels as factors
  clusters$labels <- clusters$labels%>%as.factor

  # create lookupfile
  #look <- lapply(clusters$labels, function(x) levels(x) <- c(1:length(levels(x))) )
  levels(clusters$labels) <- c(1:length(levels(clusters$labels)))
  clusters$labels <- as.integer(clusters$labels)%>%as.data.frame
  ### calculate ARI scores
  res.ari <- vector("list", length(clusters))
  names(res.ari) <- names(clusters)
  
  for (i in names(clusters)) {
    lab_i    <- clusters$labels[,1]
    
    for (j in seq_len( ncol( clusters[[i]] ) )) {
    
    pred_i   <- unlist( clusters[[i]][j] )
    res.ari[[i]]$ARI[j] <- adjustedRandIndex(pred_i,lab_i)
    res.ari[[i]]$par <-  colnames(clusters[[i]])
    
  }
    save(res.ari, file = paste0("results/run_results/ari_krange_",DATASET[h],".rda"))
    
    print(res.ari)
  }

  
  }


#### Appendix