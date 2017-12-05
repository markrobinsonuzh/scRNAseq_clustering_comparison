#######################################################
### Compute the ARI for a  range of clusters k 
### For all Methods and datasets
######################################################
#This script takes the cluster results from the k range run (cluster_krange_",DATASET[h],".rda") and computes the ARI. 
# The results are stored as a ari_krange_DATASET.rda file.

# run it per h , h for loop not working
#load libraries
library(dplyr)
library(mclust)
# load data files
DATASET <- c( "kumar2015", "trapnell2014","koh2016" , "simDataKumar", "zhengmix2016")
for (h in seq_len(length(DATASET) )) {
  print(DATASET[h])
  data_files <- as.list(paste0("results/run_results/cluster_krange_",DATASET[h],".rda"))
  
  load(data_files[[1]])
  ### calculate ARI scores
  res.ari <- vector("list", length(clusters))
  names(res.ari) <- names(clusters)
  
  for (i in names(clusters)) {
    
    lab_i    <- as.integer(clusters[[i]]$labels )
    clusters[[i]] <- subset( clusters[[i]], select = -c(labels) )
    for (j in seq_len( ncol( clusters[[i]] ) )) {
      res.ari[[i]]$par <-  colnames(clusters[[i]])
      res.ari[[i]]$ncluster <- apply(clusters[[i]], 2, function(x) length(unique( x )))
      pred_i   <- unlist( clusters[[i]][j] )
      res.ari[[i]]$ARI[j] <- tryCatch(adjustedRandIndex(pred_i,lab_i),  error=function(e) NA)
    }
    
  }
  
  save(res.ari, file = paste0("results/run_results/ari_krange_",DATASET[h],".rda"))
  
}


#### Appendix
