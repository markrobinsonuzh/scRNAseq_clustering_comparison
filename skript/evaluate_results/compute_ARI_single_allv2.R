#######################################################
### Compute the ARI for the single cluster computations
### For all Methods and datasets
######################################################

#load libraries
library(dplyr)
library(mclust)
# load data files
DATASET <- c("kumar2015","koh2016",  "trapnell2014","zhengmix2016", "simDataKumar", "simDataKumar2")
# define the datatype: "default", "filtered","unfiltered","optimalk"
datatype<-"unfiltered"

# reformat labels as factors, with integer levels

for (h in seq_len(length(DATASET)) ) {
  
  data_files <- as.list(paste0("results/run_results/cluster_single_",datatype,"_",DATASET[h],".rda"))
  
  load(data_files[[1]])

  
  clusters$labels <- lapply(clusters$labels, 
                            function(x) {
                              x <- as.factor(x)
                              levels(x) <- c(1:length(levels(x)))
                              x <- as.integer(x) %>%as.vector
                              return(x)
                              
  }
  )


  ### calculate ARI scores
  res.ari <- vector("list", length(clusters))
  names(res.ari) <- names(clusters)
  
  for (i in names(clusters)) {
    lab_i    <- clusters$labels[[i]]
    pred_i   <- unlist( clusters[[i]] )
    res.ari[[i]]$ARI <- tryCatch(adjustedRandIndex(pred_i,lab_i), error=function(e) NA )
    
  }
  save(res.ari, file = paste0("results/run_results/ari_single_",datatype,"_",DATASET[h],".rda"))
  
  print(res.ari)
  
}
### appendix

