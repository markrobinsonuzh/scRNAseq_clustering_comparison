#################################
### Compute the F1 scores    ####
#################################
# SCript loads the clustering results from the single parameter run method and calculates the the F1 score for all datasets.
# F1 scores are stored as f1_dataset.rda  file in run_results directory
#load libraries
library(plyr)
library(dplyr)
library(clue)

#load helper file
source("skript/helper_files/calc_f1_score.R")
source("skript/helper_files/Helper_functions.R")
# load data files, define which dataset to work with
DATASET <-   c("koh2016", "kumar2015", "trapnell2014", "zhengmix2016", "simDataKumar","simDataKumar2")
# which datatype
datatype <-  "unfiltered"

# 
########
for (h in seq_len(length(DATASET)) ){
  
  data_files <- as.list(paste0("results/run_results/cluster_single_",datatype,"_",DATASET[h],".rda"))
  
  load(data_files[[1]])
  # reformat labels as factors
  clusters$labels <- lapply(clusters$labels, 
                            function(x) {
                              x <- as.factor(x)
                              levels(x) <- c(1:length(levels(x)))
                              x <- as.integer(x) %>%as.vector
                              
                            }
  )  

  
  ### calculate F1 scores
  res.f1 <- vector("list", length(clusters))
  names(res.f1) <- names(clusters)
  
  
  for (j in names(res.f1) ) {
    act <- clusters$labels[[j]]
    prd <- clusters[[j]]
    res.f1[[j]]  <-  tryCatch(calc_f1_score(labels=act,cluster = prd)  , error=function(e) NA)
  }
  print(res.f1)
  
  save(res.f1, file = paste0("results/run_results/f1_single_",datatype,"_",DATASET[h],".rda"))
}

#### Appendix
