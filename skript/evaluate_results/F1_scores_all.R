#################################
### Compute the F1 scores    ####
#################################
# SCript loads the clustering results from the single parameter run method and calculates the the F1 score or all datasets.

#load libraries
library(plyr)
library(dplyr)
library(clue)

#load helper file
source("skript/helper_files/helper_calc_f1_scores.R")
source("skript/helper_files/calc_f1_score.R")
source("skript/helper_files/Helper_functions.R")
# load data files, define which dataset to work with
DATASET <-   c("koh2016", "kumar2015", "trapnell2014", "zhengmix2016")

########
for (h in seq_len(length(DATASET)) ){

data_files <- as.list(paste0("results/run_results/",DATASET[h],".rda"))

load(data_files[[1]])
# reformat labels as factors
clusters$labels <- clusters$labels%>%as.factor

# create lookupfile
look <- lapply(clusters$labels, function(x) levels(x) <- c(1:length(levels(x))) )
levels(clusters$labels) <- c(1:length(levels(clusters$labels)))
clusters$labels <- as.integer(clusters$labels)


### calculate F1 scores
res.f1 <- vector("list", length(clusters))
names(res.f1) <- names(clusters)


for (i in 1:length(clusters)){
  
calc_f1_scores(clusters[[i]], clusters$labels)
  
}


for (j in seq_len(length(res.f1))) {
  act <- clusters$labels
  prd <- clusters[[j]]
  res.f1[[j]]  <-  calc_f1_score(labels=act,cluster = prd)  
}
print(res.f1)


save(res.f1, file = paste0("results/run_results/f1_",DATASET[h],".rda"))
}

calc_f1_score(labels[[i]], cluster[[i]]) 


#### Appendix
for (i in 1:7){
  return(length(clusters[[6]]))
  
} 
names(clusters[[7]])
