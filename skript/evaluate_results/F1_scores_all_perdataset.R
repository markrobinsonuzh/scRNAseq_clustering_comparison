########################
# 
#######################
# is not workin
#load libraries
library(dplyr)
library(plyr)
library(clue)

#load helper file
source("skript/helper_files/calc_f1_score.R")
source("skript/helper_files/Helper_functions.R")


##########################
# for all Methods
##########################

METHOD <- c(SIMLR,Seurat)


for (j in length(METHOD)) {
DATA_DIR <- file.path("results",paste0(METHOD))


######################################################
### load and reformat files with the cell labels, "ground truth":
######################################################


files_labels <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_labels_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_labels_trapnell2014.txt")),
  zhengmix2016 = file.path(DATA_DIR, paste0(METHOD,"_labels_zhengmix2016.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_labels_koh2016.txt"))
)

# read in labels
labels <- read.labels(files_labels = files_labels)

# reformat labels as factors
labels <- sapply(labels,as.factor)
# create lookupfile # not necessary
look <- lapply(labels, function(x) levels(x) <- c(1:length(levels(x))) )
# 
for(i in names(labels)){
  levels(labels[[i]]) <- look[[i]]
  
}
labels <- sapply(labels, as.integer)


### load files with the cluster results:

files_clusters <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_clus_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_clus_trapnell2014.txt")),
  zhengmix2016 = file.path(DATA_DIR, paste0(METHOD,"_clus_zhengmix2016.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_clus_koh2016.txt"))
)

cluster <- read.cluster(files_clusters=files_clusters) %>% lapply(as.integer)

### Calculate f1 score for all data sets
res.f1 <- vector("list", length(cluster))
names(res.f1) <- names(cluster)

for (i in names(cluster)){ 
  {
  res.f1[[i]] <- calc_f1_score(labels[[i]], cluster[[i]]) 
    }
  

}


save(paste0("res.f1_",METHOD[j]) , file = paste0("results/run_results/f1_",METHOD[j],".rda"))


}
