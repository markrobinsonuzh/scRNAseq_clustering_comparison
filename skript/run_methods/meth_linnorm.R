
#####################
# Linnorm
#####################
# to do : include spikeinns
library("Linnorm")
library(scater)
require(dplyr)

source("skript/helper_files/Helper_functions.R")


# file paths

source("FILES.R")

# load data sets

data <- res.cluster <- vector("list", length(files))

names(data) <- names(res.cluster) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# load cell labels
labels <- load_labels(data) 


# RUN linnorm

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# Set a lower cutoff for the zhengmix data
par.minNonZeroPortion <- list(
  kumar2015 = 0.75,
  trapnell2014 =  0.75,
  zhengmix2016 =  0.25,
  koh2016 =  0.75,
  simDataKumar=0.75
)
# create store vectors
list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->transformedExp->res.cluster 

for (i in names(data)){
  sys.time [[i]] <- system.time({
    transformedExp[[i]] <- Linnorm(counts(data[[i]]), spikein=NULL, minNonZeroPortion = par.minNonZeroPortion[[i]])
    res.cluster[[i]] <- Linnorm.tSNE(transformedExp[[i]], input = "Linnorm")$k_means$cluster
  })
  
  
}
# save clusters

dir_cluster <- paste0("results/linnorm/linnorm_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/linnorm/linnorm_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/linnorm/linnorm_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/tscan/session_info_linnorm.txt")
sessionInfo()
sink()

### Appendix
plotmclust()

