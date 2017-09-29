#####################
# SIMLR
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 
# 

source("skript/helper_files/Helper_functions.R")

#load libraries

library(SIMLR)
library(igraph)
library(scater)

set.seed(1234)

# file paths

DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_SRP073808.rda")
)

# load data sets

data <- labels<- vector("list", length(files))

names(labels) <- names(data) <- names(files)



for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}


# RUN SIMLR
list <- vector("list", length(data))
names(list) <- names(data)
res.SIMLR <- sys.time <- res.cluster <-  list
# Set paramaeters
par.c <-  list(
  kumar2015 = 3,
  trapnell2014 = 3,
  xue2013 = 8,
  koh2016 = 10
)

# 
for (i in names(data)){
  sys.time[[i]] <- system.time({
  res.SIMLR[[i]] = SIMLR(X = exprs(data[[i]]), c = par.c[[i]], no.dim = NA,k=10, if.impute=FALSE, normalize = FALSE,cores.ratio = 0) # use exprs slot of SCeset; log2, normalized count_lstpm
  })
  res.cluster[[i]] <- res.SIMLR[[i]]$y$cluster  

}

# save clusters

dir_cluster <- paste0("results/SIMLR/SIMLR_clus_", names(res.cluster), ".txt")

save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SIMLR/SIMLR_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/SIMLR/SIMLR_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "results/SIMLR/session_info_SIMLR.txt")
sessionInfo()
sink()

# Appendix


