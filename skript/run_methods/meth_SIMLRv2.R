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

source("FILES.R")


# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
labels <- load_labels(data) 



# RUN SIMLR
list <- vector("list", length(data))
names(list) <- names(data)
res.SIMLR <- sys.time <- res.cluster <-  list
# Set paramaeters: c is the number of cluster expected. we define the tuning parameter k to 10 as it is the standart parameter, the number of dimension are NA
par.c <-  list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016 = 10,
  simDataKumar = 3
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
  lab_i <- as.data.frame(labels[[i]])
  write.table(lab_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "results/SIMLR/session_info_SIMLR.txt")
sessionInfo()
sink()

# Appendix


