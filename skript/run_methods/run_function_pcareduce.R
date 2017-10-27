###################
# pcaReduce #######
###################

# pcaReduce uses a PCA and hierarchical clustering to find the number of clusters in the reduced dimension given by PCA. 
# The method expects that large classes of cells ar contained in low dimension PC representation
# and more refined (subsets) of these cells types are contained in higher dimensional PC representations.
# On the latent space (nxq) a kmeans clustering with q+1 clusters is performed. And the PC with the lowest variance is deleted. this process is iteratively repeated until only one single cluster remains.
# The resulting matrix has dimension nxq with q+1 clusters.
# Parameters to define are the number the method should be repeated, as pcaReduce is stochastic. Sampling without replacement; we choose nbt = 100, if the number samples is bigger than 100.
# The number of starting principal components q, we choose 50 as the default.
# And the number clusters n which is given by the "ground truth". the stepwise merging of the clusters can be done using sampling based merging (S) or
# merging based on largest probability (M).




#load libraries
source("skript/helper_files/Helper_functions.R")

library(dplyr)
library(pcaReduce)

# import data as sceset
# file paths

source("FILES.R")


# load data sets

data <- labels<- vector("list", length(files))

names(data) <-names(labels) <-  names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}


# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}

# create store files
list <- vector("list", length(data))
names(list) <- names(data)
input_matrix <- list
# extract expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]]) # use count scaled length scaled tpms, normalized and log2 transformed
}
# set parameters, nbt is number of times to repeat pcareduce; q is number of starting dimensions, n cluster the number of clusters
par.nbt <- list(
  kumar2015 = 50,
  trapnell2014 = 50,
  zhengmix2016 = 50,
  koh2016 = 50
  
)

par.q <- list(
  kumar2015 = 30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016 = 30
)
n.cluster <- list(
  kumar2015=3,
  trapnell2014=3,
  zhengmix2016=4,
  koh2016 = 10
)
# extract k dimension,  
# run pce Reduce, vary q
run_pcareduce <- function(input_matrix, par.nbt,par.q,n.cluster){
  
  require(dplyr)
  require(pcaReduce)
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <- sys.time<- pca.red <- list
# PCAREDUCE
for (i in names(input_matrix)){
    pca.red[[i]] <- PCAreduce(t(input_matrix[[i]]), nbt = par.nbt[[i]], q = par.q[[i]], method = 'S')[[1]]
    res.cluster[[i]]  <- pca.red[[i]][, (par.q[[i]]-9):par.q[[i]] ] # improve this!
    colnames(res.cluster[[i]]) <-  c(paste0( rev(seq_len(length( (par.q[[i]]-9):par.q[[i]] ) ))) ) # this as well
}
  return( res.cluster  )
}
# run PCAreduce
res.cluster <- run_pcareduce(input_matrix,par.nbt, par.q,n.cluster)


# save clusters

dir_cluster <- paste0("results/PCAreduce/PCAreduce_krange_clus_", names(res.cluster), ".txt")
save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/PCAreduce/PCAreduce_krange_systime_",names(sys.time),".txt")
save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("results/PCAreduce/PCAreduce_krange_labels_",names(labels), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
  
}

###### Save Session Info
sink(file = "results/PCAreduce/session_info_PCAreduce_krange.txt")
sessionInfo()
sink()

# Appendix

