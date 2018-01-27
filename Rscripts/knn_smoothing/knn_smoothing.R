#############################################
# # K-nearest neighbor smoothing for UMI-filtered scRNA-Seq data
# (R implementation)

# Author: Yun Yan <yun.yan@nyumc.org>
# Copyright (c) 2017 New York University
# adapted by A.Duo for dataset smoothing
#############################################

#source  K-nearest neighbor smoothing
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/knn-smoothing-master/knn_smooth.R')

# source helper files
source("skript/helper_files/Helper_functions.R")
# source directories
source("FILES.R")
# DATA_DIR
DATA_DIR <- "data"
# load data sets
data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# knn_smoothing of all datasets
# test
knn_smoothing(mat=counts(data[[1]]), k=5)
# function smotth the data sets
return.smooth <- function(data, k){
  mat.norm <- assay(data,"normcounts")
  mat.count <- assay(data,"counts")
  mat.norm <- knn_smoothing( mat.norm, k)
  mat.raw <-  knn_smoothing( mat.count, k)
  assay(data,"normcounts") <- mat.norm
  assay(data,"counts") <- mat.raw 
  return(data)
}
# return the smoothed data
smooth.data <- lapply(data, return.smooth, k =5)
dd <- smooth.data
# save the data
DATASET <-  list(
  kumar2015 ="GSE60749-GPL13112.rda",
  trapnell2014="GSE52529-GPL16791.rda",
  zhengmix2016="zhengmix.rda",
  koh2016="SRP073808.rda",
  simDataKumar="simDataKumar.rda"
)
for (i in 1:5){
res <- dd[[i]]
save(res, file=paste0( DATA_DIR, "/sceset_smooth_", DATASET[[i]] ) )
}

# appendix
save_data <- function(DATA,DATA_DIR,DATASET){
  save(DATA, file= paste0( DATA_DIR, "/sceset_smooth_", DATASET ))
}
lapply(names(res),save_data,DATA_DIR,DATASET )


