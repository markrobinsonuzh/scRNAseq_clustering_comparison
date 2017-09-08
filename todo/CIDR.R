########################
# CIDR
########################
#####################
# tSNE + kmeans
#####################

#load libraries
source("skript/helper_files/Helper_functions.R")

library(Rtsne)
library(scater)
library(dplyr)
# file paths

DATA_DIR <- "data"
files <- list(
  
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_SRP073808.rda")
  
)

# load data sets

list<- vector("list", length(files))
names(list) <- names(files)

list->data->labels->tinput_matrix->sys.time->res.rtsne->res.cluster 

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}
# extract transposed expression data

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}
dim(tinput_matrix[[1]])

library(cidr)


sData <- scDataConstructor(t(tinput_matrix[[1]]))
sData <- determineDropoutCandidates(sData)
sData <- wThreshold(sData)
sData <- scDissim(sData)
sData <- scPCA(sData)
sData <- nPC(sData)

nCluster(sData)
sData <- scCluster(sData)


plot(sData@PC[,c(1,2)], col=cols,
     pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")


adjustedRandIndex(sData@clusters,labels[[1]])

