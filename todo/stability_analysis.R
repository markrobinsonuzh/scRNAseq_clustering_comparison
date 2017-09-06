###########################
# Stability analysis
###########################

### load libraries
source("~/Desktop/masterthesis/skript/Helper_functions.R")
library(cluster)
library(dplyr)
library(fpc)
library(scater)
library(clValid)
library(kohonen)
library(mclust)


# load data

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
  
)

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
  
)

# load data sets

data <- vector("list", length(files))
input_matrix<- data
names(data) <- names(input_matrix) <-  names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# extract transposed expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]]) # use count scaled length scaled tpms, normalized and log2 transformed
}
input_matrix[[i]] <- t(input_matrix[[i]])

##############################################££££££££££££££££££££££
#  Perform multivariate validation of cluster results with clValid
#####################################################################

valid_test <- clValid(input_matrix[[3]] , c(2:4, 15),
                      clMethods = c("hierarchical", "kmeans",  "pam" ),
                      validation = c("internal", "stability")
)


# save the results 


##############################################££££££££££££££££££££££
#  Clusterwise cluster stability assessment by resampling 
#####################################################################


km.boot2 <- clusterboot(input_matrix[[3]], B=100, bootmethod="boot",
                        clustermethod=kmeansCBI,
                        krange=2, seed=20)
km.boot3 <- clusterboot(input_matrix[[3]], B=100, bootmethod="boot",
                        clustermethod=kmeansCBI,
                        krange=3, seed=20)
km.boot4 <- clusterboot(input_matrix[[3]], B=100, bootmethod="boot",
                        clustermethod=kmeansCBI,
                        krange=4, seed=20)
# save the results
