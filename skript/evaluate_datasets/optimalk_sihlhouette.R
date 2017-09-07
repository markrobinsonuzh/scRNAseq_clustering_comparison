#############################################
# sihlouette plot
############################################

#####################################
## Within sum of squares for clusters      
#####################################

### load libraries

source("skript/helper_files/Helper_functions.R")
library(cluster)
library(dplyr) 
library(magrittr)


#pdf("results/plots/optimalk_sihlhoutte.pdf")
### set seed

set.seed(1234)

# Directories

DATA_DIR <- "data"

files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
)


# load data sets

list <- vector("list", length(files))
names(list) <- names(files)

list -> data -> tinput_matrix -> res.si


for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  rm(res)
}

# extract transposed expression data

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}

# reduce dimension




##########################
# silhouette plot
#########################
# pam
par.k= list(
  kumar2015 = 3,
  trapnell2014 = 2,
  xue2013 = 2
)


for (i in names(tinput_matrix)){
  res.si[[i]]<- silhouette(pam(tinput_matrix[[i]], k=par.k[[i]] ))
  plot(res.si[[i]])
}
# kmeans
for (i in names(tinput_matrix)){
  res.si[[i]]<- silhouette(kmeans(tinput_matrix[[i]], k=par.k[[i]] ))
  plot(res.si[[i]])
}




