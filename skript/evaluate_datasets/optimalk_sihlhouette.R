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
library(scater)


pdf("results/plots/optimalk_sihlhoutte.pdf")
#pdf("results/plots/optimalk_sihlhoutte.pdf")
### set seed

set.seed(1234)

# Directories

source("FILES.R")


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
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016= 10
  
)

for (i in names(tinput_matrix)){
  res.si[[i]]<- silhouette(pam(tinput_matrix[[i]], k=par.k[[i]] ))
  plot(res.si[[i]], main=names(data[i]))
}

# save the data 

save(res.si,file="results/number_k/ressilhouette.rda")

dev.off()




