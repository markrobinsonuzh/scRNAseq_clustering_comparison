#####################
# SIMLR
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. Other parameters are if the input data should be transposed ore normalized, the number of working dimension in kmeans steo
# 

source("skript/helper_files/Helper_functions.R")

#load libraries
library(dplyr)
library(SIMLR)
library(igraph)
library(scater)

set.seed(1234)

# file paths

source("FILES.R")


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
# Set parameter number of cluster c
par.k <-  list(
  kumar2015 = c(2:10),
  trapnell2014 = c(2:10),
  zhengmix2016 =  c(2:10),
  koh2016 =  c(2:15),
  simDataKumar =  c(2:10)
)
### SIMLR function
run_simlr <- function(data, par.k){
  require(SIMLR)
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <-  list
 # 

for (i in names(data)){
  
  df.clus <- matrix( nrow = ncol(data[[i]]), ncol = length(par.k[[i]]) ) 
  
  for  ( j in seq_len( length( par.k[[i]]) ) ){
  df.clus[,j] <- SIMLR(X = assay(data[[i]], "normcounts"), c = par.k[[i]][j], no.dim = NA,k=10, if.impute=FALSE, normalize = FALSE, cores.ratio = 1)$y$cluster # use exprs slot of SCeset; log2, normalized count_lstpm
  }
  colnames(df.clus) <-  c( paste0(par.k[[i]]) )
  res.cluster[[i]] <- df.clus
}
return(res.cluster)  
}



# run simlr
res.cluster <- run_simlr(data, par.k)
# save clusters

dir_cluster <- paste0("results/SIMLR/SIMLR_krange_clus_", names(res.cluster), ".txt")

save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SIMLR/SIMLR_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/SIMLR/SIMLR_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "results/SIMLR/session_info_SIMLR_krange.txt")
sessionInfo()
sink()

# Appendix


