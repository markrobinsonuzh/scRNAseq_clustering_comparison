###################
# pcaReduce
###################

source("~/Desktop/masterthesis/skript/helper_functions/Helper_functions.R")


#load libraries

library(pcaReduce)

# import data as sceset
# file paths

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
)

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
res.cluster <- sys.time<- input_matrix<- pca.red <- list





# extract expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]]) # use count scaled length scaled tpms, normalized and log2 transformed
}

# set parameters, nbt is number of times to repeat pcareduce; q is number of starting dimensions, n cluster the number of clusters
par.nbt <- list(
  kumar2015 <- 10,
  trapnell2014 <- 10,
  xue2013 <- 10
  
)

par.q <- list(
  kumar2015 <- 30,
  trapnell2014 <- 30,
  xue2013 <- 8
  
)
n.cluster <- list(
  kumar2015=3,
  trappnell=3,
  xue2013=8
)

# extract k dimension 

par.k <- function(i){
  (par.q[[i]]+2)-(n.cluster[[i]])
}


# run pce Reduce, vary q
for (i in names(input_matrix)){
  sys.time[[i]] <- system.time({
  pca.red[[i]] <- PCAreduce(t(input_matrix[[i]]), nbt = par.nbt[[i]], q = par.q[[i]], method = 'S')[[1]]
  pData(data[[i]])$pcaReduce <- as.character(pca.red[[i]][ ,par.k(i)])
})
  res.cluster[[i]] <- pData(data[[i]])$pcaReduce 
}

# save clusters

dir_cluster <- paste0("~/Desktop/masterthesis/results/PCAreduce/PCAreduce_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("~/Desktop/masterthesis/results/PCAreduce/PCAreduce_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/PCAreduce/PCAreduce_labels_",names(labels), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
  
}

###### Save Session Info
sink(file = "~/Desktop/masterthesis/results/PCAreduce/session_info_PCAreduce.txt")
sessionInfo()
sink()

# Appendix

