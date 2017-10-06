#####################
# DBscan
#####################
# dbscan is a density based clustering method. The k neirest neighborhoud parameter and neighborhood size (epsilon) has to be defined.
# This is usually done using a k neirest neighbour distance plot. for this strategy the k-nearest neighbor distance is plotted (k-NN distance)
# this is the distance from a choosen point to its k nearest neighbors. As a rule of thumb k >= d+1, here d is the dimension of the data.
# A sharp change in the distance plot indicates the epsilon distance. 
# If only one single cluster is obtained, then often epsilon is too large or MinPts is too small.
# similary, if epsilon is too small or MinPts is too large then every point becomes a noise point.
# note that dbscan has problems with high dim data , so we should reduce the dimension, which is not done yet!
source("skript/helper_files/Helper_functions.R")



#load libraries

library(dbscan)

# import data as sceset
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
res.cluster <- res.dbscan <- sys.time<- input_matrix<- pca.red <- list


# extract transposed expression data

for (i in names(input_matrix)){
  input_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}

# perform linear dimension reduction using PCA
for (i in names(input_matrix)) {
  if ( nrow( input_matrix[[i]])>50 ) {
    input_matrix[[i]] <- prcomp(input_matrix[[i]], center=TRUE, scale = FALSE )$x[,1:50]
  } else {
    input_matrix[[i]] <-  prcomp(input_matrix[[i]], center=TRUE, scale = FALSE )$x[,1:nrow( input_matrix[[i]])/2]
  }
}


# RUN dbscan, k is nearest neighbor

par.k <- list(
  kumar2015 = 51,
  trapnell2014 = 51,
  xue2013 = 5,
  koh2016 = 51
  
)


# run k neirest neighbour distance plot to find epsilon
par(mfrow=c(2,2))
for (i in names(input_matrix)){

  kNNdistplot(input_matrix[[i]], k = par.k[[i]]) 
  
}

# run dbscan
par.eps <- list(
  kumar2015 = 150,
  trapnell2014 = 250,
  xue2013 = 250,
  koh2016 = 150
)
par.Pts <- list(
  kumar2015 = 51,
  trapnell2014 = 51,
  xue2013 = 5,
  koh2016 = 51
  
)

for(i in names(input_matrix)){
  
res.cluster[[i]] <- dbscan(input_matrix[[i]], eps = par.eps[[i]] ,minPts = par.Pts[[i]])$cluster

}



# save clusters

dir_cluster <- paste0("results/dbscan/dbscan_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/dbscan/dbscan_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/dbscan/dbscan_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = paste0("results/dbscan/session_info_dbscan.txt"))
sessionInfo()
sink()

### Appendix


