#####################
# DBscan
#####################
# dbscan is a density based clustering method. The k neirest neighborhoud parameter and neighborhood size (epsilon) has to be defined.
# This is usually done using a k neirest neighbour distance plot. for this strategy the k-nearest neighbor distance is plotted (k-NN distance)
# this is the distance from a choosen point to its k nearest neighbors. As a rule of thumb k >= d+1, here d is the dimension of the data.
# A sharp change in the distance plot indicates the epsilon distance. 
# If only one single cluster is obtained, then often epsilon is too large or MinPts is too small.
# similary, if epsilon is too small or MinPts is too large then every point becomes a noise point.
# note that dbscan has problems with high dim data , so we should reduce the dimension, which is done using PCA and working on latent space with 50 dim.
source("skript/helper_files/Helper_functions.R")


#load libraries

library(dbscan)

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
res.cluster <- res.dbscan <- sys.time<- input_matrix<- pca.red <- list


# extract transposed expression data

for (i in names(input_matrix)){
  input_matrix[[i]] <- t(assy(data[[i]]), "normcounts") # use count scaled length scaled tpms, normalized and log2 transformed
}



# parameter k is nearest neighbor, as a rule of thumb k >= dimension+1, here we use 10 % of the cell data set
par.k <- list(
  kumar2015 = nrow(input_matrix[[1]])*0.1,
  trapnell2014 = nrow(input_matrix[[2]])*0.1,
  zhengmix2016 = nrow(input_matrix[[3]])*0.1,
  koh2016 = nrow(input_matrix[[4]])*0.1
)

# run k neirest neighbour distance plot to find epsilon
par(mfrow=c(2,2))
for (i in names(input_matrix)){
  
  kNNdistplot(input_matrix[[i]], k = par.k[[i]]) 
  
}
# RUN dbscan, k is nearest neighbor

# parameter epsilon, size of the epsilon neighborhood. Bent in kNN dist plot. choose range for.
par.eps <- list(
  kumar2015 = c(150,200,230,240,250,280,300),
  trapnell2014 = c(380,400,420,440,460),
  zhengmix2016 = c(30,32,34,36,38,40),
  koh2016 = c(380,400,420,440,460)
)
names(par.eps) <- names(files)
# parameter Pts, number of minimum points in the eps region (for core points). Default is 5 points.

par.Pts <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 20,
  koh2016 = 5
  
)

# function dbscan
run_dbscan <-  function( par.eps, par.Pts ){
  
  for(i in names(input_matrix)){
    df.clus <- matrix( nrow = nrow(input_matrix[[i]]), ncol = length(par.eps[[i]]) )
    for (j in seq_len( length( par.eps[[i]])) ){
      df.clus[,j] <- dbscan(input_matrix[[i]], eps = par.eps[[i]][j] ,minPts = par.Pts[[i]])$cluster
    }
    colnames(df.clus) <-  c( paste0(par.eps[[i]]) )
    res.cluster[[i]] <- df.clus
  }
  return(res.cluster)
}
# run dbscan

res.cluster <-  run_dbscan(par.eps, par.Pts)
#
# save clusters

dir_cluster <- paste0("results/dbscan/dbscan_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/dbscan/dbscan_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/dbscan/dbscan_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = paste0("results/dbscan/session_info_dbscan_krange.txt"))
sessionInfo()
sink()

### Appendix


