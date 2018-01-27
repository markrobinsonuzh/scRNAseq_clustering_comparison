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

run_function_dbscan <- function(data, labels, par.k, par.eps, par.Pts , datatype) {
  require(dbscan)
# create store files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- res.dbscan <- sys.time<- input_matrix<- pca.red <- list


# extract transposed expression data

for (i in names(input_matrix)){
  input_matrix[[i]] <- t((assay(data[[i]], "normcounts"))) # use count scaled length scaled tpms, normalized and log2 transformed
}


# RUN dbscan, k is nearest neighbor


# run k neirest neighbour distance plot to find epsilon
par(mfrow=c(2,2))
for (i in names(input_matrix)){

  kNNdistplot(input_matrix[[i]], k = par.k[[i]]) 
  
}

# run dbscan


for(i in names(input_matrix)){
  print(i)
  sys.time [[i]] <- system.time({
  res.cluster[[i]] <- dbscan(input_matrix[[i]], eps = par.eps[[i]] ,minPts = par.Pts[[i]])$cluster  
  })

}



# save clusters

dir_cluster <- paste0("results/",datatype,"/dbscan/dbscan_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/dbscan/dbscan_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels <-  paste0("results/",datatype,"/dbscan/dbscan_labels_",names(labels), ".txt")
save_labels(labels, dir_labels )


###### Save Session Info
sink(file = paste0("results/",datatype,"/dbscan/session_info_dbscan.txt"))
sessionInfo()
sink()

}


