########################
# CIDR
########################
# Notes ; could assign nPCs in function directly

run_function_cidr <- function( data, labels, par.k, par.nPC,datatype  ) {
  
  require(scater)
  require(dplyr)
  require(cidr)
# create containers
  list<- vector("list", length(data))
  names(list) <- names(data)
  
  list->tinput_matrix->sys.time->res.cluster 
  
# extract transposed expression data

for (i in names(data) ){
  tinput_matrix[[i]] <- t(exprs( data[[i]])) # use  length scaled tpms, normalized and log2+1 transformed
}

# RUN cidr
sData <- vector("list", length(data))
names(sData) <-  names(data)
# define number of clusters.


for  (i in names(sData)) {
print(i)
sys.time[[i]] <- system.time({
sData[[i]] <- scDataConstructor(t(tinput_matrix[[i]])) # creates a scData object with slots for the expression table, lib size, dropout candidates etc...
sData[[i]] <- determineDropoutCandidates(sData[[i]]) # determines the dropout candidates
sData[[i]] <- wThreshold(sData[[i]], plotTornado = TRUE) # sets the imputation weighting threshold
sData[[i]] <- scDissim(sData[[i]]) # computes the dissimilarity matrix for the slot dissim
sData[[i]] <- scPCA(sData[[i]]) # performs PCA on the dissimilarity matrix
sData[[i]] <- nPC(sData[[i]]) # deterimines the optimal number of PC to be used in clustering, populates nPC

# nCluster(sData) # different methods todefine the number of clusters, optional
sData[[i]] <- scCluster( object=sData[[i]], nCluster = par.k[[i]], cMethod = "ward.D2", nPC = par.nPC[[i]] )# hierarchical clustering on PC , nCluster if user defines cluster number, nPC is the number of PC used for clustering (default is 4), cMethod is hierarchical clustering method default is ward.D2
# hierarchical clustering on PC , nCluster if user defines cluster number, nPC is the number of PC used for clustering (default is 4), cMethod is hierarchical clustering method default is ward.D2

})

res.cluster[[i]] <- sData[[i]]@clusters

}


# save clusters

dir_cluster <- paste0( "results/",datatype,"/cidr/cidr_clus_", names(res.cluster), ".txt" )


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0( "results/",datatype,"/cidr/cidr_systime_",names(sys.time),".txt" )

save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels <-  paste0( "results/",datatype,"/cidr/cidr_labels_", names(labels), ".txt" )
save_labels(labels, dir_labels )

###### Save Session Info

sink( file = paste0("results/",datatype,"/cidr/session_info_cidr.txt" ) )
sessionInfo()
sink()


}

# appendx

