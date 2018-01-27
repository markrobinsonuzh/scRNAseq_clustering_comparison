###################
# pcaReduce #######
###################

# pcaReduce uses a PCA and hierarchical clustering to find the number of clusters in the reduced dimension given by PCA. 
# The method expects that large classes of cells ar contained in low dimension PC representation
# and more refined (subsets) of these cells types are contained in higher dimensional PC representations.
# On the latent space (nxq) a kmeans clustering with q+1 clusters is performed. And the PC with the lowest variance is deleted. this process is iteratively repeated until only one single cluster remains.
# The resulting matrix has dimension nxq with q+1 clusters.
# Parameters to define are the number of times the method should be repeated, as pcaReduce is stochastic. Sampling without replacement; we choose nbt = 100, if the number samples is bigger than 100.
# The number of starting principal components q, we choose 50 as the default.
# And the number clusters n which is given by the "ground truth". the stepwise merging of the clusters can be done using sampling based merging (S) or
# merging based on largest probability (M).




run_function_pcareduce <- function( data, labels, par.nbt, par.q , n.cluster , datatype ){
  require(pcaReduce )
# create store files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- input_matrix<- pca.red <- list
# extract expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]]) # use count scaled length scaled tpms, normalized and log2 transformed
}


# extract k dimension 

par.k <- function(i){
  (par.q[[i]]+2)-(n.cluster[[i]])
}

# run pce Reduce, vary q
for (i in names(input_matrix)){
  print(i)
  sys.time[[i]] <- system.time({
  pca.red[[i]] <- PCAreduce(t(input_matrix[[i]]), nbt = par.nbt[[i]], q = par.q[[i]], method = 'S')[[1]]
  res.cluster[[i]]  <- as.character(pca.red[[i]][ ,par.k(i)])
})
}

# save clusters

dir_cluster <- paste0("results/",datatype,"/PCAreduce/PCAreduce_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/PCAreduce/PCAreduce_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

dir_labels <-  paste0("results/",datatype,"/PCAreduce/PCAreduce_labels_",names(labels), ".txt")
save_labels(labels, dir_labels )


###### Save Session Info
sink(file = "results/",datatype,"/PCAreduce/session_info_PCAreduce.txt")
sessionInfo()
sink()

}

