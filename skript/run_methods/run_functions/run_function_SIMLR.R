#####################
# SIMLR
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 
# 

run_function_simlrnormal <- function( data, labels, par.c, normalize, datatype  ) {
  require(SIMLR)
  require(igraph)
  require(scater)

# RUN SIMLR
list <- vector("list", length(data))
names(list) <- names(data)
res.SIMLR <- sys.time <- res.cluster <-  list
# 
for (i in names(data)){
  print(i)
  sys.time[[i]] <- system.time({
  res.SIMLR[[i]] = SIMLR(X = exprs(data[[i]]), c = par.c[[i]], no.dim = NA,k=10, if.impute=FALSE, normalize = normalize[[i]] ) # use exprs slot of SCeset; log2, normalized count_lstpm
  })
  res.cluster[[i]] <- res.SIMLR[[i]]$y$cluster  

}

# save clusters

dir_cluster <- paste0("results/",datatype,"/SIMLR/SIMLR_clus_", names(res.cluster), ".txt")

save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/SIMLR/SIMLR_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels<-  paste0("results/",datatype,"/SIMLR/SIMLR_labels_",names(labels), ".txt")
save_labels(labels, dir_labels )

###### Save Session Info
sink(file = "results/",datatype,"/SIMLR/session_info_SIMLR.txt")
sessionInfo()
sink()

}
# Appendix


