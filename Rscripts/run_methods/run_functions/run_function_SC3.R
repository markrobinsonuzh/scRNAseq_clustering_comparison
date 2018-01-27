####################
# SC3
######################
# change sc 3 prepare
# Input in SC3 is by default the exprs slot in a SCEset object. So the data should be normalized and log transformed. AS an option one can also set the argument gene_filter, if a filtering beforhand should br performed.
# SC3 uses distance measures of the filtered and log transformed expression matrix. PCA or Laplacian graphs are used for  dimension reduction. 
# Kmeans clustering is then performed on the d different dimensions. Using the d different clustering results a consensus matrix is computed. On the distances of this consenus matrix a hierarchical
# clustering step is performed. A range of number of clusters can be used by used. Its also possible to estimate the number of clusters.

# RUN SC3
run_function_sc3 <- function( data, labels, par.ks, par.k_estimator ,par.k, pct_dropout_max, datatype ){ 
  require("scater")
  require("DESeq2")
  require("SC3")

# list to store results
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time <- list

for (i in names(data)){
 
  print(i)
  sys.time[[i]] <- system.time({
    data[[i]]<- sc3_prepare(data[[i]], ks=par.ks[i])# uses the exprs slot of SCEset ; log2transformed, non normalized data, filter data 
    data[[i]]<- sc3_estimate_k(data[[i]]) # optional estimate the number of clusters
    # estimated k (default ) or user supplied k (filtered, unfiltered)
    ifelse(par.k_estimator[[i]]==TRUE, par.k[[i]] <- metadata(data[[i]])$sc3$k_estimation, par.k[[i]] <- par.k[[i]])
    # run sc3
    data[[i]]<- sc3(data[[i]], ks = par.k[[i]] , pct_dropout_max=pct_dropout_max[[i]],
                    gene_filter = TRUE ,rand_seed= 1234, n_cores=2) # perform sc3 clustering
  })
  # store clusters
  p_data <- colData(data[[i]])
  res.cluster[[i]] <- p_data[ , grep("sc3_", colnames(p_data))]
  
}
# save clusters

dir_cluster <- paste0("results/",datatype,"/SC3/sc3_clus_", names(data), ".txt")
save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/SC3/sc3_systime_",names(data),".txt")
save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels <-  paste0("results/",datatype,"/SC3/sc3_labels_",names(data), ".txt")
save_labels(labels, dir_labels )

###### Save Session Info
sink(file = "results/SC3/",datatype,"/session_info_sc3.txt")
sessionInfo()
sink()
}


