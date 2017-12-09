######################################
# Interface functions for bootclust  #
######################################
#-----------------------------------------------------
# RtsneKmeans
rtsnekmeansCBI <- function(  data, k, perplexity ) {
  require(Rtsne)
  res.rtsne <- Rtsne(X= data , perplexity=perplexity , pca = TRUE)
  res.cluster <- as.integer( kmeans(res.rtsne$Y, centers = k )$cluster )
  #out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster==x } )
  
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "rtsnekmeansCBI"
  )
  
  return(list)
  
}

#-----------------------------------------------------
# CIDR
cidrCBI <- function ( data, par.k,par.nPC) {
  require(cidr)
  # RUN cidr
  sData <- vector("list", length(data))
  names(sData) <-  names(data)
  # define number of clusters.
      sData <- scDataConstructor(data) # creates a scData object with slots for the expression table, lib size, dropout candidates etc...
      sData <- determineDropoutCandidates(sData) # determines the dropout candidates
      sData <- wThreshold(sData, plotTornado = FALSE) # sets the imputation weighting threshold
      sData <- scDissim(sData) # computes the dissimilarity matrix for the slot dissim
      sData <- scPCA(sData) # performs PCA on the dissimilarity matrix
      sData <- nPC(sData) # deterimines the optimal number of PC to be used in clustering, populates nPC
      # nCluster(sData) # different methods todefine the number of clusters, optional
      sData <- scCluster(object = sData, nCluster = par.k, cMethod = "ward.D2", nPC = par.nPC)# hierarchical clustering on PC , nCluster if user defines cluster number, nPC is the number of PC used for clustering (default is 4), cMethod is hierarchical clustering method default is ward.D2
      res.cluster <- sData@clusters
      # out
      clusterlist <- vector("list", length( unique(res.cluster) ))
      clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
      list <- list(
        result=res.cluster,
        nc= length( unique(res.cluster) ),
        clusterlist=clusterlist,
        partition=res.cluster,
        clustermethod= "cidrCBI")
      
        return(list)
        
    
  }
#-----------------------------------------------------
# Linnorm
linnormCBI <- function ( data, par.minNonZeroPortion, par.num_center ,par.BE_strength) {
  require(Linnorm)
  transformedExp <- Linnorm(data, spikein=NULL, minNonZeroPortion = par.minNonZeroPortion, BE_strength=par.BE_strength)
  res.cluster <- Linnorm.tSNE(transformedExp, input = "Linnorm",num_center=par.num_center )$k_means$cluster
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "linnormCBI")
  
  return(list)
  
}
#-----------------------------------------------------
# SIMLR

simlrCBI <- function (data,  par.c, normalize) {
  require(SIMLR)
  require(igraph)
  require(scater)
  
  # RUN SIMLR
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.SIMLR <- res.cluster <-  list
  # run SIMLR
  res.SIMLR = SIMLR(X =data, c = par.c, no.dim = NA,k=10, if.impute=FALSE, normalize = normalize ) # use exprs slot of SCeset; log2, normalized count_lstpm
  res.cluster <- res.SIMLR$y$cluster  
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "simlrCBI")
  
  return(list)
}

#-----------------------------------------------------
# SC3

sc3CBI <- function ( data, par.ks, par.k_estimator ,par.k, pct_dropout_max ) {
  require("DESeq2")
  require("SC3")
  
  # list to store results
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <- list
  data<- sc3_prepare(data, ks=par.ks)# uses the exprs slot of SCEset ; log2transformed, non normalized data, filter data 
  data<- sc3_estimate_k(data) # optional estimate the number of clusters
  # estimated k (default ) or user supplied k (filtered, unfiltered)
  ifelse(par.k_estimator==TRUE, par.k <- metadata(data)$sc3$k_estimation, par.k <- par.k)
  # run sc3
  data<- sc3(data, ks = par.k , pct_dropout_max=pct_dropout_max,
             gene_filter = TRUE ,rand_seed=NULL, n_cores=2) # perform sc3 clustering
  
  # store clusters
  p_data <- colData(data)
  res.cluster <- p_data[ , grep("sc3_", colnames(p_data))]
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "sc3CBI")
  
  return(list)
  
}

