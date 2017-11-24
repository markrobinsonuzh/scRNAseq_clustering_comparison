######################################
# Interface functions for bootclust  #
######################################
#-----------------------------------------------------
# RtsneKmeans
rtsnekmeansCBI <- function(  data,k,perplexity=par.perp ) {
  require(Rtsne)
  res.rtsne <- Rtsne(X= data ,perplexity=par.perp , pca = TRUE)
  res.cluster <- as.integer( kmeans(res.rtsne$Y, centers = k )$cluster )
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){   res.cluster==x    } )
  
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster
  )
  
  return(list)
  
}

#-----------------------------------------------------