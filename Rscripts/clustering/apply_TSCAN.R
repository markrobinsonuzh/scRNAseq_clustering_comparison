## Apply TSCAN

suppressPackageStartupMessages({
  require(TSCAN)
})

apply_TSCAN <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      cluster <- exprmclust(dat, clusternum = k, modelNames = "VVV", reduce = TRUE)$clusterid
    })
    
    ## Determine number of clusters automatically
    est_k <- length(unique(exprmclust(dat, clusternum = params$range_clusters, 
                                      modelNames = "VVV", reduce = TRUE)$clusterid))
    
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
