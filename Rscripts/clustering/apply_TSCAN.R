## Apply TSCAN

suppressPackageStartupMessages({
  library(TSCAN)
})

apply_TSCAN <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    ## Remove genes with variance = 0
    dat <- dat[rowVars(dat) > 0, ]
    st <- system.time({
      cluster <- exprmclust(dat, clusternum = k, modelNames = "VVV", reduce = TRUE)$clusterid
    })
    
    ## Determine number of clusters automatically
    est_k <- length(unique(exprmclust(dat, clusternum = params$range_clusters, 
                                      modelNames = "VVV", reduce = TRUE)$clusterid))
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
