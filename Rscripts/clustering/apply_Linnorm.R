## Apply Linnorm

suppressPackageStartupMessages({
  library(Linnorm)
})

apply_Linnorm <- function(sce, params, k) {
  tryCatch({
    dat <- assay(sce, "normcounts")
    st <- system.time({
      transformedExp <- Linnorm(dat, spikein = NULL, Filter = FALSE, 
                                minNonZeroPortion = params$minNonZeroPortion, 
                                BE_strength = 0.5)
      
      ## Cluster with predetermined number of clusters
      cluster <- Linnorm.tSNE(transformedExp, input = "Linnorm",
                              num_PC = params$num_PC, num_center = k)$k_means$cluster
      names(cluster) = colnames(dat)
    })
    ## Determine number of clusters automatically
    est_k <- length(unique(Linnorm.tSNE(transformedExp, input = "Linnorm",
                                        num_PC = params$num_PC, 
                                        num_center = params$range_clusters)$k_means$cluster))
    
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
