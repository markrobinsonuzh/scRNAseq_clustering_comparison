## Apply monocle

suppressPackageStartupMessages({
  library(monocle)
  library(scran)
})

apply_monocle <- function(sce, params, k) {
  tryCatch({
    cds <- convertTo(sce, type = "monocle")
    cds <- tryCatch({
      estimateDispersions(cds)
    }, error = function(e) {
      cds
    })
    
    st <- system.time({
      cds <- reduceDimension(cds, max_components = params$max_components, 
                             num_dim = params$num_dim,
                             reduction_method = "tSNE", verbose = TRUE)
      cds <- clusterCells(cds, num_clusters = k + 1, method = "densityPeak")
      cluster <- cds$Cluster
      names(cluster) <- colnames(cds)
    })
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  },
  error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA),
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
