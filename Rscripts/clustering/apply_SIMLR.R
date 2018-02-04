## Apply SIMLR

suppressPackageStartupMessages({
  library(SIMLR)
  library(igraph)
})

apply_SIMLR <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      simlr = SIMLR(X = dat, c = k, no.dim = NA, k = params$k, 
                    if.impute = FALSE, normalize = FALSE, cores.ratio = 1)
      cluster <- simlr$y$cluster
      names(cluster) <- colnames(sce)
    })
    
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
