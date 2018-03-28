## Apply SIMLRlargescale

suppressPackageStartupMessages({
  library(SIMLR)
  library(igraph)
})

apply_SIMLRlargescale <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      simlr = SIMLR_Large_Scale(X = dat, c = k, k = params$k, 
                                if.impute = FALSE, normalize = FALSE)
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

