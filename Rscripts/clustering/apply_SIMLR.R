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
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
