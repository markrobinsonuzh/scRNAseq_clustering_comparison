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

