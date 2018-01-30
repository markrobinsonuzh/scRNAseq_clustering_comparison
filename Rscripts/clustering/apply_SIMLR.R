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
                    if.impute = FALSE, normalize = params$normalize)
    })
    cluster <- simlr$y$cluster  
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
