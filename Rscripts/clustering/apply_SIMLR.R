## Apply SIMLR

suppressPackageStartupMessages({
  library(SIMLR)
  library(igraph)
  library(scater)
})

apply_SIMLR <- function(sce, params, k) {
  dat <- logcounts(sce)
  st <- system.time({
    simlr = SIMLR(X = dat, c = k, no.dim = NA, k = params$k, 
                  if.impute = FALSE, normalize = params$normalize)
  })
  cluster <- simlr$y$cluster  
})
  list(st = st, cluster = cluster)
}
