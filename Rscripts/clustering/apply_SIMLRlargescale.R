## Apply SIMLRlargescale

suppressPackageStartupMessages({
  library(SIMLR)
  library(igraph)
  library(scater)
})

apply_SIMLRlargescale <- function(sce, params, k) {
  dat <- logcounts(sce)
  st <- system.time({
    simlr = SIMLR_Large_Scale(X = dat, c = k, k = params$k, 
                  if.impute = FALSE, normalize = params$normalize)
  })
  cluster <- simlr$y$cluster  
})
list(st = st, cluster = cluster)
}

