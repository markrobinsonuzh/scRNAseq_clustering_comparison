## Apply DBSCAN

suppressPackageStartupMessages({
  library(dbscan)
})

apply_DBSCAN <- function(sce, params, k) {
  dat <- t(logcounts(sce))
  st <- system.time({
    cluster <- dbscan(dat, eps = params$eps ,minPts = params$Pts)$cluster  
  })
  list(st = st, cluster = cluster)
}
