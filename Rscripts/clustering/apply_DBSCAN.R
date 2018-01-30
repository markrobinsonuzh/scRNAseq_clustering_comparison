## Apply DBSCAN

suppressPackageStartupMessages({
  library(dbscan)
})

apply_DBSCAN <- function(sce, params, k) {
  tryCatch({
    dat <- t(logcounts(sce))
    st <- system.time({
      cluster <- dbscan(dat, eps = params$eps ,minPts = params$Pts)$cluster  
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, nrow(dat)), names = rownames(dat)),
         est_k = NA)
  })
}
