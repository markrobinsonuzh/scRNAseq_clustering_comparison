## Apply TSCAN

suppressPackageStartupMessages({
  require(TSCAN)
})

apply_TSCAN <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      tscan <- preprocess(dat, minexpr_percent = params$minexpr_percent, logbase = 2)
      cluster <- exprmclust(tscan, clusternum = k)$clusterid
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
