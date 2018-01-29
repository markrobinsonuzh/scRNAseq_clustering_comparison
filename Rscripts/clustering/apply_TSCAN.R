## Apply TSCAN

suppressPackageStartupMessages({
  require(TSCAN)
  require(dplyr)
})

apply_TSCAN <- function(sce, params, k) {
  dat <- counts(sce)
  st <- system.time({
    tscan <- preprocess(dat, minexpr_percent = params$minexpr_percent, logbase = 2)
    cluster <- exprmclust(tscan, clusternum = k)$clusterid
  })
  list(st = st, cluster = cluster)
}
