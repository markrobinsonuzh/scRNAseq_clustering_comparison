## Apply pcaReduce.

suppressPackageStartupMessages({
  library(pcaReduce)
})

apply_pcaReduce <- function(sce, params, k) {
  dat <- t(logcounts(sce))
  st <- system.time({
    pca <- PCAreduce(dat, nbt = params$nbt, q = params$q, method = "S")[[1]]
    colnames(pca) <- paste0("k", (params$q + 1):2)
    cluster <- pca[, paste0("k", k)]
  })
  list(st = st, cluster = cluster)
}
