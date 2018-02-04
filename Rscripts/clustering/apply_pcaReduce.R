## Apply pcaReduce.

suppressPackageStartupMessages({
  library(pcaReduce)
})

apply_pcaReduce <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      pca <- PCAreduce(t(dat), nbt = params$nbt, q = params$q, method = "S")[[1]]
      colnames(pca) <- paste0("k", (params$q + 1):2)
      cluster <- pca[, paste0("k", k)]
    })
    
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
