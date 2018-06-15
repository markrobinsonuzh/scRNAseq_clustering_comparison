## Apply pcaReduce.

suppressPackageStartupMessages({
  library(pcaReduce)
  library(clue)
})

apply_pcaReduce <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      pca <- PCAreduce(t(dat), nbt = params$nbt, q = params$q, method = "S")
      part <- lapply(pca, function(x) {
        colnames(x) <- paste0("k", (params$q + 1):2)
        as.cl_partition(x[, paste0("k", k)])
      })
      cons <- cl_consensus(as.cl_ensemble(part), method = "SE", control = list(nruns = 50))
      cluster <- c(cl_class_ids(cons))
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
