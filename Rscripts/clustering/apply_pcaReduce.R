## Apply pcaReduce.

suppressPackageStartupMessages({
  library(pcaReduce)
})

apply_pcaReduce <- function(sce, params, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  dat <- t(logcounts(sce))
  L <- lapply(seq_len(n_rep), function(i) {  ## For each replication
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      st <- system.time({
        pca <- PCAreduce(dat, nbt = params$nbt, q = params$q, method = "S")[[1]]
        colnames(pca) <- paste0("k", (params$q + 1):2)
        cluster <- pca[, paste0("k", k)]
      })
      df <- data.frame(method = "pcaReduce", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(method = "pcaReduce",
                       run = i, 
                       k = k,
                       timing = st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"],
                       stringsAsFactors = FALSE, row.names = NULL)
      list(n_cluster = k, clusters = df, timing = tm)
    })  ## End for each k
    ## Summarize across different values of k
    assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
    timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
    list(assignments = assignments, timings = timings)
  })  ## End for each replication
  ## Summarize across different runs
  assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
  timings <- do.call(rbind, lapply(L, function(w) w$timings))
  
  list(assignments = assignments, timings = timings)
}
