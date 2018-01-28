## Apply t-SNE + K-means. 
## No automatic cluster number determination. 
## Possible to set the desired number of clusters
## Parameters: perplexity, initial_dims, k (number of clusters)

suppressPackageStartupMessages({
  library(scater)
  library(Rtsne)
  library(dplyr)
})

apply_RtsneKmeans <- function(sce, params, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  L <- lapply(seq_len(n_rep), function(i) {
    dat <- t(logcounts(sce))
    tmp <- lapply(params$range_clusters, function(k) {
      st <- system.time({
        rtsne <- Rtsne(X = dat, perplexity = params$perplexity, pca = TRUE, 
                       initial_dims = params$initial_dims, check_duplicates = FALSE)
        cluster <- structure(kmeans(rtsne$Y, centers = k)$cluster,
                             names = rownames(dat))
      })
      df <- data.frame(method = "RtsneKmeans", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(method = "RtsneKmeans",
                       run = i, 
                       k = k,
                       timing = st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"],
                       stringsAsFactors = FALSE, row.names = NULL)
      list(n_cluster = k, clusters = df, timing = tm)
    })
    assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
    timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
    list(assignments = assignments, timings = timings)
  })
  assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
  timings <- do.call(rbind, lapply(L, function(w) w$timings))
  
  list(assignments = assignments, timings = timings)
}
