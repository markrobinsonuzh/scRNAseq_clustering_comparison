## Apply t-SNE + K-means. 
## No automatic cluster number determination. 
## Possible to set the desired number of clusters
## Parameters: perplexity, initial_dims, range_clusters

suppressPackageStartupMessages({
  library(scater)
  library(Rtsne)
  library(dplyr)
})

apply_RtsneKmeans <- function(sce, params, k) {
  dat <- t(logcounts(sce))
  st <- system.time({
    rtsne <- Rtsne(X = dat, perplexity = params$perplexity, pca = TRUE, 
                   initial_dims = params$initial_dims, check_duplicates = FALSE)
    cluster <- structure(kmeans(rtsne$Y, centers = k)$cluster,
                         names = rownames(dat))
  })
  list(st = st, cluster = cluster)
}
