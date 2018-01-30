## Apply t-SNE + K-means. 

suppressPackageStartupMessages({
  library(Rtsne)
})

apply_RtsneKmeans <- function(sce, params, k) {
  tryCatch({
    dat <- t(logcounts(sce))
    st <- system.time({
      rtsne <- Rtsne(X = dat, perplexity = params$perplexity, pca = TRUE, 
                     initial_dims = params$initial_dims, check_duplicates = FALSE)
      cluster <- structure(kmeans(rtsne$Y, centers = k)$cluster,
                           names = rownames(dat))
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, nrow(dat)), names = rownames(dat)),
         est_k = NA)
  })
}
