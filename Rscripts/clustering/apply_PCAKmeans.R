## Apply PCA + K-means. 

apply_PCAKmeans <- function(sce, params, k) {
  tryCatch({
    dat <- t(logcounts(sce))
    st <- system.time({
      pca <- prcomp(dat, center = TRUE, scale. = FALSE)$x[, seq_len(params$nPC), drop = FALSE]
      cluster <- structure(kmeans(pca, centers = k)$cluster,
                           names = rownames(dat))
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, nrow(dat)), names = rownames(dat)),
         est_k = NA)
  })
}