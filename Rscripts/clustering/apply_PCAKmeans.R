## Apply PCA + K-means. 

apply_PCAKmeans <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
      pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
      cluster <- kmeans(pca, centers = k, nstart = 25)$cluster
      names(cluster) = colnames(dat)
    })
    
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}