## Apply t-SNE + K-means. 

suppressPackageStartupMessages({
  library(Rtsne)
})

apply_RtsneKmeans <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      rtsne <- Rtsne(X = t(dat), dims = params$dims, 
                     perplexity = params$perplexity, pca = TRUE, 
                     initial_dims = params$initial_dims, check_duplicates = FALSE)
      cluster <- kmeans(rtsne$Y, centers = k, nstart = 25)$cluster
      names(cluster) = colnames(dat)
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
