## Apply ascend

suppressPackageStartupMessages({
  library(ascend)
  library(BiocParallel)
  library(class)
})

apply_ascend <- function(sce, params, k) {
  register(MulticoreParam(workers = 1, progressbar = TRUE), default = TRUE)
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      emset <- NewEMSet(ExpressionMatrix = dat)
      normset <- NormaliseByRLE(emset)
      pcaset <- RunPCA(normset)
      pcaset <- ReduceDimensions(pcaset, n = params$nPC)
      ## Try to cluster all cells. If it fails (if there are outliers), cluster
      ## the ones that can be clustered and assign the remaining ones to a
      ## cluster with kNN on the extracted PCs.
      resset <- RunCORE(pcaset, conservative = TRUE, nres = 40, 
                        remove_outlier = TRUE)
      ## Select the height to use
      ks <- resset@Clusters$KeyStats
      ks <- subset(ks, ClusterCount == k)
      height <- ks$Height[which.max(ks$RandIndex)]
      cluster <- resset@Clusters$ClusteringMatrix[, height]
      names(cluster) <- rownames(resset@Clusters$ClusteringMatrix)
      ## Classify the left-out cells
      train <- names(cluster)
      test <- setdiff(rownames(pcaset@PCA$PCA), names(cluster))
      cls <- as.numeric(as.character(class::knn(train = pcaset@PCA$PCA[train, , drop = FALSE], 
                                                test = pcaset@PCA$PCA[test, , drop = FALSE], 
                                                cl = cluster, 
                                                k = 5)))
      names(cls) <- test
      cluster <- c(cluster, cls)
      cluster <- cluster[match(colnames(dat), names(cluster))]
      est_k = resset@Clusters$NumberOfClusters
    })

    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)), 
         est_k = NA)
  })
}
