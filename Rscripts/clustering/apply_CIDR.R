## Apply CIDR

suppressPackageStartupMessages({
  require(cidr)
})

apply_CIDR <- function(sce, params, k) {
  tryCatch({
    dat <- assay(sce, "normcounts")
    st <- system.time({
      sData <- scDataConstructor(dat)
      sData <- determineDropoutCandidates(sData)
      sData <- wThreshold(sData, plotTornado = TRUE)
      sData <- scDissim(sData)
      sData <- scPCA(sData)
      sData <- nPC(sData)
      
      # nCluster(sData) # different methods todefine the number of clusters, optional
      sData <- scCluster(object = sData, nCluster = k, 
                         cMethod = "ward.D2", nPC = params$nPC)
      cluster <- sData@clusters
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
