## Apply CIDR
## Input: normalized count matrix
## No automatic cluster number determination ???
## Possible to set the desired number of clusters
## Parameters: nPC

suppressPackageStartupMessages({
  require(scater)
  require(dplyr)
  require(cidr)
})

apply_CIDR <- function(sce, params, k) {
  ## Run repeatedly with a range of cluster numbers
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
  list(st = st, cluster = cluster)
}
