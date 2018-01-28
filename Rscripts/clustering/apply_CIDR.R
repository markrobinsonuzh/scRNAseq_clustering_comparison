## Apply CIDR

suppressPackageStartupMessages({
  require(scater)
  require(dplyr)
  require(cidr)
})

apply_CIDR <- function(sce, params, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  dat <- assay(sce, "normcounts")
  L <- lapply(seq_len(n_rep), function(i) {  ## For each replication
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      st <- system.time({
        sData <- scDataConstructor(dat)
        sData <- determineDropoutCandidates(sData)
        sData <- wThreshold(sData, plotTornado = TRUE)
        sData <- scDissim(sData)
        sData <- scPCA(sData)
        sData <- nPC(sData)
          
        # nCluster(sData) # different methods todefine the number of clusters, optional
        sData <- scCluster(object = sData, nCluster = params$k, 
                           cMethod = "ward.D2", nPC = params$nPC)
        cluster <- sData@clusters
      })
      df <- data.frame(method = "CIDR", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(method = "CIDR",
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
