## Apply RaceID

suppressPackageStartupMessages({
  library(dplyr)
})
source("Rscripts/clustering/RaceID_class.R")

apply_RaceID <- function(sce, params, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  dat <- as.data.frame(counts(sce))
  L <- lapply(seq_len(n_rep), function(i) {  ## For each replication
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      st <- system.time({
        sc <- SCseq(dat)
        sc <- filterdata(sc, mintotal = params$mintotal, minexpr = params$minexprs, 
                         minnumber = params$minnumber, maxexpr = params$maxexpr, 
                         downsample = FALSE, dsn = 1, rseed = 1234)
        cluster <- clustexp(sc, metric = "pearson", cln = params$cln, 
                            do.gap = params$do.gap, clustnr = 20, B.gap = 50,
                            SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
                            bootnr = 50, rseed = 1234)@kmeans$kpart
      })
      df <- data.frame(method = "RaceID", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(method = "RaceID",
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
