## Apply RaceID

source("Rscripts/clustering/RaceID_class.R")

apply_RaceID <- function(sce, params, k) {
  (seed <- round(1e6*runif(1)))
  tryCatch({
    dat <- as.data.frame(counts(sce))
    st <- system.time({
      sc <- SCseq(dat)
      sc <- filterdata(sc, mintotal = params$mintotal, minexpr = params$minexpr, 
                       minnumber = params$minnumber, maxexpr = params$maxexpr, 
                       downsample = FALSE, dsn = 1, rseed = seed)

      ## Cluster with predetermined number of clusters
      cluster <- clustexp(sc, metric = "pearson", cln = k, 
                          do.gap = FALSE, clustnr = max(params$range_clusters), B.gap = 50,
                          SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
                          bootnr = 50, rseed = seed)@kmeans$kpart
    })
    ## Determine number of clusters automatically
    est_k <- length(unique(clustexp(sc, metric = "pearson", cln = 0, 
                                    do.gap = TRUE, clustnr = max(params$range_clusters), B.gap = 50,
                                    SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
                                    bootnr = 50, rseed = seed)@kmeans$kpart))
    
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
