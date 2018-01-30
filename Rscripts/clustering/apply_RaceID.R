## Apply RaceID

source("Rscripts/clustering/RaceID_class.R")

apply_RaceID <- function(sce, params, k) {
  tryCatch({
    dat <- as.data.frame(counts(sce))
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
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
