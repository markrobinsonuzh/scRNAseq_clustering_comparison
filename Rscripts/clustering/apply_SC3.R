## Apply SC3

suppressPackageStartupMessages({
  library(scater)
  library(SC3)
})

apply_SC3 <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st1 <- system.time({
      data <- sc3_prepare(dat, ks = params$ks)
    })
    data <- sc3_estimate_k(data)
    st2 <- system.time({
      data <- sc3(data, ks = k, pct_dropout_max = params$pct_dropout_max,
                  gene_filter = TRUE, rand_seed = 1234, n_cores = 1)
      cluster <- colData(data)[, grep("sc3_", colnames(colData(data)))]
    })
    st <- st1 + st2
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = metadata(data)$sc3$k_estimation)
  },
  error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
