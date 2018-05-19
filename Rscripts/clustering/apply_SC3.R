## Apply SC3

suppressPackageStartupMessages({
  library(scater)
  library(SC3)
})

apply_SC3 <- function(sce, params, k) {
  (seed <- round(1e6*runif(1)))
  tryCatch({
    rowData(sce)$feature_symbol <- rownames(counts(sce))
    st1 <- system.time({
      dat <- sc3_prepare(sce, gene_filter = params$gene_filter, 
                         pct_dropout_min = params$pct_dropout_min, 
                         pct_dropout_max = params$pct_dropout_max, 
                         svm_max = 1e6, n_cores = 1, rand_seed = seed)
    })
    est_k <- metadata(sc3_estimate_k(dat))$sc3$k_estimation
    st2 <- system.time({
      dat <- sc3(dat, ks = k, pct_dropout_min = params$pct_dropout_min,
                 pct_dropout_max = params$pct_dropout_max,
                 gene_filter = params$gene_filter, rand_seed = seed, n_cores = 1,
                 biology = FALSE, k_estimator = FALSE, svm_max = 1e6)
      cluster <- as.numeric(colData(dat)[, paste0("sc3_", k, "_clusters")])
      names(cluster) <- rownames(colData(dat))
    })
    
    st <- st1 + st2
    st <- st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"]
    list(st = st, cluster = cluster, est_k = est_k)
  },
  error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
