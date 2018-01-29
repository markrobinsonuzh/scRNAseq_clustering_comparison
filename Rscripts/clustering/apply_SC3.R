## Apply SC3

suppressPackageStartupMessages({
  library(scater)
  library(DESeq2)
  library(SC3)
})

apply_SC3 <- function(sce, params, k) {
  dat <- counts(sce)
  st <- system.time({
    data <- sc3_prepare(dat, ks = params$ks)
    data <- sc3_estimate_k(data)
    ifelse(params$estimate_k ==TRUE, 
           params$k <- metadata(data)$sc3$k_estimation, 
           params$k <- k)
    data <- sc3(data, ks = k, pct_dropout_max = params$pct_dropout_max,
                gene_filter = TRUE, rand_seed = 1234, n_cores = 1)
    p_data <- colData(data)
    cluster <- p_data[, grep("sc3_", colnames(p_data))]
  })
  list(st = st, cluster = cluster)
}
