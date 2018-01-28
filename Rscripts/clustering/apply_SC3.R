## Apply Linnorm

suppressPackageStartupMessages({
  library(scater)
  library(DESeq2)
  library(SC3)
})

apply_SC3 <- function(sce, params, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  dat <- counts(sce)
  L <- lapply(seq_len(n_rep), function(i) {  ## For each replication
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      st <- system.time({
        data <- sc3_prepare(dat, ks = params$ks)
        data <- sc3_estimate_k(data)
        ifelse(params$k_estimator ==TRUE, 
               params$k <- metadata(data)$sc3$k_estimation, 
               params$k <- params$k)
        data <- sc3(data, ks = params$k, pct_dropout_max = params$pct_dropout_max,
                    gene_filter = TRUE, rand_seed = 1234, n_cores = 1)
        p_data <- colData(data)
        cluster <- p_data[, grep("sc3_", colnames(p_data))]
      })
      df <- data.frame(method = "SC3", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(method = "SC3",
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
