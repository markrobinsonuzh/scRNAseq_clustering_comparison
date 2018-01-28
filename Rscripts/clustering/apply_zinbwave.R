suppressPackageStartupMessages({
  library(zinbwave)
  library(scater)
  library(Rtsne)
  library(dplyr)
})

filter_hvg <- function(data, n_genes) {
  filter <- rowSums((assay(data, "normcounts")) > 5) > 5 
  data <- data[filter, ]
  vars <- counts(data) %>% log1p %>% rowVars
  names(vars) <- rownames(data)
  vars <- sort(vars, decreasing = TRUE)
  data <- data[names(vars)[seq_len(n_genes)], ]
  return(data)
}

apply_zinbwave <- function(sce, parameters, n_rep) {
  ## Run repeatedly with a range of cluster numbers
  L <- lapply(seq_len(n_rep), function(i) {
    dat <- filter_hvg(sce, parameters$n_genes)
    tmp <- lapply(parameters$range_clusters, function(k) {
      st <- system.time({
        res <- zinbFit(round(assay(dat, "counts")), 
                       K = 2, epsilon = parameters$n_genes,
                       verbose = TRUE, nb.repeat.initialize = 2,
                       maxiter.optimize = 25, stop.epsilon.optimize = 1e-4)
        d <- dist(getW(res))
        cluster <- kmeans(d, centers = k)$cluster
      })
      df <- data.frame(method = "zinbwave", 
                       cell = names(cluster),
                       run = i,
                       k = k,
                       cluster = cluster,
                       stringsAsFactors = FALSE)
      tm <- data.frame(method = "zinbwave",
                       run = i, 
                       k = k,
                       timing = st["user.self"] + st["sys.self"] + st["user.child"] + st["sys.child"],
                       stringsAsFactors = FALSE)
      list(n_cluster = k, clusters = df, timing = tm)
    })
    assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
    timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
  })
  assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
  timings <- do.call(rbind, lapply(L, function(w) w$timings))
  
  list(assignments = assignments, timings = timings)
}
