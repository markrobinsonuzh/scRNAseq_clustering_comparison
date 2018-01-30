## Apply ZINB-WaVE

suppressPackageStartupMessages({
  library(zinbwave)
  library(dplyr)
})

filter_hvg <- function(data, n_genes) {
  filter <- rowSums(counts(data)) > 5 
  data <- data[filter, ]
  vars <- counts(data) %>% log1p %>% rowVars
  names(vars) <- rownames(data)
  vars <- sort(vars, decreasing = TRUE)
  data <- data[names(vars)[seq_len(n_genes)], ]
  return(data)
}

apply_zinbwave <- function(sce, parameters, k) {
  tryCatch({
    dat <- round(counts(filter_hvg(sce, parameters$n_genes)))
    st <- system.time({
      res <- zinbFit(dat, K = 2, epsilon = parameters$n_genes,
                     verbose = TRUE, nb.repeat.initialize = 2,
                     maxiter.optimize = 25, stop.epsilon.optimize = 1e-4)
      d <- dist(getW(res))
      cluster <- kmeans(d, centers = k)$cluster
    })
    list(st = st, cluster = cluster)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
