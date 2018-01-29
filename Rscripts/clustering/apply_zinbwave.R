## Apply ZINB-WaVE
## Input: Reduced count matrix (highly variable genes)
## No automatic cluster number determination. 
## Possible to set the desired number of clusters
## Parameters: n_genes

suppressPackageStartupMessages({
  library(zinbwave)
  library(scater)
  library(Rtsne)
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
  ## Run repeatedly with a range of cluster numbers
  dat <- filter_hvg(sce, parameters$n_genes)
  st <- system.time({
    res <- zinbFit(round(counts(dat)), 
                   K = 2, epsilon = parameters$n_genes,
                   verbose = TRUE, nb.repeat.initialize = 2,
                   maxiter.optimize = 25, stop.epsilon.optimize = 1e-4)
    d <- dist(getW(res))
    cluster <- kmeans(d, centers = k)$cluster
  })
  list(st = st, cluster = cluster)
}
