## Apply Linnorm

suppressPackageStartupMessages({
  library(Linnorm)
})

apply_Linnorm <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      transformedExp <- Linnorm(dat, spikein = NULL, 
                                minNonZeroPortion = params$minNonZeroPortion, 
                                BE_strength = params$BE_strength)
      cluster <- Linnorm.tSNE(transformedExp, input = "Linnorm",
                              num_center = k)$k_means$cluster
      cluster <- structure(cluster, names = colnames(dat))
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
