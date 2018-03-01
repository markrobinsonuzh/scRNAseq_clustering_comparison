args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
filterings <- strsplit(filterings, ",")[[1]]
names(filterings) <- filterings
methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(datasets)
print(filterings)
print(methods)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

## Read clustering results
res <- do.call(rbind, lapply(datasets, function(d) {
  do.call(rbind, lapply(filterings, function(f) {
    do.call(rbind, lapply(methods, function(m) {
      x <- readRDS(paste0("results/sce_", f, "_", d, "_", m, ".rds"))
      dplyr::full_join(x$assignments %>%
                         dplyr::select(dataset, method, cell, run, k, resolution, cluster, trueclass),
                       x$k_estimates %>%
                         dplyr::select(dataset, method, run, k, resolution, est_k)
      ) %>%
        dplyr::full_join(x$timings %>% dplyr::select(dataset, method, run, k, resolution, timing))
    }))
  }))
}))

saveRDS(res, file = outrds)

date()
sessionInfo()
