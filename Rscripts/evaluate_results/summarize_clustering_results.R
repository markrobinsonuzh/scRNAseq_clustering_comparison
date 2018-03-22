args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasetssmall <- strsplit(datasetssmall, ",")[[1]]
names(datasetssmall) <- datasetssmall
datasetsbig <- strsplit(datasetsbig, ",")[[1]]
names(datasetsbig) <- datasetsbig
filterings <- strsplit(filterings, ",")[[1]]
names(filterings) <- filterings
methodssmall <- strsplit(methodssmall, ",")[[1]]
names(methodssmall) <- methodssmall
methodsbig <- strsplit(methodsbig, ",")[[1]]
names(methodsbig) <- methodsbig

methods <- union(methodssmall, methodsbig)
datasets <- union(datasetssmall, datasetsbig)

print(datasetssmall)
print(datasetsbig)
print(filterings)
print(methodssmall)
print(methodsbig)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

## Read clustering results
res <- do.call(rbind, lapply(datasets, function(d) {
  do.call(rbind, lapply(filterings, function(f) {
    do.call(rbind, lapply(methods, function(m) {
      if ((m %in% methodssmall && d %in% datasetssmall) || 
          (m %in% methodsbig && d %in% datasetsbig)) {
        x <- readRDS(paste0("results/sce_", f, "_", d, "_", m, ".rds"))
        dplyr::full_join(x$assignments %>%
                           dplyr::select(dataset, method, cell, run, k, resolution, cluster, trueclass),
                         x$k_estimates %>%
                           dplyr::select(dataset, method, run, k, resolution, est_k)
        ) %>%
          dplyr::full_join(x$timings %>% dplyr::select(dataset, method, run, k, resolution, timing))
      }
    }))
  }))
}))

saveRDS(res, file = outrds)

date()
sessionInfo()
