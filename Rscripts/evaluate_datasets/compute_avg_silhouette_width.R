args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datadir)
print(datasets)
print(ncores)
print(outrds)

## Compute silhouette widths for datasets
suppressPackageStartupMessages({
  library(cluster)
  library(scater)
  library(parallel)
})

datasets_full <- paste0(datadir, "/sce_full/sce_full_", datasets, ".rds")
names(datasets_full) <- gsub("\\.rds", "", basename(datasets_full))

full_data <- lapply(datasets_full, function(x) {
  readRDS(x)
})

## Euclidean distances from transposed count matrix
silh <- mclapply(full_data, function(x) {
  d <- dist(t(logcounts(x))) 
  s <- cluster::silhouette(
    as.integer(
      as.factor(
        colData(x)$phenoid)),
    d)
  summary(s)
}, mc.preschedule = FALSE, mc.cores = ncores)

saveRDS(silh, file = outrds)
date()
sessionInfo()
