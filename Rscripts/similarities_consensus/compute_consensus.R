## -------------------------------------------------------------------------- ##
## Compute consensus using the clue package
## -------------------------------------------------------------------------- ##

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(ncores)
print(outrds)

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(clue)
  library(multidplyr)
  library(ggplot2)
  library(parallel)
})

res <- readRDS(file = clusteringsummary)

## Help function to calculate consensus with the clue package. 
## Input: data frame with columns cell, cluster, run
cons.clue <- function(df) {
  m <- tryCatch({
    df %>% dplyr::select(cell, run, cluster) %>%
      tidyr::spread(key = run, value = cluster) %>%
      as.data.frame() %>% 
      tibble::column_to_rownames("cell") %>% as.matrix()
    }, error = function(e) NA)
  if (all(is.na(m))) {
    re <- data.frame(cell = unique(df$cell), consensus = NA, 
                     stringsAsFactors = FALSE)
  } else {
    re <- plyr::alply(m, 2, clue::as.cl_partition)
    ## Replace NAs by zeros
    re <- lapply(re, function(x) { 
      x <- as.matrix(x$.Data)
      x[is.na(x)] <- 0
      x
    }) 
    re <- clue::as.cl_ensemble(re)
    re <- clue::cl_consensus(re, method = "SE", control = list(nruns = 50))
    re <- clue::cl_class_ids(re) 
    re <- data.frame(cell = names(re), consensus = as.integer(re),
                     stringsAsFactors = FALSE)
  }
  return(re) 
}

## compute consensus
ds <- unique(res$dataset)
mt <- setdiff(unique(res$method), "Seurat")
ks <- unique(res$k)
rs <- setdiff(unique(res$resolution), NA)

res_consclue <- do.call(rbind, mclapply(ds, function(d) {
  do.call(rbind, lapply(mt, function(m) {
    do.call(rbind, lapply(ks, function(u) {
      tmp <- res %>% dplyr::filter(dataset == d & method == m & k == u)
      if (nrow(tmp) > 0)
        tmp %>% dplyr::select(dataset, method, cell, k, resolution, 
                              trueclass, est_k, elapsed) %>% 
        dplyr::left_join(cons.clue(tmp), by = "cell") %>% 
        dplyr::mutate(elapsed = median(elapsed),
                      est_k = median(est_k)) %>%
        dplyr::distinct()  ## only need to keep one of the five rows corresponding to the different runs
    }))
  }))
}, mc.preschedule = FALSE, mc.cores = ncores))

res_consclue_seurat <- do.call(rbind, mclapply(ds, function(d) {
  do.call(rbind, lapply(rs, function(u) {
    tmp <- res %>% dplyr::filter(dataset == d & method == "Seurat" & resolution == u)
    if (nrow(tmp) > 0)
      tmp %>% dplyr::select(dataset, method, cell, k, resolution, 
                            trueclass, est_k, elapsed) %>% 
      dplyr::left_join(cons.clue(tmp), by = "cell") %>% 
      dplyr::mutate(elapsed = median(elapsed),
                    est_k = median(est_k)) %>%
      dplyr::distinct()  ## only need to keep one of the five rows corresponding to the different runs
  }))
}, mc.preschedule = FALSE, mc.cores = ncores))

res_consclue <- plyr::rbind.fill(res_consclue, res_consclue_seurat)

saveRDS(res_consclue, file = outrds)

date()
sessionInfo()


