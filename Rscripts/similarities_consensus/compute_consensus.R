## -------------------------------------------------------------------------- ##
## Compute consensus using the clue and clusterExperiment packages
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
  library(clusterExperiment)
  library(clue)
  library(multidplyr)
  library(ggplot2)
})

res <- readRDS(file = clusteringsummary)

## consensus function with clue package
cons.clue <- function(cluster, run, cell) {
  m <- tryCatch({
    data.frame(cluster = as.integer(cluster), run = as.integer(run),
               id = c(rep(1:length(unique(cell)), 5))) %>%
      tidyr::spread(key = run, value = cluster) %>%
      dplyr::select(c(-id)) %>% as.matrix, 
    error = function(e) NA
  })
  if (all(is.na(m))) {
    re <- rep(NA, 5*nrow(m))
  } else {
    re <- plyr::alply(m, 2, clue::as.cl_partition)
    re <- lapply(re, function(x) { 
      x <- as.matrix(x$.Data) # some Nas in dataset, replace Nas by zeros
      x[is.na(x)] <- 0
      return(x)
    }) 
    re <- clue::as.cl_ensemble(re)
    re <- clue::cl_consensus(re, method = "SE", control = list(nruns = 50))
    re <- clue::cl_class_ids(re) 
    re <- rep(re, 5) %>% as.character
  }
  return(re) 
}

## Creating cluster 
cluster <- create_cluster(ncores)

## Register function and variable
cluster_assign_value(cluster, "cons.clue", cons.clue)

## Check registered items
cluster_ls(cluster)
## return items
cluster_get(cluster, "cons.clue")

## compute consensus,  Seurat; group by resolution
res_consclue <- res %>% dplyr::filter(!method %in% c("Seurat")) %>%
  partition(dataset, method, k, cluster = cluster) %>%
  do(dplyr::mutate(., consensus.clue = cons.clue(cluster, run, cell))) %>% collect()

## consensus for Seurat by resolution
res_consclue.seurat <- res %>% filter(method %in% c("Seurat")) %>%
  partition(dataset, method, resolution, cluster = cluster) %>%
  do(dplyr::mutate(., consensus.clue = cons.clue(cluster, run, cell))) %>% collect()
  
res_consclue2 <- bind_rows(res_consclue, res_consclue.seurat)

saveRDS(res_consclue2, file = outrds)

## Unregister function 
cluster_rm(cluster, c("cons.clue"))
#_____________________________________________

date()
sessionInfo()


