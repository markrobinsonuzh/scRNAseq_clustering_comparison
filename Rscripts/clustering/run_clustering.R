args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(sce)
print(paramfile)
print(method)
print(outrds)

suppressPackageStartupMessages({
  library(rjson)
})
source(paste0("Rscripts/clustering/apply_", method, ".R"))

sce <- readRDS(sce)
params <- fromJSON(file = paramfile)
n_rep <- 10

set.seed(123)
L <- lapply(seq_len(n_rep), function(i) {  ## For each replication
  tmp <- lapply(params$range_clusters, function(k) {  ## For each k
    res <- get(paste0("apply_", method))(sce = sce, params = params, k = k)
    df <- data.frame(method = method, 
                     cell = names(res$cluster),
                     run = i,
                     k = k,
                     cluster = res$cluster,
                     stringsAsFactors = FALSE, row.names = NULL)
    tm <- data.frame(method = method,
                     run = i, 
                     k = k,
                     timing = res$st["user.self"] + res$st["sys.self"] + 
                       res$st["user.child"] + res$st["sys.child"],
                     stringsAsFactors = FALSE, row.names = NULL)
    list(n_cluster = k, clusters = df, timing = tm)
  })  ## End for each k
  
  ## Summarize across different values of k
  assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
  timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
  list(assignments = assignments, timings = timings)
})  ## End for each run

## Summarize across different runs
assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
timings <- do.call(rbind, lapply(L, function(w) w$timings))

saveRDS(list(assignments = assignments, timings = timings), file = outrds)

sessionInfo()
