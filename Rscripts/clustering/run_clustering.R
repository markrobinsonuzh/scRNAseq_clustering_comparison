args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scefile)  ## Name of .rds file containing a SingleCellExperiment object
print(method)  ## Clustering method
print(outrds)  ## Name of .rds file where results will be written

suppressPackageStartupMessages({
  library(rjson)
  library(SingleCellExperiment)
})
source(paste0("Rscripts/clustering/apply_", method, ".R"))

## Load parameter files. General dataset and method parameters as well as
## dataset/method-specific parameters
params <- c(fromJSON(file = paste0("parameter_settings/", 
                                   gsub("\\.rds$", "_", basename(scefile)), method, ".json")),
            fromJSON(file = paste0("parameter_settings/", method, ".json")), 
            fromJSON(file = paste0("parameter_settings/", 
                                   gsub("\\.rds$", ".json", basename(scefile))))
            )
## If any parameter is repeated, take the most specific
if (any(duplicated(names(params)))) {
  warning("Possibly conflicting settings")
  params <- params[!duplicated(names(params))]
}
print(params)

## Set number of times to run clustering for each k
n_rep <- 5

## Read data
sce <- readRDS(scefile)

## Run clustering
set.seed(1234)
L <- lapply(seq_len(n_rep), function(i) {  ## For each run
  cat(paste0("run = ", i, "\n"))
  if (method == "Seurat") {
    tmp <- lapply(params$range_resolutions, function(resolution) {  ## For each resolution
      cat(paste0("resolution = ", resolution, "\n"))
      ## Run clustering
      res <- get(paste0("apply_", method))(sce = sce, params = params, resolution = resolution)
      
      ## Put output in data frame
      df <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                       method = method, 
                       cell = names(res$cluster),
                       run = i,
                       k = length(unique(res$cluster)),
                       resolution = resolution,
                       cluster = res$cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                       method = method,
                       run = i, 
                       k = length(unique(res$cluster)),
                       resolution = resolution,
                       user.self = res$st[["user.self"]],
                       sys.self = res$st[["sys.self"]],
                       user.child = res$st[["user.child"]],
                       sys.child = res$st[["sys.child"]],
                       elapsed = res$st[["elapsed"]],
                       stringsAsFactors = FALSE, row.names = NULL)
      kest <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                         method = method,
                         run = i, 
                         k = length(unique(res$cluster)),
                         resolution = resolution,
                         est_k = res$est_k,
                         stringsAsFactors = FALSE, row.names = NULL)
      list(clusters = df, timing = tm, kest = kest)
    })  ## End for each resolution
  } else {
    tmp <- lapply(params$range_clusters, function(k) {  ## For each k
      cat(paste0("k = ", k, "\n"))
      ## Run clustering
      res <- get(paste0("apply_", method))(sce = sce, params = params, k = k)
      
      ## Put output in data frame
      df <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                       method = method, 
                       cell = names(res$cluster),
                       run = i,
                       k = k,
                       resolution = NA,
                       cluster = res$cluster,
                       stringsAsFactors = FALSE, row.names = NULL)
      tm <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                       method = method,
                       run = i, 
                       k = k,
                       resolution = NA,
                       user.self = res$st[["user.self"]],
                       sys.self = res$st[["sys.self"]],
                       user.child = res$st[["user.child"]],
                       sys.child = res$st[["sys.child"]],
                       elapsed = res$st[["elapsed"]],
                       stringsAsFactors = FALSE, row.names = NULL)
      kest <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                         method = method,
                         run = i, 
                         k = k,
                         resolution = NA,
                         est_k = res$est_k,
                         stringsAsFactors = FALSE, row.names = NULL)
      list(clusters = df, timing = tm, kest = kest)
    })  ## End for each k
  }
  
  ## Summarize across different values of k
  assignments <- do.call(rbind, lapply(tmp, function(w) w$clusters))
  timings <- do.call(rbind, lapply(tmp, function(w) w$timing))
  k_estimates <- do.call(rbind, lapply(tmp, function(w) w$kest))
  list(assignments = assignments, timings = timings, k_estimates = k_estimates)
})  ## End for each run

## Summarize across different runs
assignments <- do.call(rbind, lapply(L, function(w) w$assignments))
timings <- do.call(rbind, lapply(L, function(w) w$timings))
k_estimates <- do.call(rbind, lapply(L, function(w) w$k_estimates))

## Add true group for each cell
truth <- data.frame(cell = as.character(rownames(colData(sce))),
                    trueclass = as.character(colData(sce)$phenoid),
                    stringsAsFactors = FALSE)
assignments$trueclass <- truth$trueclass[match(assignments$cell, truth$cell)]

## Save results
saveRDS(list(assignments = assignments, timings = timings,
             k_estimates = k_estimates), file = outrds)

gc()
sessionInfo()
