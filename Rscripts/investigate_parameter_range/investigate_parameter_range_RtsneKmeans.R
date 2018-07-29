## Investigate effect of the number of tSNE dimensions on timing and performance
## of RtsneKmeans
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

print(scefile)  ## Name of .rds file containing a SingleCellExperiment object
print(k) ## Number of clusters
print(outrds)  ## Name of .rds file where results will be written

method <- "RtsneKmeans"

suppressPackageStartupMessages({
    library(rjson)
    library(SingleCellExperiment)
})
source(paste0("Rscripts/clustering/apply_", method, ".R"))

params_list <- list(
    p_10_2 = list(initial_dims = 50, perplexity = 10, dims = 2),
    p_10_3 = list(initial_dims = 50, perplexity = 10, dims = 3),
    p_10_5 = list(initial_dims = 50, perplexity = 10, dims = 5),
    p_10_7 = list(initial_dims = 50, perplexity = 10, dims = 7),
    p_10_10 = list(initial_dims = 50, perplexity = 10, dims = 10),
    p_30_2 = list(initial_dims = 50, perplexity = 30, dims = 2),
    p_30_3 = list(initial_dims = 50, perplexity = 30, dims = 3),
    p_30_5 = list(initial_dims = 50, perplexity = 30, dims = 5),
    p_30_7 = list(initial_dims = 50, perplexity = 30, dims = 7),
    p_30_10 = list(initial_dims = 50, perplexity = 30, dims = 10),
    p_50_2 = list(initial_dims = 50, perplexity = 50, dims = 2),
    p_50_3 = list(initial_dims = 50, perplexity = 50, dims = 3),
    p_50_5 = list(initial_dims = 50, perplexity = 50, dims = 5),
    p_50_7 = list(initial_dims = 50, perplexity = 50, dims = 7),
    p_50_10 = list(initial_dims = 50, perplexity = 50, dims = 10)
)

## Read data
sce <- readRDS(scefile)

## Run clustering
set.seed(1234)
L <- lapply(params_list, function(params) {
    print(params)
    
    res <- get(paste0("apply_", method))(sce = sce, params = params, k = k)
    
    df <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                     method = method, 
                     cell = names(res$cluster),
                     run = 1,
                     k = k,
                     resolution = NA,
                     cluster = res$cluster,
                     initial_dims = params$initial_dims,
                     perplexity = params$perplexity,
                     dims = params$dims,
                     stringsAsFactors = FALSE, row.names = NULL)
    tm <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                     method = method,
                     run = 1, 
                     k = k,
                     resolution = NA,
                     user.self = res$st[["user.self"]],
                     sys.self = res$st[["sys.self"]],
                     user.child = res$st[["user.child"]],
                     sys.child = res$st[["sys.child"]],
                     elapsed = res$st[["elapsed"]],
                     initial_dims = params$initial_dims,
                     perplexity = params$perplexity,
                     dims = params$dims,
                     stringsAsFactors = FALSE, row.names = NULL)
    list(clusters = df, timing = tm)
})

## Summarize across different runs
assignments <- do.call(rbind, lapply(L, function(w) w$clusters))
timings <- do.call(rbind, lapply(L, function(w) w$timing))

## Add true group for each cell
truth <- data.frame(cell = as.character(rownames(colData(sce))),
                    trueclass = as.character(colData(sce)$phenoid),
                    stringsAsFactors = FALSE)
assignments$trueclass <- truth$trueclass[match(assignments$cell, truth$cell)]

## Save results
saveRDS(list(assignments = assignments, timings = timings), file = outrds)

gc()
sessionInfo()
