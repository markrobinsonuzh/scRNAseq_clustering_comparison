## Investigate effect of the number of PCs on the Linnorm timing and performance
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scefile)  ## Name of .rds file containing a SingleCellExperiment object
print(k) ## Number of clusters
print(outrds)  ## Name of .rds file where results will be written

method <- "Linnorm"

suppressPackageStartupMessages({
  library(rjson)
  library(SingleCellExperiment)
})
source(paste0("Rscripts/clustering/apply_", method, ".R"))

params_list <- list(
  npc2 = list(minNonZeroPortion = 0.75, num_PC = 2),
  npc3 = list(minNonZeroPortion = 0.75, num_PC = 3),
  npc5 = list(minNonZeroPortion = 0.75, num_PC = 5),
  npc7 = list(minNonZeroPortion = 0.75, num_PC = 7),
  npc10 = list(minNonZeroPortion = 0.75, num_PC = 10)
)

## Read data
sce <- readRDS(scefile)

## Run clustering
set.seed(1234)
L <- lapply(params_list, function(params) {
  res <- get(paste0("apply_", method))(sce = sce, params = params, k = k)
  df <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                   method = method, 
                   cell = names(res$cluster),
                   run = 1,
                   k = k,
                   resolution = NA,
                   cluster = res$cluster,
                   stringsAsFactors = FALSE, row.names = NULL)
  tm <- data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
                   method = method,
                   run = 1, 
                   k = k,
                   resolution = NA,
                   timing = res$st,
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

sessionInfo()
