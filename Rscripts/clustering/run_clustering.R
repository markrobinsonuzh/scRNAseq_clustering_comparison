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

set.seed(123)
res <- get(paste0("apply_", method))(sce = sce, params = params, n_rep = 10)

saveRDS(res, file = outrds)

sessionInfo()
