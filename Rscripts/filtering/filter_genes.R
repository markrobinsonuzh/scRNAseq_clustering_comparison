args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scefile)  ## Name of .rds file containing a SingleCellExperiment object
print(method)  ## Filtering method
print(pctkeep)  ## Percentage of genes to keep
print(outrds)  ## Name of .rds file where results will be written

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})
source(paste0("Rscripts/filtering/filter", method, ".R"))

sce <- readRDS(scefile)

sce_filt <- get(paste0("filter", method))(sce = sce, pctkeep = pctkeep)

sce_filt <- calculateQCMetrics(sce_filt)
sce_filt <- computeSumFactors(sce_filt, sizes = pmin(ncol(sce_filt), seq(20, 120, 20)), min.mean = 0.5)
table(sizeFactors(sce_filt) < 0)
sce_filt <- sce_filt[, sizeFactors(sce_filt) > 0]
sce_filt <- normalise(sce_filt, exprs_values = "counts", return_log = TRUE, 
                      return_norm_as_exprs = TRUE)
sce_filt <- normalise(sce_filt, exprs_values = "counts", return_log = FALSE, 
                      return_norm_as_exprs = FALSE)
sce_filt <- runTSNE(sce_filt, exprs_values = "logcounts", perplexity = 10)

saveRDS(sce_filt, file = outrds)

date()
sessionInfo()
