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
})

sce <- readRDS(scefile)

if (method == "Expr") {
  exprsn <- rowMeans(logcounts(sce))
  keep <- order(exprsn, decreasing = TRUE)[seq_len(pctkeep/100*length(exprsn))]
  sce_filt <- sce[keep, ]
}


sce_filt <- computeSumFactors(sce_filt, sizes = c(20, 40, 80, 150))
sce_filt <- normalise(sce_filt, exprs_values = "counts", return_log = TRUE, 
                      return_norm_as_exprs = TRUE)
sce_filt <- normalise(sce_filt, exprs_values = "counts", return_log = FALSE, 
                      return_norm_as_exprs = FALSE)

saveRDS(sce_filt, file = outrds)
