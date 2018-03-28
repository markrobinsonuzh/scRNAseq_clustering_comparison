suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(M3Drop)
})

filterM3Drop <- function(sce, pctkeep) {
  cpm <- calculateCPM(sce)
  cpm <- cpm[rowMeans(cpm) > 0.05, ]
  res <- M3DropDifferentialExpression(expr_mat = cpm, 
                                      mt_method = "none", mt_threshold = 1, 
                                      suppress.plot = TRUE)
  res <- res[order(res$p.value), ]
  stopifnot(nrow(res) > round(pctkeep/100*nrow(sce)))
  keep <- rownames(head(res, n = round(pctkeep/100*nrow(sce))))
  sce[keep, ]
}