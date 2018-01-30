## Apply SC3

suppressPackageStartupMessages({
  library(Seurat)
})

apply_Seurat <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      data <- CreateSeuratObject(raw.data = dat, min.cells = params$min.cells,
                                 min.genes = params$min.genes, project = "scRNAseq") 
      data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 1e4)
      data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR)
      data <- ScaleData(object = data)
      data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = FALSE)
      data <- FindClusters(object = data, reduction.type = "pca", save.SNN = TRUE, 
                           dims.use = params$dims.use, k.param = k,
                           resolution = params$resolution, print.output = 0)
      cluster <-  data@ident
    })
    list(st = st, cluster = cluster, est_k = NA)
  }, error = function(e) {
    list(st = NA, cluster = structure(rep(NA, ncol(dat)), names = colnames(dat)),
         est_k = NA)
  })
}
