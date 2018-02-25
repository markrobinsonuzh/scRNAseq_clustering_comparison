suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
})

filterHVG <- function(sce, pctkeep) {
  seurat <- CreateSeuratObject(raw.data = counts(sce), meta.data = as.data.frame(colData(sce)),
                               min.cells = 0, min.genes = 0, project = "scRNAseq")
  seurat <- NormalizeData(seurat)
  seurat <- ScaleData(seurat, vars.to.regress = "nUMI", display.progress = FALSE)
  seurat <- FindVariableGenes(seurat, do.plot = FALSE)
  keep <- rownames(head(seurat@hvg.info, n = round(pctkeep/100*nrow(sce))))
  sce[keep, ]
}