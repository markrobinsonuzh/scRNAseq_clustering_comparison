args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

filterings <- strsplit(filterings, ",")[[1]]
datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(filterings)
print(pctkeep)
print(datasets)
print(datadir)
print(outrds)

#------------------------------------------------
# Venn diagrams of genes for the filtered datasets
#------------------------------------------------
suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(gridExtra)
  library(SingleCellExperiment)
})

datasets_filtered <- c(outer(paste0(datadir, "/sce_filtered", filterings, pctkeep, "/sce_filtered", 
                                    filterings, pctkeep), 
                             paste0(datasets, ".rds"),
                             FUN = paste, sep = "_"))
names(datasets_filtered) <- gsub("\\.rds", "", basename(datasets_filtered))

gene_ids <- lapply(datasets_filtered, function(x) {
  rownames(readRDS(x))
})

gene_ids_by_ds <- lapply(
  lapply(datasets, function(ds) purrr::keep(gene_ids, grepl(paste0("_", ds, "$"), names(gene_ids)))), 
  function(x) {
    names(x) <- stringr::str_extract(names(x), paste(filterings, collapse = "|"))
    x
  })

venn_diagrams <- lapply(gene_ids_by_ds ,function(x){
 VennDiagram::venn.diagram(x, filename = NULL, 
                           fill = c("darkmagenta", "darkblue", "darkred"),
                           alpha = c(0.3, 0.3, 0.3),
                           cex = 0.8, cat.fontface = 1,
                           margin = 0.22, cat.cex = c(1, 1, 1),
                           euler.d = TRUE,
                           ext.text = TRUE,
                           lwd = 1,
                           cat.default.pos = "outer",
                           cat.dist = c(0.05, 0.05, 0.05)
                           )
  })

## Reformat gList to gTree
venn_tree <-  lapply(venn_diagrams, function(x)gTree(children=x))

## Plot grid
plots <- cowplot::plot_grid(plotlist = venn_tree, ncol = 3, labels = names(venn_tree))

pdf(gsub("rds$", "pdf", outrds), width = 15, height = 10)
cowplot::plot_grid(plots)
dev.off()

date()
sessionInfo()

