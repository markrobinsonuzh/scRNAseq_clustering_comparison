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

# ------------------------------------------------
# Venn diagrams of genes for the filtered datasets
# ------------------------------------------------
suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(gridExtra)
  library(SingleCellExperiment)
  library(cowplot)
  library(UpSetR)
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
                           cat.dist = c(0.05, 0.05, 0.05),
                           print.mode = "percent"
                           )
  })

## Reformat gList to gTree
plots <-  lapply(venn_diagrams, function(x) gTree(children = x))

pdf(gsub("rds$", "pdf", outrds), width = 15, height = 15)
cowplot::plot_grid(plotlist = plots, ncol = 3, labels = names(plots))
dev.off()

## UpSet plots
plots_upset <- lapply(gene_ids_by_ds, function(w) {
  allgenes <- unique(unlist(w))
  dat <- do.call(cbind, lapply(w, function(i) allgenes %in% i))
  upset(data.frame(dat + 0), nintersects = NA, order.by = "degree", 
        empty.intersections = TRUE)
  grid.edit('arrange', name = 'arrange2')
  grid.grab()
})
pdf(gsub("\\.rds$", "_upset.pdf", outrds), height = 15, width = 15)
cowplot::plot_grid(plotlist = plots_upset, ncol = 3, labels = names(plots), 
                   label_x = 0.05, hjust = 0)
dev.off()


saveRDS(list(plots = plots, plots_upset), file = outrds)
date()
sessionInfo()

