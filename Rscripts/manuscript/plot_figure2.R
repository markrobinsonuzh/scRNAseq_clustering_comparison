args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(stabilityrds)
print(entropyrds)
print(diffrds)
print(timerds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

stability <- readRDS(stabilityrds)
entropy <- readRDS(entropyrds)
differences <- readRDS(diffrds)
timing <- readRDS(timerds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 15)
cowplot::plot_grid(
  stability[["stability_heatmap_truek"]], 
  cowplot::plot_grid(
    entropy[["normentropy_allds_allk"]],
    differences[["diff_kmax_ktrue"]],
    timing[["time_normalized_by_tsne"]], 
    labels = c("B", "C", "D"), nrow = 1
  ), 
  labels = c("A", ""), ncol = 1
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()