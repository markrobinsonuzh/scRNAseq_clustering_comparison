args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(stabilityrds)
print(entropyrds)
print(diffrds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

stability <- readRDS(stabilityrds)
entropy <- readRDS(entropyrds)
differences <- readRDS(diffrds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 16)
cowplot::plot_grid(
  stability[["stability_heatmap_truek"]] + ggtitle(""), 
  cowplot::plot_grid(
    entropy[["deltanormentropy_at_truth"]],
    differences[["diff_kmax_ktrue"]],
    labels = c("B", "C"), nrow = 1, rel_widths = c(1, 1.3)
  ), 
  labels = c("A", ""), ncol = 1, rel_heights = c(1, 0.8)
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()