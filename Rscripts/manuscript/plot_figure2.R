args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(performancerds)
print(timerds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

performance <- readRDS(performancerds)
timing <- readRDS(timerds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 10)
cowplot::plot_grid(
  timing[["time_normalized_by_tsne"]], 
  performance[["scatter_time_vs_ari_truek"]], 
  labels = c("A", "B"), nrow = 1, rel_widths = c(1, 3)
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()