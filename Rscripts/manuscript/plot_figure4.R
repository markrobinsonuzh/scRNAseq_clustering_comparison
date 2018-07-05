args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(consensusrds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

consensus <- readRDS(consensusrds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 10)
cowplot::plot_grid(plotlist = sapply(consensus[1], "[", "ann_tree"))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()