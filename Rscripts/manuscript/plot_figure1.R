args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(performancerds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

performance <- readRDS(performancerds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 9)
performance[["median_ari_heatmap_truek"]]
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()