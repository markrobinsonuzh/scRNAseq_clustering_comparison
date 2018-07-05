args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(ensemblerds)
print(outrds)

suppressPackageStartupMessages({
  library(cowplot)
})

ensemble <- readRDS(ensemblerds)

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 10)
cowplot::plot_grid(
  ensemble[["ensembl_vs_bestworst_truek"]],
  ensemble[["ensembl_vs_first_truek"]],
    labels = c("A", "B"), nrow = 1, label_size = 20
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()