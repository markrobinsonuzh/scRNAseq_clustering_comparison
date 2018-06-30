args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(cowplot)
})

ensemble <- readRDS("plots/ensemble/ensemble_vs_individual.rds")

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 15)
cowplot::plot_grid(
  ensemble[["ensembl_vs_bestworst_truek"]],
  ensemble[["ensembl_vs_first_truek"]],
    labels = c("A", "B"), nrow = 1
)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()