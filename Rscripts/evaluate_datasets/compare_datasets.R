args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datasets)
print(filtering)
print(outrds)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
})


data <- lapply(datasets, function(d) {
  readRDS(paste0("data/sce_", filtering, "/sce_", filtering, "_", d, ".rds"))
})

comparison <- splatter::compareSCEs(data)

pdf(gsub("\\.rds$", "_comppanel.pdf", outrds), width = 10, height = 12)
makeCompPanel(comparison)
dev.off()

date()
sessionInfo()
