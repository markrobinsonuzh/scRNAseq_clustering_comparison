args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
filterings <- strsplit(filterings, ",")[[1]]
names(filterings) <- filterings
methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(datasets)
print(filterings)
print(methods)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(mclust)
})

## Read clustering results
res <- do.call(rbind, lapply(datasets, function(d) {
  do.call(rbind, lapply(filterings, function(f) {
    do.call(rbind, lapply(methods, function(m) {
      readRDS(paste0("results/sce_", f, "_", d, "_", m, ".rds"))$timings %>%
        dplyr::select(dataset, method, run, k, resolution, timing)
    }))
  }))
}))

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 6)

## Calculate performance indices for each method and clustering run
res_summary <- res %>% 
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce)

print(ggplot(res_summary, aes(x = k, y = timing, group = method, color = method)) + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset) + 
        scale_color_discrete(name = "") + 
        scale_y_log10())

dev.off()

date()
sessionInfo()
