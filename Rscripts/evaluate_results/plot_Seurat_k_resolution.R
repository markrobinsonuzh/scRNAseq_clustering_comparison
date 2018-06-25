args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(ggthemes)
  library(viridis)
})

# Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

## Compute ARI, true number of clusters, estimated number of clusters, 
## elapsed time
res_summary <- res %>% dplyr::group_by(dataset, method, run, k, resolution) %>% 
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   elapsed = median(elapsed)) %>%
  filter(method == "Seurat") %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# --------------------------------------
# Plot Seurat k vs. resolution
# --------------------------------------
plots[["seurat_k_vs_res"]] <- 
  ggplot(res_summary, 
         aes(x = resolution, y = k)) + 
  geom_point() + 
  facet_grid(filtering ~ dataset) +
  theme_bw() +
  scale_color_brewer(palette = "Set3")

# --------------------------------------
# Plot timing in dependence of resolution
# --------------------------------------
plots[["seurat_res_vs_time"]] <- 
  ggplot(res_summary, aes(x = resolution, y = elapsed, color = method)) + 
  geom_smooth() + 
  facet_grid(filtering ~ dataset, scales = "free") + 
  scale_color_discrete(name = "") + 
  scale_color_brewer(palette = "Set3") +
  scale_y_log10()

pdf(gsub("rds$", "pdf", outrds), width = 15, height = 10)
print(plots[["seurat_k_vs_res"]])
print(plots[["seurat_res_vs_time"]])
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()


