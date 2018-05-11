args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

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
## Read clustering results
res <- readRDS( file="output/clustering_summary/clustering_summary.rds" )
pdf("plots/performance/res_Seurat_k_resolution.pdf", width=15, height = 10)

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset,method, run, k, resolution) %>% 
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   timing = median(timing)) %>%
  filter(method=="Seurat") %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# --------------------------------------
# ## plot Seurat k vs. resolution
# --------------------------------------
print(ggplot( res_summary, 
             aes(x = resolution, y = k ) ) +
        geom_smooth() + 
        facet_grid(filtering ~ dataset)+
        theme_bw()+
        scale_color_brewer(palette = "Set3" )
)
print(ggplot(res_summary, 
             aes(x = resolution, y = k ) ) + 
        geom_point() + 
        facet_grid(filtering ~ dataset)+
        theme_bw()+
        scale_color_brewer(palette = "Set3" )
)
# --------------------------------------
# ## Plot timing in dependence of resolution
# --------------------------------------
print(ggplot(res_summary, aes(x = resolution, y = timing, color = method)) + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        scale_color_discrete(name = "") + 
        scale_color_brewer(palette = "Set3")+
        scale_y_log10())

dev.off()

date()
sessionInfo()


