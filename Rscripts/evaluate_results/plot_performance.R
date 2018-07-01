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
  library(ggalluvial)
})

## Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

## Compute ARI, true number of clusters, estimated number of clusters, 
## elapsed time
res_summary <- res %>% 
  dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   elapsed = median(elapsed)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

## ARI plots
## Median ARI vs k
plots[["median_ari_vs_k"]] <- 
  ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
           dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust)) %>%
           dplyr::ungroup(),
         aes(x = k, y = medianARI, group = method, color = method)) + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  geom_line(size = 1) + 
  theme_bw() +
  manual.scale +
  facet_grid(filtering ~ dataset, scales = "free_x") +
  labs(title = "", x = "Number of clusters", y = "Median Adjusted Rand Index") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

## Heatmap of median ARI at true k
plots[["median_ari_heatmap_truek"]] <- 
  ggplot(res_summary %>% dplyr::filter(k == truenclust) %>%
          dplyr::group_by(dataset, filtering, method, k) %>%
          dplyr::summarize(medianARI = median(ARI)),
         aes(x = reorder(method, medianARI, FUN = mean, na.rm = TRUE), 
             y = reorder(dataset, medianARI, FUN = mean, na.rm = TRUE), 
             fill = medianARI)) +
  geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
  facet_wrap(~ filtering) +
  scale_fill_viridis(name = "Median ARI", direction = -1) +
  theme_tufte(base_family = "Helvetica") +
  labs(x = NULL, y = NULL, title = "") +
  coord_equal() +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.title.align = 1,
        legend.text = element_text(size = 16),
        legend.position = "right",
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 20))

plots[["median_ari_heatmap_estnclust"]] <- 
  ggplot(res_summary %>% dplyr::filter(k == estnclust) %>%
           dplyr::group_by(dataset, filtering, method, k) %>%
           dplyr::summarize(medianARI = median(ARI)),
         aes(x = reorder(method, medianARI, FUN = mean, na.rm = TRUE), 
             y = reorder(dataset, medianARI, FUN = mean, na.rm = TRUE), 
             fill = medianARI)) +
  geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
  facet_wrap(~ filtering) +
  scale_fill_viridis(name = "Median ARI", direction = -1) +
  theme_tufte(base_family = "Helvetica") +
  labs(x = NULL, y = NULL, title = "") +
  coord_equal() +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.title.align = 1,
        legend.text = element_text(size = 16),
        legend.position = "right",
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 20))

pdf(gsub("\\.rds$", "_medianARIvsk.pdf", outrds), width = 20, height = 15)
print(plots[["median_ari_vs_k"]])
dev.off()

pdf(gsub("\\.rds$", "_medianARIheatmap_truek.pdf", outrds), width = 20, height = 9)
print(plots[["median_ari_heatmap_truek"]])
dev.off()

pdf(gsub("\\.rds$", "_medianARIheatmap_estk.pdf", outrds), width = 15, height = 9)
print(plots[["median_ari_heatmap_estnclust"]])
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()

