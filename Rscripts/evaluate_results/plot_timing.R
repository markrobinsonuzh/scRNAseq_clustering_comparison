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

## Elapsed time, one boxplot per dataset, over all ks and runs as they are similar
plots[["time_boxplot_perds"]] <- 
  ggplot(res_summary, aes(x = method, y = elapsed, group = method, color = method)) + 
  geom_boxplot() + 
  facet_grid(filtering ~ dataset, scales = "free") + 
  scale_y_log10() +
  theme_bw() +
  manual.scale +
  labs(title = "", x = "", y = "Elapsed time (s)") + 
  theme(axis.text.x = element_text(size = rel(0.8), angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

## Elapsed time, normalized by median time for RtsneKmeans
median.tsne <- res_summary %>% 
  dplyr::select(filtering, dataset, method, k, run, truenclust, elapsed) %>%
  dplyr::group_by(dataset, filtering, k) %>%
  filter(method == "RtsneKmeans") %>% dplyr::summarise(med.t = median(elapsed)) %>%
  dplyr::ungroup()
res.time <- res_summary %>% 
  group_by(filtering, dataset, method, k) %>% 
  dplyr::summarise(median.elapsed = median(elapsed))%>%
  dplyr::ungroup()
res.time <- dplyr::full_join(res.time, median.tsne, by = c("dataset", "filtering", "k")) %>% 
  dplyr::mutate(norm.time = median.elapsed/med.t)

plots[["time_normalized_by_tsne"]] <- 
  ggplot(res.time, 
         aes(x = reorder(method, norm.time, FUN = median, order = TRUE, na.rm = TRUE) , 
             y = norm.time, group = method, color = method)) +
  manual.scale +
  geom_boxplot(size = 1.1) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "", y = "Run time, normalized by RtsneKmeans", size = 16)+
  theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16))

plots[["time_by_k"]] <- 
  ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
           dplyr::summarize(medianelapsed = median(elapsed)) %>%
           dplyr::ungroup(), 
         aes(x = k, y = medianelapsed, group = method, color = method)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  facet_grid(filtering ~ dataset, scales = "free") + 
  manual.scale +
  scale_y_log10() +
  labs(title = "", y = "Elapsed time (s)", x = "Number of clusters") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

pdf(gsub("\\.rds$", "_boxplot_perds.pdf", outrds), width = 20, height = 15)
print(plots[["time_boxplot_perds"]])
dev.off()

pdf(gsub("\\.rds$", "_normalizedtime.pdf", outrds), width = 20, height = 15)
print(plots[["time_normalized_by_tsne"]])
dev.off()

pdf(gsub("\\.rds$", "_time_vs_k.pdf", outrds), width = 20, height = 15)
print(plots[["time_by_k"]])
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()
