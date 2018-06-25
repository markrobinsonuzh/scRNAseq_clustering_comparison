args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(consensusrds)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(ggthemes)
})

## Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = consensusrds)

## Initialize list to hold plots
plots <- list()

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset, method, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(consensus, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# --------------------------------------
# ## Calculate performance indices for each method
# --------------------------------------
plots[["consensus_ari"]] <- 
  ggplot(res_summary, 
         aes(x = k, y = ARI, group = method, color = method)) + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  geom_line(size = 1) + 
  theme_bw() +
  scale_color_brewer(palette = "Set3" )  +
  facet_grid(filtering ~ dataset, scales = "free_x")

# --------------------------------------
# ## Heatmap  ARI of truenclust, https://github.com/hrbrmstr/facetedcountryheatmaps
# --------------------------------------
# on true k , ARI, ordered by median
plots[["consensus_heatmap_truek"]] <- 
  ggplot(res_summary %>% dplyr::filter(k == truenclust) %>%
           dplyr::group_by(dataset, filtering, method, k),
         aes(x = reorder(method, ARI, median, na.rm=FALSE), y = dataset, fill = ARI)) +
  geom_tile(color = "white", size = 0.1, na.rm = FALSE) +
  facet_wrap(~ filtering) +
  scale_fill_viridis(name = "ARI", direction = -1 ) +
  theme_tufte(base_family = "Helvetica") +
  labs(x = NULL, y = NULL, title = "consensus ARI, k = truenclust") +
  coord_equal() +
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.title.align = 1,
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_blank())

# on estimated k,  ARI
plots[["consensus_heatmap_estk"]] <- 
  ggplot(res_summary %>% dplyr::filter(k == estnclust) %>%
           dplyr::group_by(dataset, filtering, method, k),
         aes(x = reorder(method, ARI, median, na.rm = FALSE), y = dataset, fill = ARI)) +
  geom_tile(color = "white", size = 0.1) +
  facet_wrap(~ filtering) +
  scale_fill_viridis(name = "ARI", direction = -1, na.value = "grey") +
  theme_tufte(base_family = "Helvetica") +
  labs(x = NULL, y = NULL, title = "Consensus ARI, k = estnclust") +
  coord_equal() +
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        panel.border = element_blank(),
        legend.title = element_text(size = 12),
        legend.title.align = 1,
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_blank())

pdf(gsub("rds$", "df", outrds), width = 20, height = 12)
plots[["consensus_ari"]]
plots[["consensus_heatmap_truek"]]
plots[["consensus_heatmap_est,"]]
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()

