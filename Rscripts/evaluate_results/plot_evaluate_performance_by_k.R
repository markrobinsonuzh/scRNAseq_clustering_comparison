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

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 15)

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% 
  dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   elapsed = median(elapsed)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# --------------------------------------
# ## Calculate performance indices for each method and clustering run
# --------------------------------------
# by k
print(ggplot(res_summary,
             aes(x = k, y = ARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_smooth() + 
        theme_bw() +
        manual.scale +
        facet_grid(filtering ~ dataset, scales = "free_x") +
        labs(title = "ARI by k")
)

# by median ARI
print(ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
               dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust)) %>%
               dplyr::ungroup(),
             aes(x = k, y = medianARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_line(size = 1) + 
        theme_bw() +
        manual.scale +
        facet_grid(filtering ~ dataset, scales = "free_x") +
        labs(title = "median ARI by k")
)

# --------------------------------------
# Plot timing one boxplot per dataset, over all ks and runs as they are similar
# --------------------------------------
print(ggplot(res_summary, aes(x = method, y = elapsed, group = method, color = method)) + 
        geom_boxplot() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        scale_y_log10() +
        theme_bw() +
        manual.scale +
        theme(axis.text.x = element_text(size = rel(0.8), angle = 90, hjust = 1, vjust = 1))+
        labs(title = "Runtime per method")
)

# --------------------------------------
# Plot timing 
# --------------------------------------
# by k
print(ggplot(res_summary, aes(x = k, y = elapsed, group = method, color = method)) + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        manual.scale +
        scale_y_log10() +
        labs(title = "Runtime by k")
)
  
# time by rank, k = truenclust, NAs removed
print(ggplot(data = res_summary %>% dplyr::group_by(dataset, method) %>% 
               filter(k == truenclust) %>% 
               dplyr::summarize(med.elapsed = median(elapsed, na.rm = TRUE)) %>% 
               filter(!is.na(med.elapsed)),
             aes(x = dataset, weight = med.elapsed, stratum = dataset, 
                 alluvium = method, na.rm = TRUE)) +
        geom_alluvium(aes(fill = method, colour = method), alpha = 0.75, 
                      decreasing = TRUE, na.rm = TRUE, reverse = TRUE) +
        scale_fill_brewer(type = "qual", palette = "Set3") +
        scale_color_brewer(type = "qual", palette = "Set3") +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(size = 10, angle = 90))+
        labs(title = "Runtime by rank")
)

# time by rank, k==truenclust
print(ggplot(data = res_summary %>% dplyr::group_by(dataset, method) %>% 
               filter(k == truenclust, !method %in% c("RaceID", "CIDR")) %>% 
               dplyr::summarize(med.elapsed = median(elapsed, na.rm = TRUE)) %>%
               dplyr::mutate(rank = (rank(med.elapsed, na.last = TRUE))) %>%
               dplyr::mutate(rank = factor(rank, levels = c(1:12))),
             aes(x = dataset, y = rank, stratum = rank, alluvium = method, fill = method)) +
        geom_flow(aes(fill = method)) +
        geom_stratum(aes(fill = method)) +
        geom_lode(aes(fill = method)) +
        scale_fill_brewer(type = "qual", palette = "Set3") +
        scale_color_brewer(type = "qual", palette = "Set3") +
        geom_text(stat = "stratum", label.strata = TRUE) +
        theme(axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size = 10, angle = 90)) +
        labs(title = "Runtime by rank")
)

# ---------------
# normalized by median RtsneKmeans , time combined
median.tsne <- res_summary%>%select(filtering, dataset, method, k, run, truenclust, elapsed) %>%
  dplyr::group_by(dataset, filtering, k) %>%
  filter(method == "RtsneKmeans") %>% dplyr::summarise(med.t = median(elapsed))%>%
  dplyr::ungroup()
res.time <- res_summary %>% 
  group_by(filtering, dataset, method, k) %>% 
  dplyr::summarise(median.elapsed = median(elapsed))%>%
  dplyr::ungroup()

res.time <- full_join(res.time, median.tsne, by = c("dataset", "filtering", "k")) %>% 
  dplyr::mutate(norm.time = median.elapsed/med.t)
print(ggplot(res.time, 
             aes(x = reorder(method, norm.time, FUN = median, order = TRUE, na.rm = TRUE) , 
                 y = norm.time, group = method, color = method)) +
        manual.scale +
        geom_boxplot() +
        scale_y_log10() +
        labs(y = "normalised time", x = "")+
        theme_bw() +
        labs(title = "Runtime, normalized by RtsneKmeans", size = 16)+
        theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1)) +
        theme(legend.position = "none") +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.text = element_text(size = 15)) +
        theme(legend.position = "none") +
        theme(strip.text = element_text(size = 16))
)

# --------------------------------------
# ## Heatmap median ARI of truenclust, from https://github.com/hrbrmstr/facetedcountryheatmaps
# --------------------------------------
# on true k , median ARI, ordered by median
print(res_summary %>% dplyr::filter(k == truenclust) %>%
        dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI)) %>% 
        ggplot2::ggplot(
          aes(x = reorder(method,  medianARI, FUN = mean, na.rm = TRUE), 
              y = reorder(dataset, medianARI, FUN = mean, na.rm = TRUE), 
              fill = medianARI)) +
          geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
          facet_wrap(~ filtering) +
          scale_fill_viridis(name = "medianARI", direction = -1) +
          theme_tufte(base_family = "Helvetica") +
          labs(x = NULL, y = NULL, title = "median ARI, k = truenclust") +
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
)

# on estimated k, median ARI
print(res_summary %>% dplyr::filter(k == estnclust) %>%
        dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI)) %>%
        ggplot(
          aes(x = reorder(method, medianARI, FUN = mean, na.rm = TRUE),
              y = reorder(dataset, medianARI, FUN = mean, na.rm = TRUE),
              fill = medianARI)) +
        geom_tile(color = "white", size = 0.1) +
        facet_wrap(~ filtering) +
        scale_fill_viridis(name = "medianARI", direction = -1, na.value = "grey") +
        theme_tufte(base_family = "Helvetica") +
        labs(x = NULL, y = NULL, title = "median ARI, k = estnclust") +
        coord_equal() +
        theme(axis.text.x = element_text(size = 18, angle = 90),
              axis.text.y = element_text(size = 16),
              panel.border = element_blank(),
              legend.title = element_text(size = 16),
              legend.title.align = 1,
              legend.text = element_text(size = 16),
              legend.position = "right",
              legend.key.size = unit(2, "cm"),
              legend.key.width = unit(0.5, "cm"),
              axis.ticks = element_blank(),
              strip.text = element_text(size = 18))
)

## Excerpts
# excerpt Range of k,  Koh and Zhengmix4eq filtered10Expr
print(ggplot(res_summary %>% filter(dataset %in% c("Koh", "Zhengmix4eq"), 
                                    filtering %in% c("filteredExpr10")),
             aes(x = k, y = ARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_smooth() + 
        theme_bw() +
        manual.scale +
        facet_grid(~ dataset, scales = "free_x") +
        labs(title = "ARI by k") +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              strip.text = element_text(size = 16))
)

# timings, excerpt for Ko and Zheng
print(ggplot(res_summary %>% filter(dataset %in% c("Koh", "Zhengmix4eq"), 
                                    filtering %in% c("filteredExpr10")), 
             aes(x = reorder(method, elapsed, FUN = mean, na.rm = TRUE), 
                 y = elapsed, group = method, color = method)) + 
        geom_boxplot() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        scale_y_log10() +
        theme_bw() +
        manual.scale +
        labs(title="Runtime per method for Koh and Zheng", y = "timing (seconds)", x = "method") +
        theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 1),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.position = "none",
              strip.text = element_text(size = 16))
)

print(ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
               dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust)) %>%
               filter(dataset %in% c("Kumar", "Trapnell"), filtering %in% c("filteredExpr10")) %>%
               dplyr::ungroup(),
             aes(x = k, y = medianARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_line(size = 1) + 
        theme_bw() +
        manual.scale +
        facet_grid(filtering ~ dataset, scales = "free_x") +
        labs(title = "median ARI by k for Kumar and Trapnell") +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              strip.text = element_text(size = 16))
)

dev.off()

date()
sessionInfo()
