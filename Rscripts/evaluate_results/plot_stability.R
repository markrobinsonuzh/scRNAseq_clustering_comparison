args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(outrds)

# ------------------------------------------------------------------------
# Stability analysis per method, by computing ARI for each run partition
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(reshape2)
  library(viridis)
  library(ggthemes)
})

# Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

## Compute stability
## nest df
res_summary <- res %>% dplyr::group_by(dataset, method, k) %>% nest() 

res_summary <- res_summary %>% 
  mutate(truenclust = purrr::map_int(data, function(x) {
    length(unique(x$trueclass))
  }))
  
## wide format
cast.map <-  function(x) {
  reshape2::dcast(x, cell ~ run, value.var = "cluster")
}

res_nested <- res_summary %>% mutate(data.wide = purrr::map(data, cast.map)) 

## Help function for computing ARI 
ARI_df <- function(x) {
  stopifnot(class(x) == "data.frame")
  stopifnot(class(x[, 1]) == "character")
  
  x <- select(x, -cell)
  columns <- combn(ncol(x), 2)
  ari.nk <- array(NA, ncol(columns))
  for (i in seq_len(length(ari.nk))) {
    ari.nk[i] <- mclust::adjustedRandIndex(x[, columns[1, i]], x[, columns[2, i]])
  }
  stab <- as.data.frame(cbind(ari.stab = ari.nk, run1 = columns[1, ], run2 = columns[2, ]))
  return(stab)
}

## compute ARI
res_stab.tmp <- res_nested %>% mutate(stability = purrr::map(data.wide, ARI_df)) 

## unnest
res_stab <- res_stab.tmp %>% select(dataset, method, k, stability, truenclust) %>% unnest() %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) 
res_stab$k <- as.integer(res_stab$k)

# ------------------------------------
# Plot stabilities
# ------------------------------------
# methods combined
plots[["stability_allmethods"]] <- 
  ggplot(res_stab %>% dplyr::group_by(dataset, method, filtering, k, truenclust) %>%
           dplyr::summarize(ari.stab = median(ari.stab, na.rm = TRUE)),
         aes(x = k, y = ari.stab, group = method, color = method)) + 
  geom_line(size = 1) + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale +
  facet_grid(filtering ~ dataset, scales = "free_x") +
  ylim(NA, 1) +
  labs(y = "Stability (ARI)", x = "Number of clusters") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "right")

# methods separated
plots[["stability_sepmethods"]] <- 
  ggplot(res_stab,
         aes(x = k, y = ari.stab, group = method, color = method)) + 
  geom_smooth() + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset + filtering ~ method, scales = "free_x") +
  ylim(NA, 1) +
  labs(y = "Stability (ARI)")

# stability at truenclust
plots[["stability_truek"]] <- 
  ggplot(res_stab %>% filter(k == truenclust),
        aes(x = method, y = ari.stab, group = method, color = method)) + 
  geom_boxplot() + 
  theme_bw() +
  manual.scale +
  facet_grid(filtering ~ dataset, scales = "free_x") +
  ylim(NA, 1) +
  labs(y = "Stability (ARI)", title = "k==truenclust") +
  theme(axis.text.x = element_text(size = 10, angle = 90))

# plot heat map on median stability with truenclust
plots[["stability_heatmap_truek"]] <- 
  ggplot(res_stab %>% filter(k == truenclust) %>% 
           dplyr::group_by(filtering, dataset, method, k) %>%
           dplyr::summarise(median.stability = median(ari.stab)),
         aes(x = reorder(method, median.stability, FUN = mean, na.rm = TRUE),
             y = reorder(dataset, median.stability, FUN = mean, na.rm = TRUE),
             fill = median.stability)) +
  geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
  facet_wrap(~ filtering) +
  scale_fill_viridis(name = "Median \nstability \n(ARI)", direction = -1,limits=c(0,1),na.value = "grey") +
  theme_tufte(base_family = "Helvetica") +
  labs(x = NULL, y = NULL, title = "median stability (ARI), k = truenclust") +
  coord_equal() +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.title.align = 0,
        legend.text = element_text(size = 16),
        legend.position = "right",
        legend.key.size = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 20))

pdf(gsub("\\.rds$", "_allmethods.pdf", outrds), width = 20, height = 10)
print(plots[["stability_allmethods"]])
dev.off()

pdf(gsub("\\.rds$", "_truek.pdf", outrds), width = 12, height = 7)
print(plots[["stability_sepmethods"]])
print(plots[["stability_truek"]])
print(plots[["stability_heatmap_truek"]])
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()


