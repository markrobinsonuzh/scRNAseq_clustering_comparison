args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}
# ------------------------------------------------------------------------
# Comparing filterings by ARI , for all k 
# ------------------------------------------------------------------------

print(clusteringsummary)
print(outrds)

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(reshape2)
  library(RColorBrewer)
})

res <-  readRDS(file = clusteringsummary)

# ------------------------------------
# Compute the stability, based on ARI
# ------------------------------------
res_sep <- res %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

## nest df
res_summary <- res_sep %>% dplyr::group_by(dataset, method, k) %>% nest() 

## extract truenclust
res_summary <- res_summary %>% mutate(truenclust = purrr::map_int(data, function(x) {
  length(unique(x$trueclass))
}))

## wide format
cast.map <- function(x){
  reshape2::dcast(x, cell ~ run + filtering, value.var = "cluster")
}

res_nested <- res_summary %>% mutate(data.wide = purrr::map(data, cast.map)) 

## function for computing ARI 
ARI_df <- function(x){
  stopifnot(class(x) == "data.frame")
  stopifnot(class(x[, 1]) == "character")
  
  x <- dplyr::select(x, -cell)
  cols <- combn(ncol(x), 2)
  ari.nk <- matrix(NA, ncol(cols), 5)
  
  for (i in 1:ncol(cols)){
    ari.nk[i, 1] <- mclust::adjustedRandIndex(x[, cols[1, i]], x[, cols[2, i]])
    ari.nk[i, 2] <- cols[1, i]
    ari.nk[i, 3] <- cols[2, i]
    ari.nk[i, 4] <- colnames(x[1:ncol(x)])[cols[1, i]]
    ari.nk[i, 5] <- colnames(x[1:ncol(x)])[cols[2, i]]
  }
  
  ari.nk <- as.data.frame(ari.nk)
  names(ari.nk) <- c("ARI", "run1", "run2", "filtering_1", "filtering_2")
  return(ari.nk)
}

## compute ARI
res_stab.tmp <- res_nested %>% mutate(stability = purrr::map(data.wide, ARI_df)) 

## unnest
res_stab <- res_stab.tmp %>% select(dataset, method, k, stability, truenclust ) %>% unnest()  
res_stab <- res_stab %>% tidyr::separate(filtering_1, sep = "_", into = c("run1", "filtering1")) %>%
  tidyr::separate(filtering_2, sep = "_", into = c("run2", "filtering2")) %>%
  tidyr::unite(filtering1, filtering2, col = "filtering") %>%
  dplyr::mutate(filtering = dplyr::recode_factor(
    filtering,
    filteredHVG10_filteredExpr10 = "filteredExpr10_filteredHVG10",
    filteredM3Drop10_filteredExpr10 = "filteredExpr10_filteredM3Drop10",
    filteredHVG10_filteredM3Drop10 = "filteredM3Drop10_filteredHVG10")) %>%
  dplyr::ungroup()

res_stab$ARI <- as.numeric(res_stab$ARI)

# ------------------------------------
# Plot stability by k
# ------------------------------------
# load colors
# color set , from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
colors <- c("red", "blue", "green", 
            "orange", "pink", "brown")

filters <- unique(res_stab$filtering)
# manual scale for ggplot
names(colors) <- filters
manual.scale <- scale_colour_manual(name = "method", values = colors)

pdf("plots/performance/plot_ARI_between_filterings.pdf", width = 20)
## ARI between filterings at truenclust, by method and dataset
ggplot(res_stab %>% dplyr::group_by(dataset, method, k, filtering) %>%
         dplyr::filter(k == truenclust),
       aes(x = filtering, y = ARI, group = filtering, color = filtering)) + 
  geom_boxplot(na.rm = TRUE) + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset ~ method) +
  labs(y = "ARI", title = "ARI between filterings, k = truenclust") +
  theme(axis.text.x = element_blank())

## ARI between filterings by k , by method and dataset
ggplot(res_stab %>% dplyr::group_by(dataset, method, k, filtering),
       aes(x = k, y = ARI, color = filtering)) + 
  geom_smooth() + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset ~ method , scales = "free") +
  ylim(NA, 1) +
  labs(y = "ARI", title = "ARI between filterings, by k")

## ARI between filterings at truenclust, by method 
ggplot(res_stab %>% dplyr::group_by(dataset, method, k, filtering) %>%
         dplyr::filter(k == truenclust),
       aes(x = filtering, y = ARI, group = filtering, color = filtering)) + 
  geom_boxplot(na.rm = TRUE) +
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset ~ method ) +
  labs(y = "ARI", title = "ARI between filterings, by method") +
  theme(axis.text.x = element_blank())

ggplot(res_stab,
       aes(x = k, y = ARI, group = method, color = method)) + 
  geom_smooth() + 
  theme_bw() +
  facet_grid(filtering ~ dataset, scales = "free_x") +
  ylim(NA, 1) +
  labs(y = "ARI", title = "ARI by method")
dev.off()

## for Seurat 

res_summary_seurat <- res %>% dplyr::filter(method == "Seurat") %>%
  dplyr::group_by(dataset, method, resolution) %>% nest() 

res_summary_seurat <- res_summary_seurat %>% 
  dplyr::mutate(truenclust = purrr::map_int(data, function(x) {
  length(unique(x$trueclass))
}))

## compute ARI
res_nested <- res_summary_seurat %>% dplyr::mutate(data.wide = purrr::map(data, cast.map)) 
res_stab.tmp <- res_nested %>% dplyr::mutate(stability = purrr::map(data.wide, ARI_df)) 

## unnest
res_stab <- res_stab.tmp %>% 
  dplyr::select(dataset, method, resolution, stability, truenclust) %>%
  unnest() %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce)

ggplot(res_stab,
       aes(x = resolution, y = ari.stab, group = method, color = method)) + 
  geom_line() + 
  theme_bw() +
  manual.scale +
  facet_grid(filtering ~ dataset) +
  ylim(NA, 1) +
  labs(y = "Stability (ARI)")


