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

## Read clustering results
res <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

# ------------------------------------
# Compute the stability, based on ARI
# ------------------------------------
res_sep <- res %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

## nest df
res_summary <- res_sep %>% dplyr::group_by(dataset, method, k) %>% nest() 

## extract truenclust
res_summary <- res_summary %>% 
  mutate(truenclust = purrr::map_int(data, function(x) {
    length(unique(x$trueclass))
  }))

## wide format
cast.map <- function(x) {
  reshape2::dcast(x, cell ~ run + filtering, value.var = "cluster")
}

res_nested <- res_summary %>% mutate(data.wide = purrr::map(data, cast.map)) 

## function for computing ARI 
ARI_df <- function(x) {
  stopifnot(class(x) == "data.frame")
  stopifnot(class(x$cell) == "character")
  
  x <- dplyr::select(x, -cell)
  cols <- combn(ncol(x), 2)

  do.call(rbind, lapply(seq_len(ncol(cols)), function(i) {
    data.frame(ARI = mclust::adjustedRandIndex(x[, cols[1, i]], x[, cols[2, i]]),
               run1 = cols[1, i],
               run2 = cols[2, i],
               filtering1 = colnames(x)[cols[1, i]],
               filtering2 = colnames(x)[cols[2, i]],
               stringsAsFactors = FALSE)
  }))
}

## compute ARI
res_stab.tmp <- res_nested %>% mutate(stability = purrr::map(data.wide, ARI_df)) 

## unnest
res_stab <- res_stab.tmp %>% select(dataset, method, k, stability, truenclust) %>% unnest()  
res_stab <- res_stab %>% tidyr::separate(filtering1, sep = "_", into = c("run1", "filtering1")) %>%
  tidyr::separate(filtering2, sep = "_", into = c("run2", "filtering2")) %>%
  tidyr::unite(filtering1, filtering2, col = "filtering") %>%
  dplyr::mutate(filtering = dplyr::recode_factor(
    filtering,
    filteredHVG10_filteredExpr10 = "filteredExpr10_filteredHVG10",
    filteredM3Drop10_filteredExpr10 = "filteredExpr10_filteredM3Drop10",
    filteredHVG10_filteredM3Drop10 = "filteredM3Drop10_filteredHVG10")) %>%
  dplyr::ungroup()

res_stab$ARI <- as.numeric(res_stab$ARI)

## ------------------------------------
## Plot stability by k
## ------------------------------------
## load colors
## color set , from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
colors <- c("red", "blue", "green", 
            "orange", "pink", "brown")

filters <- unique(res_stab$filtering)
## manual scale for ggplot
names(colors) <- filters
manual.scale <- scale_colour_manual(name = "", values = colors)

## ARI between filterings at truenclust, by method and dataset
plots[["filterings_truek"]] <- 
  ggplot(res_stab %>% dplyr::group_by(dataset, method, k, filtering) %>%
         dplyr::filter(k == truenclust),
       aes(x = filtering, y = ARI, group = filtering, color = filtering)) + 
  geom_boxplot(na.rm = TRUE) + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset ~ method) +
  labs(y = "ARI between filterings", title = "", x = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

## ARI between filterings by k, by method and dataset
plots[["filterings_byk"]] <- 
  ggplot(res_stab %>% dplyr::group_by(dataset, method, k, filtering) %>%
           dplyr::mutate(ARI = median(ARI)),
       aes(x = k, y = ARI, color = filtering)) + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  geom_line(size = 1) + 
  theme_bw() +
  manual.scale +
  facet_grid(dataset ~ method , scales = "free") +
  ylim(NA, 1) +
  labs(y = "ARI", title = "", x = "Number of clusters")

pdf(gsub("\\.rds$", "_truek.pdf", outrds), width = 20, height = 12)
plots[["filterings_truek"]]
dev.off()

pdf(gsub("\\.rds$", "_byk.pdf", outrds), width = 20, height = 12)
plots[["filterings_byk"]]
dev.off()

saveRDS(plots, file = outrds)
date()
sessionInfo()

