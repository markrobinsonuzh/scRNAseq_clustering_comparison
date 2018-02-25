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
  library(fmsb)
  library(networkD3)
})

## Read clustering results
res <- do.call(rbind, lapply(datasets, function(d) {
  do.call(rbind, lapply(filterings, function(f) {
    do.call(rbind, lapply(methods, function(m) {
      x <- readRDS(paste0("results/sce_", f, "_", d, "_", m, ".rds"))
      dplyr::full_join(x$assignments %>%
                         dplyr::select(dataset, method, cell, run, k, resolution, cluster, trueclass),
                       x$k_estimates %>%
                         dplyr::select(dataset, method, run, k, resolution, est_k)
      ) %>%
        dplyr::full_join(x$timings %>% dplyr::select(dataset, method, run, k, resolution, timing))
    }))
  }))
}))

pdf(gsub("rds$", "pdf", outrds), width = 10, height = 6)

## -------------------------------------------------------------------------- ##
## Calculate performance indices for each method and clustering run
res_summary <- res %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   timing = median(timing)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

print(ggplot(res_summary, aes(x = k, y = ARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset, scales = "free_x") + 
        scale_color_discrete(name = ""))

print(ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
               dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust)),
             aes(x = k, y = medianARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_path() + 
        facet_grid(filtering ~ dataset, scales = "free_x") + 
        scale_color_discrete(name = ""))

## -------------------------------------------------------------------------- ##
## Sankey diagrams
# tmp <- res %>% dplyr::filter(run == 1) %>% 
#   dplyr::filter(dataset == "sce_filtered_Koh") %>% dplyr::filter(method == "CIDR")
# nodes <- unique(paste0("k", tmp$k, "_", tmp$cluster))
# links <- c()
# for (k0 in sort(unique(tmp$k))) {
#   k1 <- k0 + 1
#   if (k1 %in% tmp$k) {
#     for (c0 in unique(tmp$cluster[tmp$k == k0])) {
#       for (c1 in unique(tmp$cluster[tmp$k == k1])) {
#         print(c(k0, k1, c0, c1))
#         l <- length(intersect(tmp$cell[tmp$k == k0 & tmp$cluster == c0],
#                               tmp$cell[tmp$k == k1 & tmp$cluster == c1]))
#         links <- rbind(links, c(match(paste0("k", k0, "_", c0), nodes) - 1, 
#                                 match(paste0("k", k1, "_", c1), nodes) - 1,
#                                 l))
#       }
#     }
#   }
# }
# links <- data.frame(links)
# colnames(links) = c("source", "target", "value")
# links <- subset(links, value != 0)
# nodes <- data.frame(name = nodes, stringsAsFactors = FALSE) %>%
#   dplyr::mutate(group = gsub("k", "", sapply(strsplit(name, "_"), .subset, 1)))
# sankeyNetwork(Links = links, Nodes = nodes,
#               Source = "source", Target = "target",
#               Value = "value", NodeID = "name",
#               fontSize = 0, nodeWidth = 30, NodeGroup = "group", iterations = 1000)

## -------------------------------------------------------------------------- ##
## Plot performance vs timing
print(ggplot(res_summary %>% dplyr::filter(k == truenclust), 
             aes(x = ARI, y = timing, color = method)) +
        geom_point(alpha = 0.2) + 
        geom_point(data = res_summary %>% dplyr::filter(k == truenclust) %>%
                     dplyr::group_by(dataset, filtering, method) %>%
                     dplyr::summarize(ARI = mean(ARI), timing = mean(timing)),
                   size = 3) + 
        facet_grid(filtering ~ dataset, scales = "free"))

## -------------------------------------------------------------------------- ##
## Calculate median performance indices for each method, for the right number of
## cluster and for the estimated number of clusters
print(res_summary %>% dplyr::filter(k == truenclust) %>%
        dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI)) %>%
        ggplot(aes(x = method, y = dataset, fill = medianARI)) + 
        geom_raster() + facet_wrap(~ filtering) + 
        theme(rect = element_blank(), line = element_blank()))

print(res_summary %>% dplyr::filter(k == estnclust) %>%
        dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI)) %>%
        ggplot(aes(x = method, y = dataset, fill = medianARI)) + 
        geom_raster() + facet_wrap(~ filtering) + 
        theme(rect = element_blank(), line = element_blank()))

par(mfrow = c(2, 2), oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1), xpd = NA)
tmp <- res_summary %>% dplyr::filter(k == truenclust) %>%
  dplyr::group_by(dataset, filtering, method, k) %>%
  dplyr::summarize(medianARI = median(ARI))
for (f in unique(tmp$filtering)) {
  tmpp <- tmp %>% dplyr::filter(filtering == f) %>% dplyr::ungroup() %>%
    dplyr::select(method, dataset, medianARI) %>%
    tidyr::spread(dataset, medianARI) %>% as.data.frame()
  rownames(tmpp) <- tmpp$method
  tmpp$method <- NULL
  tmpp <- rbind(rep(1, ncol(tmpp)), rep(0, ncol(tmpp)), tmpp)
  fmsb::radarchart(tmpp, maxmin = TRUE, pcol = 1:8, plty = 1, plwd = 2, 
                   cglcol = "grey", cglty = 1, axislabcol = "grey", 
                   caxislabels = seq(0, 1, length.out = 5), cglwd = 0.8, 
                   axistype = 1, vlcex = 1, centerzero = FALSE, title = f)
  legend(x = 2.5, y = 1.0, legend = rownames(tmpp[-c(1, 2), ]), bty = "n", 
         pch = 20, col = 1:8, text.col = "black", cex = 1, pt.cex = 2, y.intersp = 0.9)
}

print(res_summary %>% dplyr::filter(k == estnclust | k == truenclust) %>% 
        dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust),
                         estnclust = unique(estnclust)) %>%
        dplyr::mutate(n = length(medianARI)) %>%
        dplyr::filter(n == 2) %>% 
        dplyr::summarize(diffARI = medianARI[k == truenclust] - medianARI[k == estnclust]) %>%
        ggplot(aes(x = method, y = dataset, fill = diffARI)) + 
        geom_raster() + facet_wrap(~ filtering) + 
        theme(rect = element_blank(), line = element_blank()))

## -------------------------------------------------------------------------- ##
## Plot the estimated number of clusters
print(ggplot(res_summary %>% dplyr::filter(!is.na(estnclust)), 
             aes(x = method, y = estnclust, color = method)) + 
        geom_hline(aes(yintercept = truenclust), linetype = "dashed") + 
        geom_boxplot(outlier.size = -1) + 
        geom_point(position = position_jitter(width = 0.2, height = 0)) + 
        facet_grid(filtering ~ dataset) + 
        expand_limits(y = 0) + 
        scale_y_continuous(breaks = seq_len(max(res_summary$estnclust, na.rm = TRUE))))

## -------------------------------------------------------------------------- ##
## Calculate stability scores for each method and k, across runs
res <- res %>% dplyr::filter(method != "Seurat") ## Problem since different resolutions can give the same number of clusters, so there may be many "k clusters, run 1" lines. 

tbl <- res %>% dplyr::mutate(ds = paste0(dataset, method, k))
tbl <- split(tbl, f = tbl$ds)
tbl <- do.call(rbind, lapply(tbl, function(w) {
  setNames(data.frame(t(combn(unique(w$run), 2))), c("run1", "run2")) %>%
    dplyr::mutate(dataset = unique(w$dataset), method = unique(w$method), k = unique(w$k)) %>%
    dplyr::select(dataset, method, k, run1, run2)
}))
tbl$ARI <- sapply(seq_len(nrow(tbl)), function(i) {
  tmp <- res %>% dplyr::select(-resolution) %>% 
    dplyr::filter(dataset == tbl$dataset[i] & method == tbl$method[i] & 
                    k == tbl$k[i] & run %in% c(tbl$run1[i], tbl$run2[i])) %>%
    tidyr::spread(run, cluster)
  mclust::adjustedRandIndex(tmp[, as.character(tbl$run1[i])], 
                            tmp[, as.character(tbl$run2[i])])
})
tbl <- tbl %>% tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset"))

print(ggplot(tbl, aes(x = k, y = ARI, group = method, color = method)) + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset, scales = "free_x") + 
        scale_color_discrete(name = ""))

dev.off()

date()
sessionInfo()

