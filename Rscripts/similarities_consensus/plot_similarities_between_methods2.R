args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(consensusrds)
print(methodchart)
print(outrds)

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(clusterExperiment)
  library(clue)
  library(multidplyr)
  library(ggplot2)
  library(viridis)
  library(ggthemes)
  library(pheatmap)
  library(reshape2)
  library(mclust)
  library(RColorBrewer)
  library(ggtree)
  library(purrr)
  library(cowplot)
  library(ape)
})

# Load files
res <- readRDS(file = consensusrds)

### Compare clustering between methods by consensus ARI, k = truenclust
#------------------------------------------------------------------
plot_crossmethod_concordance <- function(res, ncluster) {
  # group data
  set.seed(42)
  res.cons <- res %>% dplyr::group_by(dataset, method, k) %>% 
    dplyr::mutate(truenclust = length(unique(trueclass))) %>% 
    dplyr::filter(k == ifelse(ncluster == 0, truenclust, ncluster), !method %in% c("Seurat")) %>% 
    dplyr::select(dataset, method, consensus, k, cell)
  # for Seurat, pick resolution param 
  res.cons.seurat <- res %>% dplyr::group_by(dataset, method, k) %>% 
    dplyr::filter(method %in% c("Seurat")) %>% 
    dplyr::mutate(truenclust = length(unique(trueclass))) %>% 
    dplyr::filter(k == ifelse(ncluster == 0, truenclust, ncluster)) %>%
    filter(resolution == sample(resolution, 1)) %>% 
    dplyr::select(dataset, method, consensus, k, cell)
  res.cons <- dplyr::bind_rows(res.cons , res.cons.seurat)
  
  # compute  matrix with pairwise ARIs
  l <- vector("list", length = length(unique(res.cons$dataset)))
  names(l) <- unique(res.cons$dataset)
  
  for (i in unique(res.cons$dataset)) {
    print(i)
    x <- subset(res.cons, dataset == i)
    m <- matrix(NA, nrow = length(unique(x$method)), ncol = length(unique(x$method)))
    rownames(m) <- unique(x$method)
    colnames(m) <- unique(x$method)
    for (j in unique(x$method)) {
      print(j)
      c1 <- subset(x, method == j) 
      for (u in unique(x$method)) {
        print(u)
        c2 <- subset(x, method == u) 
        ## get shared cells, calculate ARI
        scl <- intersect(c1$cell, c2$cell)
        m[j, u] <- mclust::adjustedRandIndex(c1$consensus[match(scl, c1$cell)], 
                                             c2$consensus[match(scl, c2$cell)])
      }
    }
    l[[i]] <- m
  }
  
  # average ARI over datasets
  df <- reshape2::melt(l, value.name = "ARI")
  
  df.median <- df %>% dplyr::group_by(Var1, Var2) %>% 
    dplyr::summarise(median = median(ARI, na.rm = TRUE)) %>%
    dplyr::arrange(Var1, Var2) %>% dplyr::ungroup()
  
  ## Calculate average area under concordance curve across all data set
  ## instances, for each pair of methods
  
  m.median <- reshape2::acast(df.median, Var1 ~ Var2, value.var = "median") 
  stopifnot(all(rownames(m.median) == colnames(m.median)))
  
  # plot all ARIs, facet by filter and dataset
  if (ncluster == 0) {
    ttl <- ""
  } else {
    ttl <- paste0("Number of clusters = ", ncluster)
  }
  heatmap.dataset <- 
    ggplot(df %>% tidyr::separate(L1, sep = "_", into = c("sce", "filtering", "dataset")) %>%
             dplyr::select(-sce),
           aes(x = Var1, y = Var2, fill = ARI)) +
    geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
    facet_grid(filtering ~ dataset, drop = FALSE) +
    scale_fill_viridis(name = "ARI", direction = -1, na.value = "white") +
    theme_tufte(base_family = "Helvetica") +
    labs(x = NULL, y = NULL, title = ttl) +
    coord_equal() +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.title.align = 0,
          legend.text = element_text(size = 6),
          legend.position = "right",
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.2, "cm"),
          axis.ticks = element_blank())
  
  # plot median ARIs from datasets
  heatmap.median <- ggplot(df.median, aes(x = Var1, y = Var2, fill = median)) +
    geom_tile(color = "white", size = 0.5, na.rm = FALSE) +
    scale_fill_viridis(name = "ARI", direction = -1) +
    theme_tufte(base_family = "Helvetica") +
    labs(x = NULL, y = NULL, title = paste0("Median similarity ", ttl)) +
    coord_equal() +
    theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.title.align = 0,
          legend.text = element_text(size = 6),
          legend.position = "right",
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.2, "cm"),
          axis.ticks = element_blank())
  
  #-------------------------------------------------------------------
  
  # print tree
  ## Get all subclusters from an hclust object, from  https://github.com/csoneson/conquer_comparison/blob/master/scripts/help_function_crossmethod_concordance.R
  get_subclusters <- function(hcl) {
    m <- hcl$merge
    labs <- hcl$labels
    L <- list()
    for (i in seq_len(nrow(m))) {
      tmp <- c()
      if (m[i, 1] < 0) tmp <- c(tmp, labs[-m[i, 1]])
      else tmp <- c(tmp, L[[m[i, 1]]])
      if (m[i, 2] < 0) tmp <- c(tmp, labs[-m[i, 2]])
      else tmp <- c(tmp, L[[m[i, 2]]])
      L[[i]] <- sort(tmp)
    }
    L
  }
  
  ## Hierarchical clustering based on 1 - m.median, 
  hcl_average <- hclust(as.dist(1 - m.median)) 
  
  ## Get all subclusters
  subclusters_average <- get_subclusters(hcl_average)
  
  ## Get all subclusters for individual data set instances
  cmtmp <- df 
  uniqvals <- unique(cmtmp$L1) ## data sets
  # per dataset i
  subclusters_all <- lapply(uniqvals, function(i) {
    print(i)
    cmtmp2 <- cmtmp %>% dplyr::filter(L1 == i) %>% dplyr::select(Var1, Var2, ARI) # for each dataset
    cmtmp2 <-  reshape2::acast(cmtmp2, Var1 ~ Var2, value.var = "ARI")
    stopifnot(all(rownames(cmtmp2) == colnames(cmtmp2)))
    stopifnot(all((cmtmp2 == t(cmtmp2))[!is.na(cmtmp2 == t(cmtmp2))]))
    cmtmp2 <- cmtmp2[apply(cmtmp2, 1, function(x){!all(is.na(x))}), ]
    cmtmp2 <- cmtmp2[, apply(cmtmp2, 2, function(x){!all(is.na(x))})]
    get_subclusters(hclust(as.dist(1 - cmtmp2)))
  })
  
  ## Get stability values for each subcluster in subclusters_average
  ## find datasets with a missing method element in subclusters_all 
  which_missing <- function(x, invert, toMatch) {
    toMatch <- c(paste0(toMatch))
    m <- purrr::map(x, function(x) {
      grep(paste(toMatch, collapse = "|"), x = x, value = FALSE, 
           invert = FALSE, ignore.case = TRUE) %>%
        is_empty
    })
    m <- which(m == invert)
    return(m)
  }
  
  ## for this particular datasets (column in w), add the missing method in the
  ## particular branch, such that the we have a TRUE in w (if the rest is the
  ## same)
  x <- sapply(
    subclusters_all, function(w) {
      subclusters_average %in% w 
    })
  
  p <- rownames(m.median)
  # which missing
  for (i in p) {
    x[which_missing(subclusters_average, invert = FALSE, toMatch = paste(i)), 
      which_missing(subclusters_all, invert = TRUE, toMatch = paste(i))] <- TRUE
  }
  stability_scores <- rowMeans(x)
  
  subcl <- subclusters_average
  stab <- stability_scores
  
  ggt <- ggtree(as.phylo(hcl_average))
  for (m in seq_len(length(subcl))) {
    tryCatch({
      i <- MRCA(ggt, subcl[[m]])
      if (stab[m] >= 0)
        ggt$data[i, "label"] <- round(stab[m], 2)
    }, error = function(e) NULL)
  }
  
  trees <- ggt + geom_label2(aes(subset = !isTip, label = label), size = 5) + 
    geom_tiplab(aes(angle = 90), hjust = 1, size = 10) + 
    ggplot2::scale_x_reverse() + ggplot2::coord_flip() + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "mm")) + 
    xlim_tree(0.75) +
    ggtitle(ttl)
  
  # add annotation to tree
  annot <- readr::read_delim(methodchart, ";", escape_double = FALSE, trim_ws = TRUE)
  tiporder <- trees$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
  annot <- annot[match(tiporder$label, annot$method), 
                 c("dimension", "clustering", "counts"), drop = FALSE]
  
  # colors for heatmap
  colors <- list(
    dimension = c(PCA = "#C10359", tSNE = "#fd9ec9", Various = "#FF33FF", None = "#fed8e9"),
    clustering = c(Hierarchical = "#01368C", Graph = "#2F6CCE", Kmeans = "#93B8F2",
                   SOM = "#d1e1fa", ModelBased = "#030a17", Various = "#9933FF"),
    counts = c(Raw = "#2EA801", LogNorm = "#96E878", Various = "#DBFC07")
  )
  
  g1 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = dimension)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$dimension) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 1, 0, 1), "mm"))
  g2 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = clustering)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$clustering) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 1, 0, 1), "mm"))
  g3 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = counts)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$counts) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 1, 0, 1), "mm"))
  
  
  # legends
  gds <- plot_grid(get_legend(g1 + 
                                theme(legend.position = "bottom") + 
                                guides(fill = 
                                         guide_legend(ncol = 4,
                                                      title = "Dimension reduction", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 20,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 16,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   get_legend(g2 + 
                                theme(legend.position = "left") + 
                                guides(fill = 
                                         guide_legend(ncol = 4,
                                                      title = "Clustering method", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 20,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 16,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   get_legend(g3 + 
                                theme(legend.position = "left") + 
                                guides(fill = 
                                         guide_legend(ncol = 4,
                                                      title = "Input", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 20,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 16,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   
                   nrow= 1)
  
  ann.tree <- plot_grid(trees + theme(plot.margin = unit(c(0, 0, -150, 0), "mm")), 
                        g1 + theme(legend.position = "none"),
                        g2 + theme(legend.position = "none"),
                        g3 + theme(legend.position = "none"),
                        gds + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")), 
                        rel_heights = c(3, 0.3, 0.3, 0.3, 0.3), ncol = 1)
  
  # store plots in plot.list 
  plot.list <- list()
  plot.list[["heatmap_dataset"]] <- heatmap.dataset
  plot.list[["heatmap_median"]] <- heatmap.median
  plot.list[["ann_tree"]] <- ann.tree
  return(plot.list)
}

# all
list.k <- as.list(c(0, 3:10))

# plot 
plots <- lapply(list.k, function(x) {plot_crossmethod_concordance(res, ncluster = x)})

pdf(gsub("\\.rds$", "_median_all.pdf", outrds), width = 15, height = 10)
cowplot::plot_grid(plotlist = sapply(plots , "[", "heatmap_median"))
dev.off()

pdf(gsub("\\.rds$", "_median_all_truenclust.pdf", outrds), width = 7, height = 7)
cowplot::plot_grid(plotlist = sapply(plots[1] , "[", "heatmap_median"))
dev.off()

pdf(gsub("\\.rds$", "_dataset_all.pdf", outrds), width = 50, height = 50)
cowplot::plot_grid(plotlist = sapply(plots , "[", "heatmap_dataset"), ncol = 2)
dev.off()

pdf(gsub("\\.rds$", "_dataset_all_truenclust.pdf", outrds), width = 15, height = 10)
cowplot::plot_grid(plotlist = sapply(plots[1] , "[", "heatmap_dataset"))
dev.off()

pdf(gsub("\\.rds$", "_tree_all.pdf", outrds), width = 15, height = 10)
cowplot::plot_grid(plotlist = sapply(plots , "[", "ann_tree"))
dev.off()

pdf(gsub("\\.rds$", "_tree_all_truenclust.pdf", outrds), width = 15, height = 10)
cowplot::plot_grid(plotlist = sapply(plots[1], "[", "ann_tree"))
dev.off()

saveRDS(plots, file = outrds)
date()
sessionInfo()

