args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

filterings <- strsplit(filterings, ",")[[1]]
datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datadir)
print(filterings)
print(pctkeep)
print(datasets)
print(clusteringsummary)
print(outrds)

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
})

res_summary <- readRDS(file = clusteringsummary)

datasets_filtered <- c(outer(paste0(datadir, "/sce_filtered", filterings, pctkeep, "/sce_filtered", 
                                    filterings, pctkeep), 
                             paste0(datasets, ".rds"),
                             FUN = paste, sep = "_"))
names(datasets_filtered) <- gsub("\\.rds", "", basename(datasets_filtered))

full_data <- lapply(datasets_filtered, function(x) {
  readRDS(x)
})

datasets <- c("sce_filteredExpr10_Zhengmix4eq", "sce_filteredExpr10_Zhengmix4eq",
              "sce_filteredExpr10_Zhengmix4uneq", "sce_filteredExpr10_Trapnell")
nrun <- rep(1, 4)
methods <- c("FlowSOM", "CIDR", "TSCAN", "SC3")

plot_tSNE <- function(res_summary, full_data, meth, dat, nrun){
  res <- res_summary %>% 
    dplyr::filter(dataset %in% dat, 
                  method %in% meth, run == nrun) %>%
    group_by(dataset) %>%
    filter(k == length(unique(trueclass))) %>%
    select(cell, dataset, cluster, trueclass)
  sel.data <- full_data[[dat]]
  df_to_plot <- as.data.frame(SingleCellExperiment::reducedDim(sel.data, "TSNE")) %>%
    tibble::rownames_to_column("cell") %>% 
    dplyr::left_join(res, by = "cell")
  
  ggplot(df_to_plot, 
         aes(x = V1, y = V2)) + 
    xlab("dimension 1") + 
    ylab("dimension 2") + 
    geom_rug(colour = "gray20", 
             alpha = 0.65) + 
    theme_bw() +
    geom_point(aes(colour = cluster, 
                   shape = trueclass), 
               alpha = 0.65)+
    guides(shape = guide_legend(title = "True class"),
           colour = guide_legend(title = "Clusters")) +
    labs(title = paste(dat, meth))
}

plots <- lapply(c(1:4), function(i) {
  plot_tSNE(res_summary, full_data, methods[i], datasets[i], nrun[i])
})

pdf(gsub("rds$", "pdf", outrds), width = 20, height = 15)
cowplot::plot_grid(plotlist = plots, labels = "auto")
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()