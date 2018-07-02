args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(datadir)
print(datasets)
print(outrds)

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
})

datasets_full <- paste0(datadir, "/sce_full/sce_full_", datasets, ".rds")
names(datasets_full) <- gsub("\\.rds", "", basename(datasets_full))

full_data <- lapply(datasets_full, function(x) {
  readRDS(x)
})

plots <- lapply(names(full_data), function(nm) {
  x <- full_data[[nm]]
  df_to_plot <- as.data.frame(SingleCellExperiment::reducedDim(x, "TSNE")) %>%
    tibble::rownames_to_column("cell") %>% 
    dplyr::left_join(as.data.frame(colData(x)) %>% tibble::rownames_to_column("cell"), 
                     by = "cell") %>%
    dplyr::mutate(dataset = nm)
  
  ggplot(df_to_plot, 
         aes(x = V1, y = V2)) + 
    xlab("dimension 1") + 
    ylab("dimension 2") + 
    geom_rug(colour = "gray20", 
             alpha = 0.65) + 
    theme_bw() +
    geom_point(aes(colour = phenoid), 
               alpha = 0.65) +
    scale_color_discrete(name = "") + 
    labs(title = gsub("sce_full_", "", nm))
})

pdf(gsub("rds$", "pdf", outrds), width = 16, height = 16)
cowplot::plot_grid(plotlist = plots, labels = "", ncol = 3)
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()