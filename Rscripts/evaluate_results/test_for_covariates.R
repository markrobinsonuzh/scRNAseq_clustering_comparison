args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

filterings <- strsplit(filterings, ",")[[1]]
datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(filterings)
print(pctkeep)
print(datasets)
print(datadir)
print(clusteringsummary)
print(outrds)

#------------------------------------------------
# Test covariates
#------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(SingleCellExperiment)
  library(ggplot2)
})

## load colors
source("Rscripts/Colorscheme.R")

## load data
datasets_filtered <- c(outer(paste0(datadir, "/sce_filtered", filterings, pctkeep, "/sce_filtered", 
                                    filterings, pctkeep), 
                             paste0(datasets, ".rds"),
                             FUN = paste, sep = "_"))
names(datasets_filtered) <- gsub("\\.rds", "", basename(datasets_filtered))

datasets_full <- paste0(datadir, "/sce_full/sce_full_", datasets, ".rds")
names(datasets_full) <- gsub("\\.rds", "", basename(datasets_full))

all_data <- lapply(c(datasets_filtered, datasets_full), function(x) {
  readRDS(x)
})

## load summaries
res <- readRDS(file = clusteringsummary)
res <- res %>% dplyr::group_by(dataset, method, run) %>%
  dplyr::mutate(cluster = as.numeric(cluster), truenclust = length(unique(trueclass)))

## remove Seurat, k == truenclust
## Select one run
res_filt <- res %>% dplyr::filter(k == truenclust, run == 1, !method %in% c("Seurat"))

## combine dataframes 
## select what
sel <- c("total_counts", "total_features", "pct_counts_top_50_features")

tbl <- do.call(rbind, lapply(names(all_data), function(nm) {
  x <- all_data[[nm]]
  colData <- colData(x) %>%
    as.data.frame() %>% 
    dplyr::select(total_counts, total_features, pct_counts_top_50_features) %>%
    dplyr::mutate(cell = row.names(.), dataset = nm)
}))

## nest dataframe, remove NA groups
comb.tbl <- dplyr::inner_join(res_filt, tbl, by = c("cell", "dataset")) %>%
  dplyr::select(-resolution) %>%
  dplyr::group_by(dataset,method, k, run) %>%
  dplyr::filter(any(!is.na(cluster))) %>%
  nest()

## test on total_features, total_counts
## functions for test with KW
kw.fun <- function(df, formula){
  kruskal.test(formula = formula, data = df)
}

## extract pval
p_kw <- function(mod) {
  mod$p.value
}

res.pval <- comb.tbl %>% 
  dplyr::mutate(mod.counts = map(data, kw.fun, formula = as.formula(total_counts~cluster)),
                mod.features = map(data, kw.fun, formula = as.formula(total_features~cluster))) %>% 
  transmute(dataset,method, p.counts = map_dbl(mod.counts, p_kw), 
            p.features =  map_dbl(mod.features, p_kw)) %>%
  mutate(p.counts = p.adjust(p.counts, "BH"),
         p.features = p.adjust(p.features, "BH"))

## plot pvals
## format
res.pval <- res.pval %>% 
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce)

## plot
pdf(gsub("rds$", "pdf", outrds), width = 12, height = 8)
ggplot2::ggplot(res.pval,
                aes(x = method, y = -log10(p.features))) +
  geom_point(aes(color = method, shape = filtering), size = 3, 
             position = position_jitter(width = 0.1)) +
  expand_limits(y = 0) + manual.scale + 
  theme_bw() + facet_wrap(~ dataset, scales = "free_y") + 
  labs(y = "-log10(adjusted p-value)", x = "", 
       title = "Association between number of detected features and cluster assignment") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 16))

ggplot2::ggplot(res.pval,
                aes(x = method, y = -log10(p.counts))) +
  geom_point(aes(color = method, shape = filtering), size = 3, 
             position = position_jitter(width = 0.1)) +
  expand_limits(y = 0) + manual.scale + 
  theme_bw() + facet_wrap(~ dataset, scales = "free_y") + 
  labs(y = "-log10(adjusted p-value)", x = "",
       title = "Association between number of counts and cluster assignment") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 16))
dev.off()

saveRDS(NULL, outrds)
date()
sessionInfo()





