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
# load data
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

# load summaries
res <- readRDS(file = clusteringsummary)
res <-  res %>%group_by(dataset, method,run) %>%
  mutate( cluster=as.numeric(cluster), truenclust=length(unique(trueclass) ))

# remove Seurat, k == truenclust
# what to do with the different runs(test all, conensus)? choose only one as starting point
res_filt <- res %>% filter( k==truenclust, run==1, !method %in% c("Seurat"))
# combine dataframes 
## select what
select <- as.name(c)
sel <- c("total_counts", "total_features", "pct_counts_top_50_features")
as.name(sel[2] )

tbl <- do.call(rbind, lapply( names(all_data), function(nm) {
  x <- all_data[[nm]]
  colData <- colData(x) %>%
    as.data.frame() %>% 
    select( total_counts, total_features, pct_counts_top_50_features) %>%
    mutate(cell=row.names(.), dataset=nm )
}))
# nest dataframe, remove NA groups
comb.tbl <- left_join(res_filt,tbl, by=c("cell", "dataset")) %>%
  select(-resolution) %>%
  group_by(dataset,method, k, run) %>%
  filter(any(!is.na(cluster))) %>%
  nest()

## test on total_features, total_counts
# functions for test with KW
kw.fun <- function(df, formula){
  kruskal.test(formula=formula, data=df)
}
# extract pval
p_kw <- function(mod){
  mod$p.value
}

res.pval <- comb.tbl %>% mutate( mod.counts= map(data, kw.fun, formula=as.formula(total_counts~cluster)),
                          mod.features=map(data, kw.fun, formula=as.formula(total_features~cluster)) ) %>% 
  transmute(dataset,method, p.counts = map_dbl(mod.counts, p_kw), 
                            p.features =  map_dbl(mod.features, p_kw))  %>%
  mutate(p.counts=p.adjust(p.counts, "BH"),
         p.features=p.adjust(p.features, "BH"))

# plot pvals
# format
res.pval <- res.pval%>%tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce)
# plot
ggplot2::ggplot(res.pval,
                aes(x=method,y=p.features)) +
  geom_point(aes(color=dataset),position="jitter") +
  #facet_grid(filtering~dataset, scales = "free" ) +
  geom_hline(yintercept=0.01 ) +
  #facet_grid(.~dataset) +
  theme(axis.text.x = element_text(angle=90))

ggplot2::ggplot(res.pval,
                aes(x=method,y=p.counts)) +
  geom_point(aes(color=dataset),position="jitter") +
  #facet_grid(filtering~dataset, scales = "free" ) +
  geom_hline(yintercept=0.01 ) +
  #facet_grid(.~dataset) +
  theme(axis.text.x = element_text(angle=90))

# To do , add mt genes etc. as covariates
lapply(all_data,function(x){
  sum(grepl("mt", rownames(x)), na.rm = TRUE) 
})
lapply(all_data,function(x){
  sum(grepl("rpl", rownames(x)), na.rm = TRUE) 
})





