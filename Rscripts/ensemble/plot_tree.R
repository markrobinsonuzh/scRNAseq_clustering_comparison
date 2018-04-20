
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(ggthemes)
  library(viridis) 
  library(data.tree)
  
})
## Read clustering results
res_ensmbl <- readRDS(file="output/ensemble/clustering_ensemble_allmethods2.rds")
res_summary <- readRDS(file="output/clustering_summary/clustering_summary.rds")


# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
sum_ensmbl <- res_ensmbl %>% dplyr::group_by(dataset, method, run) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce)  %>%
  tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo", "methtree"), remove=FALSE) %>% 
  dplyr::ungroup()

sum_summary <- res_summary %>% dplyr::group_by(dataset, method, run,k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
# join datasets
huhu <- inner_join(sum_ensmbl , sum_summary, by=c("filtering", "dataset", "methone"="method" ))
huhu2 <- left_join(sum_ensmbl , sum_summary, by=c("filtering", "dataset" ))
huhu
res_median <- res_summary%>% group_by(dataset, filtering, method) %>% dplyr::summarize(medianARI=median(ARI))%>%
tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo", "methtree"), remove=FALSE)

res_median.sub <- filter(res_median, dataset=="Kumar")

# work with one dataset

res_median$pathString <- paste("dataset", 
                               res_median$dataset,
                               res_median$filtering,
                               res_median$methone,res_median$methtwo,
                               sep="/")

ARI <- as.Node(res_median)
print(ARI, "medianARI", limit=20)
plot(ARI)
