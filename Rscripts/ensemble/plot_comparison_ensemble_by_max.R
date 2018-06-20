# ------------------------------------------------------------------------
# Comparing ensemble clusterings, 
# against best performing method per dataset and k
# ------------------------------------------------------------------------

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
  library(ggtree)
  
})

## Read clustering and ensemble results
res_ensmbl3 <- readRDS(file="output/ensemble/clustering_ensemble_allmethods3_by_k.rds")
res_ensmbl2 <- readRDS(file="output/ensemble/clustering_ensemble_allmethods2_by_k.rds")
res_summary <- readRDS(file="output/clustering_summary/clustering_summary.rds")

# ------------------------------------
# compute ARI, no of unique clusters
# ------------------------------------
sum_ensmbl2 <- res_ensmbl2 %>% dplyr::group_by(dataset, method, run,k) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo"), remove=FALSE) %>% 
  dplyr::ungroup()

sum_ensmbl3 <- res_ensmbl3 %>% dplyr::group_by(dataset, method, run,k) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo", "methtree"), remove=FALSE) %>% 
  dplyr::ungroup()

sum_summary <- res_summary %>% dplyr::group_by(dataset, method, run,k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass))) 
# format
sum_ensmbl2$run <- as.integer( sum_ensmbl2$run )
sum_ensmbl3$run <- as.integer( sum_ensmbl3$run )
sum_summary$k <- as.factor(sum_summary$k )
# join datasets
sum_all2 <-  left_join(sum_ensmbl2 , sum_summary, by=c( "dataset" ,"methone"="method", "run", "k"))
colnames(sum_all2) <- c(
  "dataset",
  "method.ensmbl",
  "methone",
  "methtwo",
  "run",
  "k",
  "ARIensmbl",
  "truenclustz",
  "ARIone",
  "truenclust.y")
sum_all3 <- left_join(sum_all2 , sum_ensmbl3, by=c( "dataset" ,"methone"="methone","methtwo"="methtwo", "run", "k"))
colnames(sum_all3) <- c(
  "dataset",
  "method.ensmbl2",
  "methone",
  "methtwo",
  "run",
  "k",
  "ARIensmbl2",
  "truenclustx",
  "ARIone",
  "truenclust",
  "method.ensembl3",
  "methtree",
  "ARIensmbl3",
  "truenclustyy")

#rename, compute difference 1 based on the best method (best method defined as method with max in median per run per k) per dataset
# median ARI of method one

sum_max <-  sum_all3 %>% select(dataset,methone,methtwo , method.ensmbl2, k, run,ARIone,ARIensmbl2, truenclust   )%>%
  group_by( dataset, methone, k) %>% 
  mutate( med.ARIone= median(ARIone, na.rm=TRUE) )%>% ungroup()
# best method per k and dataset 

sum_max <-  sum_max  %>% group_by( dataset, k) %>% 
  mutate( max.ARI= max(med.ARIone, na.rm=TRUE) )%>% ungroup()
#difference between Ensembl and best method
res_median <-  sum_max %>% group_by( dataset, methone,k) %>% 
  mutate(ARIdiff1=ARIensmbl2-max.ARI) %>% distinct()%>%
  ungroup()
# ------------------------------------
# PLot
# ------------------------------------
pdf("plots/ensemble/plot_comparison_ensemble_max.pdf", width=15,height=8)

# plot histo of differences in ARI, by the two methods 
ggplot( res_median , aes(ARIdiff1))+
  geom_histogram() +
  facet_grid(methone~ methtwo)+
  geom_vline(xintercept = 0, linetype = "dashed", colour=2) +
  labs(title="Difference ensemble to best performing method, all k")

# plot histo of differences in ARI, by filtering

ggplot( res_median %>%
          tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
          dplyr::select(-sce) %>% dplyr::ungroup(), 
        aes(ARIdiff1, fill=dataset) )+
  geom_histogram(aes(ARIdiff1, fill=dataset)) +
  scale_fill_discrete()+
  facet_grid(dataset~filtering)+
  geom_vline(xintercept = 0, linetype = "dashed", colour=1)+
  labs(title="Difference ensemble to best performing method, by dataset")
dev.off()

# plot histo of differences in ARI, seperated by datasets

ggplot( res_median %>%filter(k==truenclust)%>% dplyr::ungroup()
        , aes(ARIdiff1, fill=dataset))+
  geom_dotplot() +
  facet_grid(methone~ methtwo)+
  geom_vline(xintercept = 0, linetype = "dashed", colour=2) +
  labs(title="Difference ensemble to best performing method, by dataset")

