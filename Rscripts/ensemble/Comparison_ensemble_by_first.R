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
## Read clustering results
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
#rename, compute difference
sum_all <- sum_all3 %>% mutate(ARIdiff1=ARIensmbl2-ARIone, ARIdiff2=ARIensmbl3-ARIone)


# compute median
res_median <- sum_all%>% group_by(dataset, method.ensmbl2,method.ensembl3, k, truenclust, methone, methtwo, methtree ) %>% 
  dplyr::summarize(ARIensmbl2=median(ARIensmbl2), ARIensmbl3=median(ARIensmbl3),ARIone=median(ARIone), ARIdiff1=median(ARIdiff1), ARIdiff2=median(ARIdiff2) )


pdf("plots/ensemble/plot_diff_ensemble_allk.pdf", width = 12, height = 12)
ggplot( res_median , aes(ARIdiff1))+
  geom_histogram() +
  facet_grid(methone~ methtwo, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed", colour=2) +
  theme(axis.text.x = element_text(size=6)) +
  labs(title="Difference single method vs ensemble, all k")
dev.off()
pdf("plots/ensemble/plot_diff_ensemble_truenclust.pdf", width = 12, height = 12)
ggplot( res_median%>%filter(k==truenclust) , aes(ARIdiff1))+
  geom_histogram() +
  facet_grid(methone~ methtwo, scales = "free")+
  geom_vline(xintercept = 0, linetype = "dashed", colour=2) +
  theme(axis.text.x = element_text(size=6)) +
  labs(title="Difference single method vs ensemble, truenclust")
dev.off()

