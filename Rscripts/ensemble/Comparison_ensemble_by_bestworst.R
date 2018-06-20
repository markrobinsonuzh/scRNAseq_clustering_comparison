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

sum_summary <- res_summary %>% dplyr::group_by(dataset, method, run,k) %>%
               dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass), truenclust = length(unique(trueclass))) 
# format
sum_ensmbl2$run <- as.integer( sum_ensmbl2$run )

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
sum_all3 <-  left_join(sum_all2 , sum_summary, by=c( "dataset" ,"methtwo"="method", "run", "k"))
colnames(sum_all3 ) <- c(
  "dataset",
  "method.ensmbl",
  "methone",
  "methtwo",
  "run",
  "k",
  "ARIensmbl",
  "truenclustz",
  "ARIone",
  "truenclust.y",
  "ARItwo")
sum_all3 <- select(sum_all3,1:11)

#rename, compute difference 1 based on the best method (best method defined as method with max in median per run per k) per dataset
# median ARI of method one

sum_max <-  sum_all3 %>% select( dataset, methone, methtwo ,k, ARIensmbl,run,ARIone,ARItwo,truenclust.y   )%>%
  group_by( dataset, methone,methtwo,  k) %>% 
  summarise( med.ARIone= median(ARIone, na.rm=TRUE),
             med.ARItwo=median(ARItwo, na.rm=TRUE),
             med.ensemble=median(ARIensmbl, na.rm=TRUE),
             truenclust=unique(truenclust.y))%>% ungroup()

# best method per ensemble

# diff
sum_max2 <- sum_max %>% dplyr::mutate(best = pmax(med.ARIone, med.ARItwo, na.rm=TRUE),
                                      worst= pmin(med.ARIone, med.ARItwo, na.rm=TRUE) ,
                                      diff.best=med.ensemble-best,
                                      diff.worst=med.ensemble-worst)
df.long <- reshape2::melt(sum_max2 , id=c( "dataset","methone","methtwo" ,"k","truenclust", "med.ARIone" ,"med.ARItwo", "med.ensemble","best" ,"worst"))
df.long2 <- df.long %>% mutate(ARI=case_when(variable=="diff.best"~best, variable=="diff.worst"~worst))
# best method per k and dataset 

# ------------------------------------
# Plot
# ------------------------------------
pdf("plots/ensemble/plot_ensembles_differences_bestworst.pdf", width = 10, height = 7)

# plot histo of differences in ARI, by the two methods 
ggplot( df.long2%>%filter() , aes(variable, value))+
  geom_point(position="jitter", alpha=0.2, aes(color=ARI) )+
  geom_boxplot(alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept=0, linetype="dashed")+
  viridis::scale_color_viridis(option="C", direction=-1)+
  labs(color="ARI best/worst",x="difference to ensemble", y="ARI difference", 
       title="difference between ensemble and best/worst method, for all k")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title = element_text(size=15) )+
  scale_x_discrete(labels=c("diff.best" = "ensemble-best", "diff.worst" = "ensemble-worst"))

ggplot( df.long2%>%filter(k==truenclust) , aes(variable, value))+
  geom_point(position="jitter", alpha=0.2, aes(color=ARI) )+
  geom_boxplot(alpha=0.5)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+
  viridis::scale_color_viridis(option="C", direction=-1)+
  labs(color="ARI best/worst",x="difference to ensemble", y="ARI difference", 
       title="difference between ensemble and best/worst method, for truenclust")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title = element_text(size=15) )+
  scale_x_discrete(labels=c("diff.best" = "ensemble-best", "diff.worst" = "ensemble-worst"))

dev.off()
