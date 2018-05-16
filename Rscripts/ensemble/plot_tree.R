
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
res_ensmbl3 <- readRDS(file="output/ensemble/clustering_ensemble3_truenclust.rds")
res_ensmbl2 <- readRDS(file="output/ensemble/clustering_ensemble2_truenclust.rds")
res_summary <- readRDS(file="output/clustering_summary/clustering_summary.rds")

# ------------------------------------
# compute ARI, no of unique clusters
# ------------------------------------

sum_ensmbl2 <- res_ensmbl2 %>% dplyr::group_by(dataset, method, run) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo"), remove=FALSE) %>% 
  dplyr::ungroup()

sum_ensmbl3 <- res_ensmbl3 %>% dplyr::group_by(dataset, method, run) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(method, sep="[[:punct:]]", into=c("methone","methtwo", "methtree"), remove=FALSE) %>% 
  dplyr::ungroup()

sum_summary <- res_summary %>% dplyr::group_by(dataset, method, run,k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass))) 
sum_ensmbl2$run <- as.integer( sum_ensmbl2$run )
sum_ensmbl3$run <- as.integer( sum_ensmbl3$run )

# join datasets
sum_all2 <- left_join(sum_ensmbl2 , sum_summary, by=c( "dataset" ,"methone"="method", "run"))
colnames(sum_all2) <- c(
  "dataset",
  "method.ensmbl",
  "methone",
  "methtwo",
  "run",
  "ARIensmbl",
  "truenclustz",
  "k",
  "ARIone",
  "truenclust.y")
sum_all3 <- left_join(sum_all2 , sum_ensmbl3, by=c( "dataset" ,"methone"="methone","methtwo"="methtwo", "run"))
colnames(sum_all3) <- c(
  "dataset",
  "method.ensmbl2",
  "methone",
  "methtwo",
  "run",
  "ARIensmbl2",
  "truenclustx",
  "k",
  "ARIone",
  "truenclust",
  "method.ensembl3",
  "methtree",
  "ARIensmbl3",
  "truenclustyy")
#rename, compute difference
sum_all <- sum_all3 %>% mutate(ARIdiff1=ARIensmbl2-ARIone, ARIdiff2=ARIensmbl3-ARIensmbl2)

# compute median
res_median <- sum_all%>% group_by(dataset, method.ensmbl2,method.ensembl3, k, truenclust, methone, methtwo, methtree ) %>% 
  dplyr::summarize(ARIensmbl2=median(ARIensmbl2), ARIensmbl3=median(ARIensmbl3),ARIone=median(ARIone), ARIdiff1=median(ARIdiff1), ARIdiff2=median(ARIdiff2) )
list <- vector("list",length( unique((res_median$dataset)) ))
names(list) <-  unique(names(res_median$dataset) )

for(i in unique(res_median$dataset) ){ 
  print(i)
  res_median.sub <- filter(res_median, k==truenclust, dataset==i)
  
  # work with one dataset
  res_median.sub$pathString <- paste("dataset", 
                                     res_median.sub$dataset,
                                     res_median.sub$methone,
                                     res_median.sub$method.ensmbl2,
                                     res_median.sub$method.ensembl3,
                                     sep="/")
  # create tree using data.tree
  ARItree <- as.Node(res_median.sub)
  annotation <- ToDataFrameTypeCol(ARItree,"ARIone", "ARIdiff1","ARIdiff2")
  # add layers
  anno1 <- select(annotation, level_3, ARIone)
  anno2 <- select(annotation, level_4, ARIdiff1)
  anno3 <- select(annotation, level_5, ARIdiff2)
  # format as phylo object
  gh <- as.phylo.Node(ARItree)
  p <- ggtree(gh, layout = 'rectangular')
  
  p <- p %<+% anno1
  p <- p %<+% anno2
  p <- p %<+% anno3
  # plot tree
  # change names for 2nd layer and 3rd layer
  # p$data$label <- gsub('.*\\.',"", p$data$label)
  p2 <- p+ geom_text(aes(color=ARIone , label=round( ARIone ,2) ), size=8, hjust=1, na.rm = TRUE)+ # median ARi
    geom_text(size=8,hjust=-1,aes(label=round(ARIdiff2, 1)))+ # tip difference
    geom_text(size=8,hjust=2,aes(label=round(ARIdiff1, 1)))+ # node difference
    geom_nodelab(geom="label",aes(x=branch), size=8)+ #  node labels
    geom_tiplab(geom="label",aes(x=branch), size=8) # tip labels
    
    list[[i]] <- p2
  
}

pdf("plots/ensemble/tree_ensemble.pdf", width=20,height=20)
lapply(list, print)
dev.off()






