
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
res_ensmbl <- readRDS(file="output/ensemble/clustering_ensemble_allmethods2.rds")
res_summary <- readRDS(file="output/clustering_summary/clustering_summary.rds")


# ------------------------------------
# compute ARI, no of unique clusters
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
sum_ensmbl$run <- as.integer( sum_ensmbl$run )
# join datasets
sum_all<- left_join(sum_ensmbl , sum_summary, by=c("filtering", "dataset" ,"methone"="method", "run"))
#
#rename, compute difference
colnames(sum_all) <- c("filtering",
                       "dataset",
                       "method.ensmbl",
                       "methone",
                       "methtwo",
                       "methtree",
                       "run",
                       "ARIensmbl",
                       "truenclust",
                       "k",
                       "ARIone",
                       "truenclust.y")

sum_all <- sum_all %>% mutate(ARIdiff=ARIone-ARIensmbl)

# compute median
res_median <- sum_all%>% group_by(dataset, filtering, method.ensmbl, k, truenclust, methone, methtwo, methtree ) %>% dplyr::summarize(ARIensmbl=median(ARIensmbl), ARIone=median(ARIone), ARIdiff=median(ARIdiff) )
list <- vector("list",length( unique((res_median$dataset)) ))

for(i in unique(res_median$dataset) ){ 
res_median.sub <- filter(res_median, k==truenclust, dataset==i)

# work with one dataset
res_median.sub$pathString <- paste("dataset", 
                                   res_median.sub$dataset,
                                   res_median.sub$filtering,
                                   res_median.sub$methone,
                                   res_median.sub$method.ensmbl,
                                   sep="/")
# create tree using data.tree
ARItree <- as.Node(res_median.sub)
#first simple tree

#color=c("green","blue" )
#useRtreeList <- ToListExplicit(ARItree, unname = TRUE)
#radialNetwork( useRtreeList, opacity=0.74, nodeColour=color, nodeStroke = "blue")
#?radialNetwork

# convert to dataframes
#x <- ToDataFrameTree(ARItree , "ARIdiff")
#y <-  ToDataFrameTree(ARItree , "ARIone" )

annotation <- ToDataFrameTypeCol(ARItree, "ARIdiff", "ARIensmbl", "ARIdiff")

# format as phylo object
gh <- as.phylo.Node(ARItree)
# plot tree
# add annotataion
annotation <-  annotation%>%select("level_5", "ARIdiff")
p <- ggtree(gh, layout = 'unrooted', branch.length = "branch.length") %<+% annotation
# change names for 2nd layer and 3rd layer
p$data$label <- gsub('.*\\.',"", p$data$label)
p$data$branch.length <- 10

p2 <- p  +
  geom_text2(aes(color=ARIdiff, label=round( ARIdiff,2)), size=3, nudge_x = -2.0, na.rm = TRUE) +
  geom_tiplab2(size=2, hjust=1.0 , geom='text', align=FALSE)+
  geom_nodelab( size=2, hjust=1, geom='label')
print(p2)
list[[i]] <- p2
}
pdf("tree_ensemble_dataset.pdf")
lapply(list, print)
dev.off()


# Stuff
q()

#________________________________________________________________________________________________
#adjust branch length

p.new <- p$data %>% filter(isTip==TRUE)
p.new$label==annotation$level_5

p$data$branch.length <- annotation$ARIdiff[ match(p$data$label, annotation$level_5) ]
p$data$branch.length <- p$data$branch.length*100


#-------------------------------------------------------------------------------------

# creat a facet wrap

plot.tree <- function(dataset, res_median){
res_median.sub <- filter( res_median, dataset==dataset, k==truenclust )

# work with one dataset
res_median.sub$pathString <- paste("dataset", 
                                   res_median.sub$dataset,
                                   res_median.sub$filtering,
                                   res_median.sub$methone,
                                   res_median.sub$method.ensmbl,
                                   sep="/")
# create tree using data.tree
ARItree <- as.Node(res_median.sub)
#first simple tree
return(ARItree)

}
# as facet
list <- list(Koh="Koh", Kumar="Kumar")

p.list <- lapply(list, plot.tree, res_median=res_median)
class(p.list) <- "multiPhylo"

ggtree(p.list) + facet_wrap() + ggtitle("Many trees. Such phylogenetics. Wow.")


set.seed(42)
trees <- lapply(rep(c(10, 25, 50, 100), 3), rtree)
class(trees) <- "multiPhylo"
class(trees)
ggtree(trees) + facet_wrap(~.id, scale="free", ncol=4) + ggtitle("Many trees. Such phylogenetics. Wow.")
library("ape")
data(chiroptera)
groupInfo <- split(chiroptera$tip.label, gsub("_\\w+", "", chiroptera$tip.label))
chiroptera <- groupOTU(chiroptera, groupInfo)
ggtree(chiroptera, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))

