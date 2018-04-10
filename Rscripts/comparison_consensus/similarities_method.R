#------------------------------------------------------------
# Method similarities; based on the consensus, k=truenclust
#____________________________________________________________

suppressPackageStartupMessages({
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(clusterExperiment)
  require(clue)
  require(multidplyr)
  require(ggplot2)
  require(viridis)
  require(ggthemes)
  require(pheatmap)
  require(reshape2)
  require(mclust)
  require(RColorBrewer)
  require(ggtree)
  
})
### Compare clustering between methods by consensus ARI, k = truenclust
#------------------------------------------------------------------
# load files

res <- readRDS(file = "../../output/consensus/consensus_clue.rds")
# group data
res.cons <- res %>% dplyr::group_by(dataset, method, k)%>% 
  dplyr::mutate( truenclust=  length( unique(trueclass) )) %>% 
  dplyr::filter(k==truenclust) %>% dplyr::select(dataset, method, consensus.clue,k)

plot_crossmethod_concordance <- function(res.cons){
# compute pairwise ARIs
l <- vector("list", length=length(unique(res.cons$dataset)))
names(l) <- unique(res.cons$dataset)
for ( i in unique(res.cons$dataset) ) {
  x <- subset(res.cons,dataset == i )
  m <- matrix(NA, nrow = length(unique(x$method)), ncol = length(unique(x$method)) )
  rownames(m) <- unique(x$method)
  colnames(m) <- unique(x$method)
  for( j in unique(x$method) ) {
    c1 <- subset(x, method==j) 
    c1 <-c1$consensus.clue
    for(u in unique(x$method)){
      c2 <- subset(x, method==u) 
      c2 <-c2$consensus.clue
      m[j,u] <- mclust::adjustedRandIndex(c1,c2)
    }
  }
  l[[i]] <- m
  
}
# average ARI over datasets
df <- reshape2::melt(l, value.name="ARI")

df.median <- df%>%group_by(Var1, Var2)%>%dplyr::summarise(median=median(ARI))%>%ungroup()

## Get all subclusters from an hclust object, from  https://github.com/csoneson/conquer_comparison/blob/master/scripts/help_function_crossmethod_concordance.R
get_subclusters <- function(hcl) {
  m <- hcl$merge
  labs <- hcl$labels
  L <- list()
  for (i in seq_len(nrow(m))) {
    tmp <- c()
    if (m[i, 1] < 0) tmp <- c(tmp, labs[-m[i, 1]])
    else tmp <- c(tmp, L[[m[i, 1]]])
    if (m[i, 2] < 0) tmp <- c(tmp, labs[-m[i, 2]])
    else tmp <- c(tmp, L[[m[i, 2]]])
    L[[i]] <- sort(tmp)
  }
  L
}
  ## Calculate average area under concordance curve across all data set
  ## instances, for each pair of methods
  
  m.median <- reshape2::acast(df.median, Var1 ~ Var2, value.var = "median")
  stopifnot(all(rownames( m.median ) == colnames( m.median )))
  ## Hierarchical clustering based on 1 - cmcons, 
  
  hcl_average <- hclust(as.dist(1-m.median))
  
  ## Get all subclusters
  subclusters_average <- get_subclusters(hcl_average)
  
  ## Get all subclusters for individual data set instances
  cmtmp <- df 
  uniqvals <- unique(cmtmp$L1)
  subclusters_all <- lapply(uniqvals, function(i) {
    cmtmp2 <- cmtmp %>% dplyr::filter(L1 == i) %>% dplyr::select(Var1, Var2, ARI) # for each dataset
   
    cmtmp2 <-  reshape2::acast(cmtmp2, Var1 ~ Var2, value.var = "ARI")
    stopifnot(all(rownames(cmtmp2) == colnames(cmtmp2)))
    stopifnot(all((cmtmp2 == t(cmtmp2))[!is.na(cmtmp2 == t(cmtmp2))]))
    get_subclusters(hclust(as.dist(1-cmtmp2)))
  })
  
  ## Get stability values for each subcluster in subclusters_average
  stability_scores <- rowMeans( x <- sapply(subclusters_all, function(w) {
    subclusters_average %in% w
  }))
  # plot all ARIs, facet by filter and dataset
  print( ggplot( df %>% tidyr::separate(L1, sep = "_", into = c("sce", "filtering", "dataset")) %>%
                  dplyr::select(-sce),
                aes(x = Var1, y = Var2, fill =ARI))+
    geom_tile(color="white", size=0.5, na.rm =FALSE)+
    facet_grid(filtering ~ dataset, drop=FALSE)+
    scale_fill_viridis(name="medianARI", direction=-1 )+
    theme_tufte(base_family="Helvetica")+
    labs(x=NULL, y=NULL, title=" ARI, k = truenclust") +
    coord_equal() +
    theme(strip.text.x = element_text(size = 8))+
    theme(axis.text.x=element_text(size=8, angle=90))+
    theme(axis.text.y=element_text(size=8))+
    theme(legend.title=element_text(size=8))+
    theme(legend.title.align=1)+
    theme(legend.text=element_text(size=6))+
    theme(legend.position="right")+
    theme(legend.key.size=unit(2, "cm"))+
    theme(legend.key.width=unit(0.5, "cm"))+
    theme(axis.ticks=element_blank())
  )
  # plot median ARIs from datasets
  print( ggplot(df.median,aes(x = Var1, y = Var2, fill = median))+
    geom_tile(color="white", size=0.5, na.rm =FALSE)+
    scale_fill_viridis(name="medianARI", direction=-1 )+
    theme_tufte(base_family="Helvetica")+
    labs(x=NULL, y=NULL, title="median ARI, k = truenclust") +
    coord_equal() +
    theme(axis.text.x=element_text(size=8, angle=90))+
    theme(axis.text.y=element_text(size=8))+
    theme(legend.title=element_text(size=8))+
    theme(legend.title.align=1)+
    theme(legend.text=element_text(size=6))+
    theme(legend.position="right")+
    theme(legend.key.size=unit(2, "cm"))+
    theme(legend.key.width=unit(0.5, "cm"))+
    theme(axis.ticks=element_blank())
    )
  subcl <- subclusters_average
  stab <- stability_scores
  
   ggt <- ggtree(as.phylo(hcl_average))
   for (m in seq_len(length(subcl))) {
     tryCatch({
       i <- MRCA(ggt, subcl[[m]])
       if (stab[m] >= 0.001)
         ggt$data[i, "label"] <- round(stab[m], 2)
     }, error = function(e) NULL)
   }

   ggt <- ggt + geom_label2(aes(subset = !isTip, label = label), size = 2) + 
     geom_tiplab(aes(angle = 90), hjust = 1) + 
     ggplot2::scale_x_reverse() + ggplot2::coord_flip() + 
     theme(plot.margin = unit(c(0, 0, 10, 0), "mm")) + 
     xlim_tree(0.85)
   tiporder <- ggt$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
   
   

}
# plot 
pdf("../../plots/consensus/method_similarities.pdf", width=15, height=10)
plot_crossmethod_concordance(res.cons)
dev.off()

date()
sessionInfo()


