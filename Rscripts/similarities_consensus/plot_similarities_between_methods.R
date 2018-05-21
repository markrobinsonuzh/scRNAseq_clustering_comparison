# TODO: in x either change column length or change entry to TRUE
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
  require(purrr)
  require(cowplot)
  
})
### Compare clustering between methods by consensus ARI, k = truenclust
#------------------------------------------------------------------
# load files
res <- readRDS(file = "output/consensus/consensus_clue.rds")
plot_crossmethod_concordance <- function(res, ncluster){
  # group data,
  res.cons <- res %>% dplyr::group_by(dataset, method, k) %>% 
    dplyr::mutate( truenclust =  length( unique(trueclass) )) %>% 
    dplyr::filter(k == ifelse(ncluster==0, truenclust, ncluster) , !method %in% c( "Seurat" ) ) %>% 
    dplyr::select(dataset, method, consensus.clue,k, cell)
  # for Seurat, pick resolution param 
  res.cons.seurat <- res %>% dplyr::group_by(dataset, method, k) %>% 
    dplyr::filter(method %in% c( "Seurat" ) ) %>% 
    dplyr::mutate( truenclust=  length( unique(trueclass) )) %>% 
    dplyr::filter(k== ifelse(ncluster == 0, truenclust, ncluster) )%>%filter(resolution==sample(resolution,1) )%>% 
    dplyr::select(dataset, method, consensus.clue, k, cell)
  res.cons <- bind_rows(res.cons , res.cons.seurat )
  
  # compute  matrix with pairwise ARIs
  l <- vector("list", length=length(unique(res.cons$dataset)))
  names(l) <- unique(res.cons$dataset)
  
  for ( i in unique(res.cons$dataset) ) {
    print(i)
    x <- subset(res.cons, dataset == i )
    m <- matrix(NA, nrow = length(unique(x$method)), ncol = length(unique(x$method)) )
    rownames(m) <- unique(x$method)
    colnames(m) <- unique(x$method)
    for( j in unique(x$method) ) {
      print(j)
      c1 <- subset(x, method==j) 
      c1 <-c1$consensus.clue
      for(u in unique(x$method)){
        print(u)
        c2 <- subset(x, method==u) 
        c2 <-c2$consensus.clue
        m[j,u] <- mclust::adjustedRandIndex(c1,c2)
      }
    }
    l[[i]] <- m
    
  }
  
  # average ARI over datasets
  df <- reshape2::melt(l, value.name="ARI")
  
  df.median <- df%>%group_by(Var1, Var2)%>%dplyr::summarise(median=median(ARI, na.rm = TRUE))%>%arrange(Var1, Var2)%>%
    ungroup()
  
  ## Calculate average area under concordance curve across all data set
  ## instances, for each pair of methods
  
  m.median <- reshape2::acast(df.median, Var1 ~ Var2, value.var = "median") 
  stopifnot(all(rownames( m.median ) == colnames( m.median )))
  
  
  # plot all ARIs, facet by filter and dataset
  heatmap.dataset <-  ggplot( df %>% tidyr::separate(L1, sep = "_", into = c("sce", "filtering", "dataset")) %>%
                                dplyr::select(-sce),
                              aes(x = Var1, y = Var2, fill =ARI))+
    geom_tile(color="white", size=0.5, na.rm =FALSE)+
    facet_grid(filtering ~ dataset, drop=FALSE)+
    scale_fill_viridis(name="ARI", direction=-1 )+
    theme_tufte(base_family="Helvetica")+
    labs(x=NULL, y=NULL, title=paste0("k =", ncluster) ) +
    coord_equal() +
    theme(strip.text.x = element_text(size = 8))+
    theme(axis.text.x=element_text(size=6, angle=90))+
    theme(axis.text.y=element_text(size=6))+
    theme(legend.title=element_text(size=6))+
    theme(legend.title.align=1)+
    theme(legend.text=element_text(size=6))+
    theme(legend.position="right")+
    theme(legend.key.size=unit(0.5, "cm"))+
    theme(legend.key.width=unit(0.2, "cm"))+
    theme(axis.ticks=element_blank())
  
  # plot median ARIs from datasets
  heatmap.median <-  ggplot( df.median,aes(x = Var1, y = Var2, fill = median))+
    geom_tile(color="white", size= 0.5, na.rm =FALSE)+
    scale_fill_viridis(name="ARI", direction=-1 )+
    theme_tufte(base_family="Helvetica")+
    labs(x=NULL, y=NULL, title=paste0("median ARI, k = ", ncluster)) +
    coord_equal() +
    theme(axis.text.x=element_text(size=6, angle=90))+
    theme(axis.text.y=element_text(size=6))+
    theme(legend.title=element_text(size=6))+
    theme(legend.title.align=1)+
    theme(legend.text=element_text(size=6))+
    theme(legend.position="right")+
    theme(legend.key.size=unit(0.5, "cm"))+
    theme(legend.key.width=unit(0.2, "cm"))+
    theme(axis.ticks=element_blank())
  
  #-------------------------------------------------------------------
  # print tree
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
  ## Hierarchical clustering based on 1 - m.median, 
  
  hcl_average <- hclust(as.dist(1-m.median)) 
  
  ## Get all subclusters
  subclusters_average <- get_subclusters(hcl_average)
  
  ## Get all subclusters for individual data set instances
  cmtmp <- df 
  uniqvals <- unique(cmtmp$L1)
  #per dataset i
  subclusters_all <- lapply( uniqvals, function(i) {
    print(i)
    cmtmp2 <- cmtmp %>% dplyr::filter(L1 == i) %>% dplyr::select(Var1, Var2, ARI) # for each dataset
    
    cmtmp2 <-  reshape2::acast(cmtmp2, Var1 ~ Var2, value.var = "ARI")
    stopifnot(all(rownames(cmtmp2) == colnames(cmtmp2)))
    stopifnot(all((cmtmp2 == t(cmtmp2))[!is.na(cmtmp2 == t(cmtmp2))]))
    cmtmp2 <- cmtmp2[apply(cmtmp2,1, function(x){!all(is.na(x) ) } ),]
    cmtmp2 <- cmtmp2[,apply(cmtmp2,2, function(x){!all(is.na(x) ) } )]
    get_subclusters(hclust(as.dist(1-cmtmp2)))
  })
  
  
  ## Get stability values for each subcluster in subclusters_average
  # find datasets with a missing method elemnt in subclusters_all 
  ## dataset with missing seurat
  # finf missings
  
  which_missing <- function(x, invert,  toMatch){
    
    toMatch <- c(paste0( toMatch))
    m <- purrr::map(x, function(x){
      grep( paste(toMatch,collapse="|"), x=x, value=FALSE, invert=FALSE, ignore.case=TRUE) %>%
        is_empty
    }  )
    m <- which(m==invert)
    
    return(m)
  }
  # for this particular datasets (column in w), add the missing method in the particular branch, suchthat the we have a TRUE in w (if the rest is the same)

  x <- sapply(
    subclusters_all, function(w) {
    subclusters_average %in% w })

  p <- rownames(m.median)
  # which missing
  for (i in  p){
    
    x[  which_missing(subclusters_average,invert=FALSE,toMatch = paste(i) ), 
        which_missing(subclusters_all,invert=TRUE, toMatch = paste(i) )  ] <- TRUE
    
  }
  stability_scores <- rowMeans( x )
  
  subcl <- subclusters_average
  stab <- stability_scores
  
  ggt <- ggtree(as.phylo(hcl_average))
  for (m in seq_len(length(subcl))) {
    tryCatch({
      i <- MRCA(ggt, subcl[[m]])
      if ( stab[m] >= 0 )
        ggt$data[i, "label"] <- round(stab[m], 2)
    }, error = function(e) NULL)
  }
  
  trees <- ggt + geom_label2(aes(subset = !isTip, label = label), size = 2) + 
    geom_tiplab(aes(angle = 90), hjust = 1) + 
    ggplot2::scale_x_reverse() + ggplot2::coord_flip() + 
    theme( plot.margin = unit(c(0, 0, 10, 0), "mm") ) + 
    xlim_tree(0.85)+
    ggtitle(paste0("number of cluster=",ncluster) )
  
  # add annotation to tree
  annot <- readr::read_delim("methodchart.csv",";", escape_double = FALSE, trim_ws = TRUE)
  tiporder <- trees$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
  annot <- annot[match(tiporder$label, annot$method), 
                 c("dimension", "clustering", "counts"), drop = FALSE]
  
  #colors for heatmap
  colors <- list(
    dimension = c(PCA = "#C10359", tSNE = "#E0B1C6"),
    clustering = c(hc = "#01368C", graph = "#2F6CCE", kmeans = "#93B8F2", svm = "#D6E4F9"),
    counts = c(counts = "#2EA801", logcounts = "#96E878")
  )
  
  g1 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = dimension)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$dimension) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 7, 0, 7), "mm"))
  g2 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = clustering)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$clustering) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 7, 0, 7), "mm"))
  g3 <- ggplot(annot) + geom_tile(aes(x = 1:nrow(annot), y = 1, fill = counts)) + 
    geom_vline(xintercept = (1:(nrow(annot) - 1)) + 0.5, linetype = "solid", color = "white", size = 0.25) + 
    scale_fill_manual(values = colors$counts) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 7, 0, 7), "mm"))
  
  gds <- plot_grid(get_legend(g1 + 
                                theme(legend.position = "left") + 
                                guides(fill = 
                                         guide_legend(ncol = 2,
                                                      title = "dimension", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 12,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 10,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   get_legend(g2 + 
                                theme(legend.position = "left") + 
                                guides(fill = 
                                         guide_legend(ncol = 1,
                                                      title = "clustering", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 12,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 10,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   get_legend(g3 + 
                                theme(legend.position = "left") + 
                                guides(fill = 
                                         guide_legend(ncol = 1,
                                                      title = "counts", title.position = "top",
                                                      override.aes = list(size = 1.5),
                                                      title.theme = element_text(size = 12,
                                                                                 angle = 0),
                                                      label.theme = element_text(size = 10,
                                                                                 angle = 0),
                                                      keywidth = 1, default.unit = "cm"))),
                   
                   nrow = 1, rel_widths = c(1.25, 1, 1, 1))
  
  ann.tree <- plot_grid(trees + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")), 
            g1 + theme(legend.position = "none"), 
            g2 + theme(legend.position = "none"), 
            g3 + theme(legend.position = "none"), 
            
            gds + theme(plot.margin = unit(c(0, 0, 0, 15), "mm")), 
            rel_heights = c(6, 0.75, 0.75 ,0.75 ,0.75, 2), ncol = 1)
  
  # store plots in plot.list 
  plot.list <- list()
  plot.list[[1]] <- heatmap.dataset
  plot.list[[2]] <- heatmap.median
  plot.list[[3]] <- ann.tree
  return(plot.list)
  
}
# no Zheng
res1 <- res%>%filter(!method %in% c("RaceID"))%>%filter( !dataset %in% grep("Zheng",unique(dataset), value=TRUE) )
l<- as.list(c(0,3:10))
# plot 
list.trees <- lapply(l, function(x) {plot_crossmethod_concordance(res1, ncluster=x)} )
pdf("plots/similarities_between_methods/similarities_median_nozheng.pdf", width=15, height=10)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',2))
dev.off()
pdf("plots/similarities_between_methods/similarities_dataset_nozheng.pdf", width=30, height=30)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',1), ncol=2)
dev.off()
pdf("plots/similarities_between_methods/similarities_tree_nozheng.pdf", width=15, height=10)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',3))
dev.off()


# Zheng
res2 <- res%>%filter(!method %in% c("RaceID"))%>%filter( dataset %in% grep("Zheng",unique(dataset), value=TRUE) )
l<- as.list(c(0, 3:10))
# plot 
list.trees <- lapply(l, function(x) {plot_crossmethod_concordance(res2, ncluster=x)} )
pdf("plots/similarities_between_methods/similarities_median_zheng.pdf", width=15, height=10)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',2))
dev.off()

pdf("plots/similarities_between_methods/similarities_dataset_zheng.pdf", width=30, height=30)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',1), ncol=2)
dev.off()

pdf("plots/similarities_between_methods/similarities_tree_zheng.pdf", width=15, height=10)
cowplot::plot_grid(plotlist = sapply(list.trees ,'[',3))
dev.off()

