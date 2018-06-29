#tSNE plots

DIR <- "data/"
FILTERINGS  <-  c("Expr", "M3Drop" ,"HVG")
PCTKEEP  <-  c("10")
DATASETS <-  c("Kumar", "Trapnell" ,"Koh" ,"SimKumar4easy" ,"SimKumar4hard", "SimKumar8hard", "KohTCC", "TrapnellTCC" ,"KumarTCC","Zhengmix4eq" ,"Zhengmix4uneq" ,"Zhengmix8eq")

DATASETS_FILTERED <- outer(paste0("sce_filtered",FILTERINGS,PCTKEEP), paste0(DATASETS, ".rds"), FUN = paste,sep = "_")
list <- as.list(paste(paste0(DIR,"sce_filtered",FILTERINGS ,PCTKEEP),DATASETS_FILTERED, sep = "/"))
full_data <- lapply(list, function(x) readRDS(x)  )
names(full_data) <- DATASETS_FILTERED
# Clustering results
res_summary <- readRDS(file="output/clustering_summary/clustering_summary.rds")

#plot tSNE
plots.tsne <- plyr::laply(res, function(x)scater::plotTSNE(x, colour_by = "phenoid", title=NULL) )

data <- c("sce_filteredExpr10_Zhengmix4eq", "sce_filteredExpr10_Zhengmix4eq","sce_filteredExpr10_Zhengmix4uneq", "sce_filteredExpr10_Trapnell")
nrun <- c(1)
meth <- c("FlowSOM", "CIDR","TSCAN", "SC3")

plot_tSNE <- function(res_summary, full_data, meth,data,run){
  require(dplyr)
  require(ggplot2)

  res <- res_summary%>%dplyr::filter( dataset%in%data, method%in%meth, run==nrun )%>%
    group_by(dataset)%>%
    filter(k==length( unique(trueclass)) )%>%
    select(cluster, trueclass)
  sel.data <- full_data[[paste0(data,".rds")]]
  df_to_plot <- data.frame(SingleCellExperiment::reducedDim(sel.data,"TSNE"),
                         colnames(sel.data),res$cluster,res$trueclass )
  p <- ggplot(df_to_plot, 
        aes(x = X1 , y = X2)) + 
        xlab("dimension 1") + 
        ylab("dimension 2") + 
        geom_rug(colour = "gray20", 
        alpha = 0.65) + 
        theme_bw() +
        geom_point(aes(colour = res.cluster, 
                 shape =  res.trueclass), 
                 alpha = 0.65)+
        guides(shape = guide_legend(title = "true class"),
               colour= guide_legend(title = "method clustering")) +
        labs(title=paste(data,meth))
      
  return(p)
}
p.list <- lapply(c(1:4), function(x){
    plot_tSNE(res_summary, full_data, meth[x],data[x],nrun)
})
#save plots
plot.grid <- cowplot::plot_grid(plotlist=p.list ,labels = "auto")
cowplot::ggsave(plot.grid, filename = "plot_facets_clustering.pdf", width=20, height=15)
