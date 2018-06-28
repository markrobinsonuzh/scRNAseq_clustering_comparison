#------------------------------------------------
# Venn diagrams of genes for the filtered datasets
#------------------------------------------------
# libs
library('VennDiagram')
library(gridExtra)
# repos
DIR <- "" # repo of datasets
FILTERINGS  <-  c("Expr", "M3Drop" ,"HVG")
PCTKEEP  <-  c("10")
DATASETS <-  c("Kumar", "Trapnell" ,"Koh" ,"SimKumar4easy" ,"SimKumar4hard", "SimKumar8hard", "KohTCC", "TrapnellTCC" ,"KumarTCC","Zhengmix4eq" ,"Zhengmix4uneq" ,"Zhengmix8eq")

DATASETS_FILTERED <- outer(paste0("sce_filtered",FILTERINGS,PCTKEEP), paste0(DATASETS, ".rds"), FUN = paste,sep = "_")
list <- as.list(paste(paste0(DIR,"sce_filtered",FILTERINGS ,PCTKEEP),DATASETS_FILTERED, sep = "/"))
# load data
full_data <- lapply(list, function(x) readRDS(x)  )
names(full_data) <- DATASETS_FILTERED
# extract genes
red.data <- lapply(full_data , function(x){
  genes <- rownames(x)
})
# filter, rename cols
genes.list <- lapply(
  lapply( paste0(DATASETS, ".rds") , function(x){
  purrr::keep(red.data, grepl(x, names(red.data))==TRUE)
}), setNames, nm= c("Expr10","M3Drop","HVG"))
# venn diagrams
list <- lapply(genes.list ,function(x){
 VennDiagram::venn.diagram(x, filename = NULL, 
                           fill=c("darkmagenta", "darkblue", "darkred"),alpha=c(0.3,0.3,0.3),
                           cex=0.8,cat.fontface=1, margin=0.22,cat.cex=c(1,1,1),
                           euler.d=TRUE,
                           ext.text=TRUE,
                           lwd=1,
                           cat.default.pos="outer",
                           cat.dist = c(0.05, 0.05,0.05)
                           )
  })
# reformat gList to gTree
list.g3 <-  lapply(list, function(x)gTree(children=x))
# plot grid
lab <- DATASETS
grid <- cowplot::plot_grid(plotlist = list.g3, ncol=3, labels =lab )
cowplot::ggsave("plots/manuscript/filterings_shared_genes.pdf", grid, width=15, height=10)

