##################################################
### plot heatmap for F1 scores in all clusters  ##
##################################################
# Takes the f1_"dataset".rda files , changes data to matrix and plots beeswarm.
# saves plots in resuls/plots directory

## load libraries
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

## define the data directories
source("skript/helper_files/Helper_functions.R")

DATA_DIR <-  "results/run_results"

names.method <-  c("CIDR", "Linnorm", "pcaReduce", "RaceID", "tSNEkmeans", "SC3", "Seurat", "SIMLR", "SIMLR (large)", "TSCAN", "ZINBWaVE")

## function: reads in cluster results from Rdata files, create heatmaps
# Input: files_f1 = directory of files, datatype = type of runmode, name.method = names of method, as char vector
# Output: geom_point ggplot plots with F1 scores as datapoints for each cluster, size by cluster size
plot_point_f1 <- function(files_f1, datatype, names.method ){
  files_f1 <- list(
    kumar2015 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_kumar2015.rda") ),
    trapnell2014 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_trapnell2014.rda") ),
    zhengmix2016=file.path(DATA_DIR, paste0("f1_single_",datatype,"_zhengmix2016.rda") ),
    simDataKumar =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar.rda")),
    simDataKumar2 =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar2.rda")),
    koh2016 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_koh2016.rda") )
    
  )

  p <- vector("list", length(files_f1))
  names(p) <- names(files_f1)
  # reorder
  
  for ( i in names(files_f1)){
    # load dataset
    tmp <- lapply(files_f1[[i]], function(x)get(load(x)))
    
    # create table with data , remove the column with the "ground truth" (label)
    tmp <- ldply(tmp[[1]], data.frame) %>% select(.id, f1, labels, n.cluster) %>% subset(!(.id=="labels") )
   # assign a zero to missing clusters, so that they visible in plot
     tmp$f1[is.na(tmp$f1)] <- 0
   # rename method name
     tmp$.id <- factor(tmp$.id)
     levels(tmp$.id) <- names.method
    ## plot dotplot
    p[[i]] <- ggplot( tmp, aes(x=.id, y=f1))+
      geom_point(aes(size=n.cluster, color=factor(labels, exclude = NULL)))+
      labs(x="", y="F1",  color="cluster" , size="no of cells")+
      theme_gray()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      guides(fill = "none")+
      coord_cartesian(ylim = c(0, 1))

  }
 p.list <- plot_grid(plotlist = p, ncol=2, nrow=3,labels="auto" )
 save_plot(filename = paste0("results/plots/plot_f1_point_",datatype,".pdf"), plot=p.list,   base_width=13, base_height  =15)
}

# Save plots
 plot_point_f1( files_f1 ,"default", names.method)
 plot_point_f1(files_f1,"filtered",  names.method)
 plot_point_f1(files_f1, "unfiltered",  names.method)
 plot_point_f1(files_f1, "optimalk", names.method)

 