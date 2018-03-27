##################################################
### plot heatmap for ari scores in all clusters  ##
##################################################
# Takes the ari_"dataset".rda files , changes data to matrix and plots with pheatmap.
# saves plots in resuls/plots directory

## load libraries
#library(stringr)
library(tibble)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
require(cowplot)

## define the data directories
DATA_DIR <-  "results/run_results"
datatype <- "unfiltered"

# files directories per dataset
plot_pheatmap_ari <- function(files_ari, datatype){
  files_ari <- list(
    kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_kumar2015.rda") ),
    trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_trapnell2014.rda") ),
    koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_koh2016.rda")),
    zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype,"_zhengmix2016.rda")),
    simDataKumar = file.path(DATA_DIR,paste0("ari_single_",datatype,"_simDataKumar.rda")),
    simDataKumar2 = file.path(DATA_DIR,paste0("ari_single_",datatype,"_simDataKumar2.rda"))
    
  )
  
  
  tmp<- list(
    kumar2015 = NULL,
    trapnell2014 =NULL,
    koh2016 = NULL,
    zhengmix2016 =NULL,
    simDataKumar = NULL,
    simDataKumar2 = NULL
    
  )
  
  # load dataset
  tmp <- lapply(files_ari, function(x) get(load(x)))
  # create name vector
  names.methods <- lapply(tmp,names) %>% ldply( data.frame)
  names(names.methods) <- c("data", "method")
  # create table with data , remove the column with the "ground truth" (label)
  tmp <- unlist(tmp)
  tmp <- ldply(tmp, data.frame, .id=) 
  # combine data and names, remove the labels
  tmp <- cbind(tmp , names.methods)%>% subset(!(method=="labels") ) 
  # which methods are missing?
  print(table(tmp$method, tmp$data)==0)
  # table
  tbl <- acast(tmp, data~method, value.var="X..i..")
  # rename
  rownames(tbl) <- c("Koh", "Kumar", "simDataKumar", "simDataKumar2", "Trapnell", "Zheng")
  colnames(tbl) <-c("CIDR", "Linnorm", "pcaReduce", "RaceID", "tSNEkmeans", "SC3", "Seurat", "SIMLR", "SIMLR (large)", "TSCAN", "ZINBWaVE")
  # sort by colmeans
  avgs <-  colMeans(tbl)
  neworder <- order(avgs, decreasing = TRUE)
  tbl <- tbl[,  neworder ]

  ## plot pheatmap
  p <- pheatmap( tbl , color = colorRampPalette(brewer.pal(3, "YlOrRd"))(10),
                 display_numbers = TRUE, number_color = "black", fontsize_number = 11,
                 cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=30, cellheight = 30,
                 number_format="%.2f")$gtable
  return(p)
}

# Appendix
## read in cluster results from Rdata files
# which datatype:  "default", "filtered","unfiltered","optimalk"
p1 <- plot_pheatmap_ari(files_ari, "default")
p2 <- plot_pheatmap_ari(files_ari,"filtered")
p3 <- plot_pheatmap_ari(files_ari, "unfiltered")
#p4 <- plot_pheatmap_ari(files_ari, "optimalk")
p.grid <- plot_grid(p1,p2,p3,ncol = 2, nrow=2,labels="auto" )
save_plot("results/plots/plot_ari_all.pdf", p.grid,
          base_width=13, base_height  =8)
