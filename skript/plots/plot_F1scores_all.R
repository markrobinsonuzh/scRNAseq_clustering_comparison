##################################################
### plot heatmap for F1 scores in all clusters  ##
##################################################
# Takes the f1_"dataset".rda files , changes data to matrix and plots with pheatmap.
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
DATA_DIR <-  "results/run_results"

datatype <-  "default"
# datasets


## function: reads in cluster results from Rdata files, create heatmaps
# Input: files_f1 = directory of files, datatype = type of runmode
# Output: cowplot grid object with heatmaps

plot_pheatmaps_f1 <- function(files_f1,datatype ){
  files_f1 <- list(
    kumar2015 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_kumar2015.rda") ),
    trapnell2014 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_trapnell2014.rda") ),
    zhengmix2016=file.path(DATA_DIR, paste0("f1_single_",datatype,"_zhengmix2016.rda") ),
    simDataKumar =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar.rda")),
    simDataKumar2 =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar2.rda")),
    koh2016 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_koh2016.rda") )
    
  )
  p<- vector("list", length(files_f1))
  names(p) <- names(files_f1)
  # reorder
  
  for ( i in names(files_f1)){
  # load dataset
  tmp <- lapply(files_f1[[i]], function(x)get(load(x)))
  
  ## plot pheatmap
  # create table with data , remove the column with the "ground truth" (label)
  tmp <- ldply(tmp[[1]], data.frame) %>% select(.id, f1, labels, n.cluster) %>% subset(!(.id=="labels") )
  # tmp <- daply(tmp, .(labels, .id), function(x) x$f1) old way 
  tmp <- xtabs( f1~labels+.id, tmp, addNA=FALSE )
  
  
  # remove NA labels
  tmp <- tmp[!is.na( rownames(tmp) ), ]
  # assign a zero to undefined entries
  tmp <- ifelse(is.na(tmp)==TRUE, NA, tmp)
  # remove the label column from the table
  
  ##plot it

  p[[i]] <- pheatmap(tmp , color = colorRampPalette(brewer.pal(3, "YlGnBu"))(10),
           display_numbers = TRUE, number_color = "black", fontsize_number = 8, 
           cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=20, cellheight = 20,
           main = paste0(names(files_f1[i])), 
           width = 6, 
           height = 6,
           number_format="%.2f")$gtable
  }
  p <-  plot_grid(plotlist=p,ncol = 2, nrow=3,labels="auto" )
  save_plot(filename= paste0("results/plots/plot_f1_",datatype,".pdf") , plot=p, base_height=10, base_width = 8)
  
  }

# Create plots in Grid
p1 <-  plot_pheatmaps_f1(files_ari, "default")
p2 <-  plot_pheatmaps_f1(files_ari,"filtered")
p3 <-  plot_pheatmaps_f1(files_ari, "unfiltered")
p4 <-  plot_pheatmaps_f1(files_ari, "optimalk")
# save to pdf
