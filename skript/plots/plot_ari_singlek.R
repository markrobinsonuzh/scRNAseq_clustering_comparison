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

## define the data directories
DATA_DIR <-  "results/run_results"


## read in cluster results from Rdata files
# which dataset
datatype <- "unfiltered"
# files directories per dataset
files_ari <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_koh2016.rda")),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype,"_zhengmix2016.rda")),
  simDataKumar = file.path(DATA_DIR,paste0("ari_single_",datatype,"_simDataKumar.rda"))
)
tmp<- list(
  kumar2015 = NULL,
  trapnell2014 =NULL,
  koh2016 = NULL,
  zhengmix2016 =NULL,
  simDataKumar = NULL
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
  #xtabs <-  xtabs(X..i..~ data+ method, data=tmp)
  tbl <- acast(tmp, data~method, value.var="X..i..")

  ## plot pheatmap
  pheatmap( tbl , color = colorRampPalette(brewer.pal(3, "YlOrRd"))(10),
         display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
         cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=30, cellheight = 30,
         main = "ARI scores", 
         number_format="%.2f", filename=paste0("results/plots/plot_ari_",datatype,".pdf"))
  

# Appendix

  