##################################################
### plot heatmap for differences  in ari scores ##
##################################################
# Takes the ari_"dataset".rda files , changes data to matrix and plots the differences with pheatmap.
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
# files directories per dataset
files_ari_filtered <- list(
  kumar2015 = file.path(DATA_DIR, "ari_single_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "ari_single_trapnell2014.rda"),
  koh2016 = file.path(DATA_DIR, "ari_single_koh2016.rda"),
  zhengmix = file.path(DATA_DIR, "ari_single_zhengmix2016.rda"),
  simDataKumar = file.path(DATA_DIR, "ari_single_simDataKumar.rda")
)
files_ari_unfiltered <- list(
  kumar2015 = file.path(DATA_DIR, "ari_single_unfiltered_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "ari_single_unfiltered_trapnell2014.rda"),
  koh2016 = file.path(DATA_DIR, "ari_single_unfiltered_koh2016.rda"),
  zhengmix = file.path(DATA_DIR, "ari_single_unfiltered_zhengmix2016.rda"),
  simDataKumar = file.path(DATA_DIR, "ari_single_unfiltered_simDataKumar.rda")
)

#_______________________________________________________________________
# function to convert stored Ari results to matrix
ari2matrix <- function(files.ari) {
  tmp<- list(
    kumar2015 = NULL,
    trapnell2014 =NULL,
    koh2016 = NULL,
    zhengmix2016 =NULL,
    simDataKumar = NULL
  )
  tmp <- lapply(files.ari, function(x) get(load(x)))
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
  # order
  tbl[, order(colnames(tbl))]
  
  return(tbl)
  
  
}
#_______________________________________________________________________

# convert results
ari.filtered <- ari2matrix(files_ari_filtered )
ari.unfiltered <- ari2matrix(files_ari_unfiltered )
# subset methods
match <- intersect( colnames(ari.filtered), colnames(ari.unfiltered) )

ari.filter.sub <- ari.filtered[, match ] 
ari.unfilter.sub <- ari.unfiltered[, match]
# differences
ari.diff <- ari.filter.sub - ari.unfilter.sub


## plot pheatmap
pheatmap( ari.diff , color = colorRampPalette(brewer.pal(3, "RdYlGn"))(10),
          display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
          cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=30, cellheight = 30,
          main = "difference ARI by filtered vs. unfiltered", 
          number_format="%.2f", filename=paste0("results/plots/plot_ari_diff_filterd_unfilterd.pdf"))


# Appendix

