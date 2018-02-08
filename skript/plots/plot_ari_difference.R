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
library(cowplot)

## define the data directories
DATA_DIR <-  "results/run_results"
datatype1<- "filtered"
datatype2 <- "default"
datatype3 <- "unfiltered"
datatype4 <- "optimalk"


## read in cluster results from Rdata files
# files directories per dataset
files_ari_filtered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar2.rda"))
  
)
files_ari_default <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_simDataKumar2.rda"))
  
)
files_ari_unfiltered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar2.rda"))
  
)
files_ari_optimalk <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar2.rda"))
  
)

#_______________________________________________________________________
# function to convert stored Ari results to matrix
ari2matrix <- function(files.ari) {
  tmp<- list(
    kumar2015 = NULL,
    trapnell2014 =NULL,
    koh2016 = NULL,
    zhengmix2016 =NULL,
    simDataKumar = NULL,
    simDataKumar2 = NULL
    
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
  # rename
  rownames(tbl) <- c("Koh", "Kumar", "simDataKumar", "simDataKumar2", "Trapnell", "Zheng")
  colnames(tbl) <-c("CIDR", "Linnorm", "pcaReduce", "RaceID", "tSNEkmeans", "SC3", "Seurat", "SIMLR", "SIMLR (large)", "TSCAN", "ZINBWaVE")

  return(tbl)
  
  
}
#_______________________________________________________________________
plot_ari_differences <- function(files_ari_datatype1, files_ari_datatype2, main) {
# convert results
ari.datatype1 <- ari2matrix(files_ari_datatype1 )
ari.datatype2 <- ari2matrix(files_ari_datatype2 )

# subset methods
match <- intersect( colnames(ari.datatype2), colnames(ari.datatype2 ) )

ari.datatype1.sub <- ari.datatype1[, match ] 
ari.datatype2.sub <- ari.datatype2 [, match]
# differences
ari.diff <- ari.datatype1.sub  - ari.datatype2.sub


## plot pheatmap
pheatmap( ari.diff , color = colorRampPalette(brewer.pal(3, "PiYG"))(10),
          display_numbers = TRUE, number_color = "black", fontsize_number = 8, 
          cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=25, cellheight = 25,
          main = main, 
          number_format="%.2f")

}
p1 <- plot_ari_differences(files_ari_filtered , files_ari_unfiltered , main="difference ARI: filtered vs. unfiltered")
p2 <- plot_ari_differences(files_ari_filtered, files_ari_default , main="difference ARI: filtered vs.default")
p3 <- plot_ari_differences(files_ari_default, files_ari_unfiltered, main="difference ARI: default vs. unfiltered" )
p4 <- plot_ari_differences(files_ari_optimalk,files_ari_filtered , main="difference ARI: optimalk vs. filtered" )
p.filtdef <- plot_grid(p1$gtable, p2$gtable,ncol = 2, labels="auto")
save_plot("results/plots/plot_ari_diff_filtdef.pdf", p.filtdef , base_width=12, base_height = 5)

p.all <- plot_grid(p1$gtable, p2$gtable, p3$gtable, p4$gtable,ncol = 1, labels="auto")
save_plot("results/plots/plot_ari_diff_all.pdf", p.all, base_width=7, base_height = 14)

