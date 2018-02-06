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
  tmp <- cbind(tmp , names.methods)%>% subset(!(method== "labels") ) %>% subset(!(method== "SNNCliq") ) 
  # which methods are missing?
  print(table(tmp$method, tmp$data)==0)
  # table
  #xtabs <-  xtabs(X..i..~ data+ method, data=tmp)
  tbl <- acast(tmp, data~method, value.var="X..i..")
  rownames(tbl) <- c( "Koh"   ,    "Kumar" ,    "simDataKumar" ,"simDataKumar2", "Trapnell",  "Zheng"  )
  # order
  tbl[, order(colnames(tbl))]
  
  return(tbl)
  
  
}
#_______________________________________________________________________
plot_ari_ranks <- function(data, main) {
  require(GGally)
  
  # convert results
  ari.mat <- ari2matrix(data )
  rank.mat <- apply(ari.mat, 1,rank, na.last =TRUE, ties.method="min")%>%as.data.frame
  rank.mat <- cbind( method= rownames(rank.mat), rank.mat) %>% as.data.frame
  
  #rank.mat <- melt(rank.mat, id.vars = c("method"), variable.name = "dataset", value.name = ".rank")%>%as.data.frame 
  p1 <- ggparcoord(  rank.mat , columns = c(2:7),groupColumn = 'method', scale = 'globalminmax',missing= "exclude",order="anyClass",title=main )

}
p1 <- plot_ari_ranks(files_ari_filtered , main="rank ARI filtered")
p2 <- plot_ari_ranks(files_ari_unfiltered , main="rank ARI unfiltered")
p3 <- plot_ari_ranks(files_ari_default , main="rank ARI default")
p4 <- plot_ari_ranks(files_ari_optimalk , main="rank ARI optimalk")

p.all <- plot_grid(p1,p2,p3,p4,ncol = 2, labels="auto")
save_plot("results/plots/plot_ari_rank_all.pdf", p.all, base_width=16, base_height = 10)

### Appendix
