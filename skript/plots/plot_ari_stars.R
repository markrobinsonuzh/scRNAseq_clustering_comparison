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
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar2.rda"))
  
)
files_ari_unfiltered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar2.rda"))
  
)
files_ari_optimalk <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar2.rda"))
  
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
  # order
  tbl[, order(colnames(tbl))]
  
  return(tbl)
  
  
}
#_______________________________________________________________________
plot_ari_stars <- function(data, main) {
  require(GGally)
  
  # convert results
  ari.mat <- ari2matrix(data )
  p1 <- stars(t(ari.mat), key.labels = rownames(ari.mat),  key.loc = c(8,2), main=main )
  #p2 <- stars(ari.mat, labels = rownames(ari.mat),main=main)
}
pdf("results/plots/plot_ari_stars_all.pdf")
p1 <- plot_ari_stars(files_ari_filtered ,main="ARI filtered")
p2 <- plot_ari_stars(files_ari_unfiltered , main="rank ARI unfiltered")
p3 <- plot_ari_stars(files_ari_default , main="rank ARI default")
p4 <- plot_ari_stars(files_ari_optimalk , main="rank ARI optimalk")
dev.off()
plot_ari_stars <- function(data, main) {
  require(GGally)
  
  # convert results
  ari.mat <- ari2matrix(data )
  #p1 <- stars(t(ari.mat), key.labels = rownames(ari.mat),  key.loc = c(8,2), main=main )
  p2 <- stars(ari.mat, key.labels = colnames(ari.mat),  key.loc = c(5,2), main=main)
}
pdf("results/plots/plot_ari_stars_perdata.pdf")
p1 <- plot_ari_stars(files_ari_filtered ,main="ARI filtered")
p2 <- plot_ari_stars(files_ari_unfiltered , main="rank ARI unfiltered")
p3 <- plot_ari_stars(files_ari_default , main="rank ARI default")
p4 <- plot_ari_stars(files_ari_optimalk , main="rank ARI optimalk")
dev.off()

### Appendix
