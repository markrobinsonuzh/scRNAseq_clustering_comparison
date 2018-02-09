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

## define the data directories
DATA_DIR <-  "results/run_results"

# which dataset
datatype <- "filtered"
## read in cluster results from Rdata files
# files directories per dataset
files_f1 <- list(
  kumar2015 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("f1_single_",datatype,"_koh2016.rda") ),
  zhengmix2016=file.path(DATA_DIR, paste0("f1_single_",datatype,"_zhengmix2016.rda") ),
  simDataKumar =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar.rda")),
  simDataKumar2 =file.path(DATA_DIR, paste0( "f1_single_",datatype,"_simDataKumar2.rda"))
  
)

row_labels <- list(
  kumar2015 = c("Dgcr8serum","v6.52i","v6.5serum") ,
  trapnell2014 = c("T0","T24","T48") ,
  koh2016 =  c("APS","D2LtM","D3GARPpCrdcM", "D5CntrlDrmmtm", "DLL1pPXM", "ESMT","MPS","Sclrtm","D2.25_Smtmrs","H7hESC") ,
  zhengmix2016= c("b.cell","cd14.monocytes","naive.cytotoxic", "regulatory.t") ,
  simDataKumar = c("G1","G2","G3", "G4") ,
  simDataKumar2 = c("G1","G2","G3", "G4") 
  
)

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
  pheatmap(tmp , color = colorRampPalette(brewer.pal(3, "YlGnBu"))(10),
           display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
           cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=30, cellheight = 30,
           main = paste0(names(files_f1[i])), 
           labels_row = row_labels[[i]]  ,
           width = 6, 
           height = 6,
           number_format="%.2f", filename= paste0("results/plots/plot_f1_",datatype,"_",names(files_f1[i]),".pdf")
  )
  
}

# Appendix
