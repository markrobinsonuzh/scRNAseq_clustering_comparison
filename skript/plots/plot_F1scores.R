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


## read in cluster results from Rdata files
# files directories per dataset
files_f1 <- list(
  kumar2015 = file.path(DATA_DIR, "f1_single_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "f1_single_trapnell2014.rda"),
  koh2016 = file.path(DATA_DIR, "f1_single_koh2016.rda"),
  zhengmix2016=file.path(DATA_DIR, "f1_single_zhengmix2016.rda"),
  simDataKumar =file.path(DATA_DIR, "f1_single_simDataKumar.rda")
)

for ( i in names(files_f1)){
# load dataset
tmp <- lapply(files_f1[[i]], function(x)get(load(x)))

## plot pheatmap
# create table with data , remove the column with the "ground truth" (label)
tmp <- ldply(tmp[[1]], data.frame) %>% select(.id, f1, labels, n.cluster) %>% subset(!(.id=="labels") )
# tmp <- daply(tmp, .(labels, .id), function(x) x$f1) old way 
tmp <- xtabs( f1~labels+.id, tmp, addNA=TRUE )


# remove NA labels
tmp <- tmp[!is.na( rownames(tmp) ), ]
# assign a zero to undefined entries
tmp <- ifelse(is.na(tmp)==TRUE, 0, tmp)
# remove the label column from the table
##plot it
pheatmap(tmp , color = colorRampPalette(brewer.pal(3, "YlGnBu"))(10),
         display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
         cluster_rows = FALSE, cluster_cols = FALSE, cellwidth=30, cellheight = 30,
         main = paste0(names(files_f1[i])), 
         width = 6, 
         height = 6,
         number_format="%.2f", filename=paste0("results/plots/plot_f1_",names(files_f1[i]),".pdf")
)
         
}

# Appendix