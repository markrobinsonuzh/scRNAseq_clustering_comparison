##################################################
### plot heatmap for F1 scores in all clusters  ##
##################################################


## load libraries
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

## define the data directories
DATA_DIR <-  "results/run_results"


## read in cluster results from Rdata filess

# files directories per dataset
files_f1 <- list(
  kumar2015 = file.path(DATA_DIR, "f1_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "f1_trapnell2014.rda"),
  xue2013 = file.path(DATA_DIR, "f1_xue2013.rda"),
  koh2016 = file.path(DATA_DIR, "f1_koh2016.rda")
)

for ( i in names(files_f1)){
# load dataset
tmp <- lapply(files_f1[[i]], function(x) get(load(x)))
## plot pheatmap
# create table with data 
# rearrange data 
tmp <- ldply(tmp[[1]], data.frame) %>% select(.id, f1, act)
tmp <- daply(tmp, .(act, .id), function(x) x$f1)

##plot it
pheatmap(tmp , color = colorRampPalette(brewer.pal(3, "YlGnBu"))(10),
         display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         main = "f1 plot", 
         width = 6, 
         height = 6,
         number_format="%.2f",
         filename=paste0("results/plots/plot_f1_",names(files_f1[i]),".pdf")
         )



}


# Appendix
