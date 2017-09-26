###################################
### plot confusion matrix for method RtSNE 
###################################

pdf("results/plots/confusion_matrix_dbscan.pdf")

# load libraries
source("skript/helper_files/Helper_functions.R")

library(reshape)
library(caret)
library(dplyr)
library(magrittr)

###########################
# load labels and cluster #
###########################

# define method:

DATA_DIR <-  "results/dbscan/"
METHOD <- c("dbscan")

### files with the cell labels, "ground truth":
files_labels <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_labels_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_labels_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_labels_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_labels_koh2016.txt"))
)
# read in labels
labels <- read.labels(files_labels = files_labels)
# read in cluster results 
files_clusters <- list(
  kumar2015 = file.path(DATA_DIR, paste0(METHOD,"_clus_kumar2015.txt")),
  trapnell2014 = file.path(DATA_DIR, paste0(METHOD,"_clus_trapnell2014.txt")),
  xue2013 = file.path(DATA_DIR, paste0(METHOD,"_clus_xue2013.txt")),
  koh2016 = file.path(DATA_DIR, paste0(METHOD,"_clus_koh2016.txt"))
)
cluster <- read.cluster(files_clusters = files_clusters)




#####################
# confusion matrix #
#####################
#list <-  vector("list", 3) %>% names() %$% names(files_clusters)
conv.tbl <- vector("list", 3)
for (i in seq_len(length(cluster))){
  conv.tbl[[i]] <- table(cluster=cluster[[i]], label=labels[[i]])
}

##################################
# confusion matrix as ggplot ###
##################################

p <- vector(mode="list", length=length(conv.tbl))
names(conv.tbl) <- names(files_labels)


pheat.list <- vector(mode="list", length=length(conv.tbl))

for (i in seq_len(length(conv.tbl))){
  
  dd <- melt(conv.tbl[[i]])
  
  pheat.plot[[i]] <- pheatmap(conv.tbl[[i]], 
                              color = colorRampPalette(brewer.pal(7, "GnBu"))(100), 
                              display_numbers = TRUE, number_color = "black", fontsize_number = 9, 
                              cluster_rows = FALSE, cluster_cols = FALSE, 
                              main = paste0(names(conv.tbl[i])), 
                              width = 6, 
                              height = 6,
                              number_format="%.0f",
                              filename=paste0("results/plots/confusion_matrix_dbscan",names(conv.tbl[i]),".pdf")
  )
}



dev.off()
### Appendix
