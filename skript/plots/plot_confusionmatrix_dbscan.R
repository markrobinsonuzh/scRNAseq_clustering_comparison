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

for (i in seq_len(length(conv.tbl))){
  
  dd <- melt(conv.tbl[[i]])
  print(ggplot(dd, aes(as.factor(dd[,1]), dd[,2], fill=value))+geom_tile()
        +geom_text(aes(label=dd$value),colour="white")
        +ggtitle(paste0(names(conv.tbl[[i]])))+theme(axis.text.x=element_text(size = 15),axis.text.y=element_text(size = 15))
        +ylab("labels")+xlab("cluster")
        
  )
}



dev.off()
### Appendix
