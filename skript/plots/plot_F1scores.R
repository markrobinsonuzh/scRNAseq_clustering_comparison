##################################################
### plot heatmap for F1 scores in all clusters  ##
##################################################


## load libraries
library(plyr)

library(dplyr)
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

# load dataset
load(files_f1[[1]])

## plot pheatmap
# create table with data 
res.f1$pcaReduce$f1
res.f1$pcaReduce$act

unique(rapply(res.f1, function(x) head(x, 1)))
sapply(res.f1, mean)

lapply(res.f1, function(x) x[1])



x <- lapply(res.f1, function(x)cbind(f1=x$f1, label=x$act)) 
x <- x[-c(7,8)]
table(x)
?table
%>%unlist()%>%as.matrix(nrow=3,ncol=7)
?as.matrix
names(res.f1[res.f1])

cbind(f1=res.f1$pcaReduce$f1,label=res.f1$pcaReduce$act)

unlist(res.f1)
# rearrange data 

##############