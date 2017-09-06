#######################################
# Mean variance plots for all data sets
#######################################

# load libraries
source("skript/helper_files/Helper_functions.R")

library(vsn)
library(cowplot)


# file paths

DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rds")
)

# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
########## which dataset?
data <- data[[3]]


########### try different transformations
count_lstpm <- as.matrix(get_exprs(data, "counts")) # no trnasformation
count_lstpm.log <- log2(count_lstpm +1)  # log2
count_lstpm.norm <- as.matrix(exprs(data)) # standart from scater, is log2
########## mean var plots
########### plot
pdf("results/QC_data/meanvarplots_xue2013.pdf")

par(mfrow=c(1,3))
meansdplot(data=count_lstpm, title = "count_lstpm" ,ylim=c(0,2000), rank=TRUE)
meansdplot(count_lstpm.log, title="count_lstpm.log2", rank=TRUE )
meansdplot(count_lstpm.norm, title="count_lstpm.expr", rank=TRUE )

dev.off()
