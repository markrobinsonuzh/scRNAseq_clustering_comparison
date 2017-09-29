
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
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_SRP073808.rda")
  
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

name <- "xue2013"


######## function
meansd_plot <- function(data, name ){
########### try different transformations
  data <- data
count_lstpm <- as.matrix(get_exprs(data, "counts")) # no trnasformation
count_lstpm.log <- log2(count_lstpm +1)  # log2
count_lstpm.asinh <- asinh(count_lstpm)
#count_lstpm.norm <- as.matrix(exprs(data)) # standart from scater, is log2
########## mean var plots
########### plot
x <- paste0("results/QC_data/meanvarplots_", name,".pdf")
pdf(x)

par(mfrow=c(2,3))
meansdplot(data=count_lstpm, title = "count_lstpm" ,ylim=c(0,2000), rank=TRUE)
meansdplot(count_lstpm.log, title="count_lstpm.log2", rank=TRUE )
meansdplot(data=count_lstpm.asinh, title = "count_lstpm.asinh" ,ylim=c(0,5), rank=TRUE)

dev.off()
}
meansd_plot(data=data[[4]],name="koh2016")

