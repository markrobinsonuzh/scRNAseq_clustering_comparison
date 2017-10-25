
#######################################
# Mean variance plots for all data sets
#######################################

# load libraries
source("skript/helper_files/Helper_functions.R")

library(vsn)
library(cowplot)
library(DESeq2)


# file paths, defined in FILES.R

source("FILES.R")


# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

########## which dataset?
#data <- data[[3]]

#name <- "xue2013"
# vst function adapted from https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/
vst <- function(countdata){
  require(DESeq)
  countdata <- newCountDataSet((round(count_lstpm,0)), conditions = array(1,dim = ncol(count_lstpm)))
  countdata <- estimateSizeFactors( countdata)
  cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
  vstdata <- varianceStabilizingTransformation( cdsBlind )
  return(exprs(vstdata))
}

######## define function to transform and plot
meansd_plot <- function(data, name, rank ){
########### try different transformations
data <- data[[1]]

count_lstpm <- ( as.matrix(get_exprs(data, "counts")) ) # no trnasformation
count_lstpm.log <- log2(count_lstpm +1)  # log2
count_lstpm.asinh <- asinh(count_lstpm)
count.integer <- (round(count_lstpm,0))

count_lstpm.vst <- vst( count_lstpm) ## for counts, as it is tpm on count scale necessary to change to integer
#count_lstpm.norm <- as.matrix(exprs(data)) # standart from scater, is log2

#### have a look on the data
count_lstpm[1:3,1:3]
count_lstpm.log[1:3,1:3]
count_lstpm.asinh[1:3,1:3]
count_lstpm.vst[1:3,1:3]
########## mean var plots
########### plot
FILE_NAME<- paste0("results/QC_data/meanvarplots_", name,".pdf")
pdf(FILE_NAME)

par(mfrow=c(2,4))
meansdplot(data=count_lstpm, title = "count_lstpm" ,ylim=c(0,2000), rank=TRUE)
meansdplot(count_lstpm.log, title="count_lstpm.log2", rank=TRUE )
meansdplot(data=count_lstpm.asinh, title = "count_lstpm.asinh" ,ylim=c(0,5), rank=TRUE)
meansdplot(data=count_lstpm.vst, title = "count_lstpm.vst" ,ylim=c(0,5), rank=TRUE)

dev.off()

}

#### use the function, plot with and without rank
for (i in names(data)){
meansd_plot(data=data[[i]],rank = TRUE,name= paste0("rank",names(data)[i]))
}
## try http:// if https:// URLs are not suppo
for (i in names(data)){
meansd_plot(data=data[[i]],rank = FALSE,name= paste0("mean",names(data)[i]))
}

## varaince stabilizing transformation
