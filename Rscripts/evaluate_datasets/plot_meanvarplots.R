#######################################
# Mean variance plots for all data sets
#######################################

# load libraries
source("skript/helper_files/Helper_functions.R")

library(vsn)
library(cowplot)
library(DESeq2)
library(gridExtra)
library(DESeq)
library(cowplot)
library(scater)



# file paths, defined in FILES.R

source("FILESraw.R")

# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# rlog

#rld <- rlog((round(count_lstpm,0)), blind=TRUE) # takes too much time

###### define funcion  meansd_plot for transformation and plotting of count data
meansd_plot <- function(data, ranks,cofactor ){
# compute transformations
count_lstpm <- ( as.matrix(assay(data, "counts")) ) # extract count data
# log and asin transformaations
count_lstpm.log <- log2(count_lstpm +1)  # log2
cofactor <- cofactor
count_lstpm.asinh <- asinh(sqrt(count_lstpm/cofactor) )# arc sin
# vst transform
countdata <- newCountDataSet((round(count_lstpm,0)), conditions = array(1,dim = ncol(count_lstpm)))
countdata <- estimateSizeFactors( countdata)
cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
vstdata <- varianceStabilizingTransformation( cdsBlind )
count_lstpm.vst <- exprs(vstdata)

# plot it
msd.non <- meanSdPlot(count_lstpm, plot = FALSE, ranks = ranks)
msd.log <- meanSdPlot(count_lstpm.log, plot = FALSE, ranks = ranks)
msd.asin <- meanSdPlot(count_lstpm.asinh, plot = FALSE, ranks = ranks)
msd.vst <- meanSdPlot(count_lstpm.vst , plot = FALSE, ranks = ranks)
# in grid

grid.arrange( msd.non$gg+ggtitle("untransformed"), msd.log$gg+ggtitle("log transformed"), msd.asin$gg+ggtitle(paste0("arcussin cofactor=", cofactor)), msd.vst$gg+ggtitle("vst"),ncol=2)

}

## plot all data sets
for (i in seq_len(length(data))){
  
  FILE_NAME<- paste0("results/QC_data/meanvarplots_", names(data)[i],".pdf")
  pdf(FILE_NAME)
  meansd_plot(data=data[[i]],ranks= FALSE,cofactor=1)
  dev.off()
}

# Appendix
###### define vst function,adapted from DEseq2, from DEseq vignette, vignette("DESeq2")

