###########################
# Stability analysis
###########################

### load libraries
library(cluster)
library(dplyr)
library(fpc)
library(scater)
library(clValid)
library(kohonen)
library(mclust)

source("skript/helper_files/Helper_functions.R")
source("skript/stability_analysis/interface_clusterboot.R")
# source file paths: fileterd , raw etc.
source("FILES.R")


# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 
# one daataset
data <- assay(data[[1]], "normcounts")


##############################################££££££££££££££££££££££
#  Clusterwise cluster stability assessment by resampling 
#####################################################################
par.perp <- 30
### run clusterboot for Rtsnekmeans

plot(km.boot)
tsne.boot$partition
tsne.boot$bootresult
tsne.boot$bootmean
tsne.boot$result
tsne.boot$multipleboot
### run clusterboot for   
tsne.boot <- clusterboot(data=t(data), B=50, bootmethod="boot",
                         clustermethod=rtsnekmeansCBI,dissolution=0.5,
                         recover=0.75,
                         k=3, perplexity=par.perp, seed=NULL, count = TRUE )


options(digits=3)
set.seed(20000)
face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
cf1 <- clusterboot(face,B=3,bootmethod=
                     c("boot","noise","jitter"),clustermethod=kmeansCBI,
                   krange=5,seed=15555)

