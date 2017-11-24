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

# source method CIDR
source("skript/run_methods/run_functions/run_function_rtsnekmeans.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 
# one daataset
data <- assay(data[[1]], "normcounts")


##############################################££££££££££££££££££££££
#  Clusterwise cluster stability assessment by resampling 
#####################################################################

### run clusterboot for Rtsnekmeans
tsne.boot <- clusterboot(data=data, B=50, bootmethod="boot",
                        clustermethod=rtsnekmeansCBI,
                        k=3,perplexity=par.perp, seed=NULL)

plot(km.boot)

### run clusterboot for   
