#####################
# DBscan
#####################
# dbscan is a density based clustering method. The k neirest neighborhoud parameter and neighborhood size (epsilon) has to be defined.
# This is usually done using a k neirest neighbour distance plot. for this strategy the k-nearest neighbor distance is plotted (k-NN distance)
# this is the distance from a choosen point to its k nearest neighbors. As a rule of thumb k >= d+1, here d is the dimension of the data.
# A sharp change in the distance plot indicates the epsilon distance. 
# If only one single cluster is obtained, then often epsilon is too large or MinPts is too small.
# similary, if epsilon is too small or MinPts is too large then every point becomes a noise point.
# note that dbscan has problems with high dim data , so we should reduce the dimension, which is done using PCA and working on latent space with 50 dim.

# load helper files
source("skript/helper_files/Helper_functions.R")
# source file paths: fileterd , raw etc.
source("FILES.R")
# source method dbscan 
source("skript/run_methods/run_functions/run_function_dbscan.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 
# parameters k is nearest neighbor, epsilon is the the neighboorhoodsize and Pts the minim neighboors
# par. k as ten percent of dataset

par.k <- list(
  kumar2015 = 25,
  trapnell2014 = 22,
  zhengmix2016 = 196,
  koh2016 = 53,
  simDataKumar=50
)
par.eps <- list(
  kumar2015 = 40000,
  trapnell2014 = 60000,
  zhengmix2016 = 70,
  koh2016 = 50000,
  simDataKumar=40000
)
par.Pts <- list(
  kumar2015 = 5,
  trapnell2014 = 5,
  zhengmix2016 = 5,
  koh2016 = 5,
  simDataKumar=5
)
# define datatype

datatype <- "default"

# check if files, parameters and data are the same:

names(files)==names(par.Pts) 
# RUN dbscan

run_function_dbscan( data, labels, par.k, par.eps, par.Pts , datatype )
