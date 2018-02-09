#####################
# SIMLR
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 
# 
set.seed(1234)

source("skript/helper_files/Helper_functions.R")


# source file paths: fileterd , raw etc.

source("FILES.R")
# source method CIDR
source("skript/run_methods/run_functions/run_function_simlr.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# Set paramaeters: c is the number of cluster expected. we define the tuning parameter k to 10 as it is the standart parameter, the number of dimension are NA
par.c <-  list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016 = 10,
  simDataKumar = 3
)
# define datatype

datatype <- "filtered"
# check if files, parameters and data are the same:
names(data) == names(par.c)
# RUN cidr
run_function_simlr( data=data, labels=labels, par.c=par.c, datatype=datatype, normalize = TRUE  ) 
