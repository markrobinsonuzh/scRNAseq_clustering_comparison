#####################
# ZINB-WaVE
#####################
# uses reduced raw count  matrix. 
source("skript/helper_files/Helper_functions.R")
# source file paths: fileterd , raw etc.

source("FILES.R")
# source method zinbwave
source("skript/run_methods/run_functions/run_function_zinbwave.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# parameters: number of clusters k in kmeans, number of genes n.genes for ZINBWaVE
par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016 = 10,
  simDataKumar = 3
   
)

n.genes <- list(
  kumar2015 = 1000,
  trapnell2014 = 1000,
  zhengmix2016 = 200,
  koh2016 = 1000,
  simDataKumar = 1000
  
)
# define datatype

datatype <- "default"
# check if files, parameters and data are the same:
names(files)==names(data) 


# RUN ZINB-WaVE
run_function_zinbwave( data, labels, par.k, n.genes,datatype ) 


