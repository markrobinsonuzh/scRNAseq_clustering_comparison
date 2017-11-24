#####################
# SIMLR large scale
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 

analyze_simlrlargescale <- function( datatype ){ 
  
  set.seed(1234)
  
  source("skript/helper_files/Helper_functions.R")
  # source method SIMLR large scale
  source("skript/run_methods/run_functions/run_function_SIMLR_largescale.R")
  
  
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")) { source("FILES.R"); print("filtered files")
  } else {
    if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
    else {print("datatype not defined") }
  }

  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters
  par.c <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 10,
    simDataKumar = 4
  )
 
  print(par.c)
  
  
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.c) == names(data) )
  # RUN SIMLR
  
  run_function_simlr(data,labels, par.c, datatype )
  

}
