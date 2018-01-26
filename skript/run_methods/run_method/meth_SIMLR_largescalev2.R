#####################
# SIMLR large scale
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 

analyze_simlrlargescale <- function( datatype, dataset ){ 
  
  source("skript/helper_files/Helper_functions.R")
  # source method SIMLR large scale
  source("skript/run_methods/run_functions/run_function_SIMLR_largescale.R")
  
    # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  

  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters default, filtered, 
  par.c1 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4,
    simDataKumar2 = 4
    
  )
  # unfiltered
  par.c2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 10,
    simDataKumar = 4,
    simDataKumar2 = 4
    
  )
 
  # optimalk
  par.c3 <-  list(
    kumar2015 = 3,
    trapnell2014 = 2,
    zhengmix2016=4,
    koh2016 =10,
    simDataKumar=3,
    simDataKumar2=3
    
  )
  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered")| (datatype=="smooth")) { par.c <- par.c1 } 
  else if ((datatype == "unfiltered")) { par.c <- par.c2 }
  
  else if ((datatype == "default")) { par.c <- par.c1 }
  else if ((datatype == "optimalk")) { par.c <- par.c3 }
  else {print("datatype not defined") }
  

  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.c) == names(data) )
  # RUN SIMLR
  
  run_function_simlr(data[paste0(dataset)],labels[paste0(dataset)], par.c[paste0(dataset)], datatype )
  

}
