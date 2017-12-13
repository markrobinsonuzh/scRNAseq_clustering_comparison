#####################
# SIMLR normal version
#####################
# Given a gene expression matrix ( normalized , log transformed ) SIMLR learns a distance metric through multiple kernels (Gaussian) that best fits the data. 
# These similiraties are then used for a RtSNE step to represent the data in lower dimension space and clustering using kmeans. 
# Parameterrs given by the user are the number of cluster c to be estimated over the expression matrix. 

analyze_simlr <- function( datatype, dataset ){ 
  
  set.seed(1234)
  
  source("skript/helper_files/Helper_functions.R")
  # source method SIMLR large scale
  source("skript/run_methods/run_functions/run_function_SIMLR.R")
  
  
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")
  } else {
    if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
    else {print("datatype not defined") }
  }
  
  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters
  #default
  par.c1 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4
  )
  normalize1 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 =FALSE,
    koh2016 = FALSE,
    simDataKumar = FALSE
  )
  # parameters
  #filtered unfiltered
  par.c2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4
  )
  normalize2 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 = FALSE,
    koh2016 = FALSE,
    simDataKumar = FALSE
  )
  # optimalk
  par.c3 <-  list(
    kumar2015 = 3,
    trapnell2014 = 2,
    zhengmix2016=4,
    koh2016 =10,
    simDataKumar=3
  )
  normalize3 <-  list(
    kumar2015 = TRUE,
    trapnell2014 = TRUE,
    zhengmix2016 =TRUE,
    koh2016 = TRUE,
    simDataKumar = TRUE
  )
  
  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered")) { par.c <- par.c2 ; normalize <-  normalize2
  } else {
    if ((datatype == "default")) {  par.c <- par.c1 ; normalize <-  normalize1 }
    if ((datatype == "optimalk")) { par.c <- par.c3 ; normalize <-  normalize3 }
    else {print("datatype not defined") }
  }
  print(par.c)
  print(normalize)
  
  
  
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.c) == names(data) )
  # RUN SIMLR
  
  run_function_simlrnormal(data[paste0(dataset)],labels[paste0(dataset)], par.c[paste0(dataset)], normalize[paste0(dataset)],datatype )
  
  
}
