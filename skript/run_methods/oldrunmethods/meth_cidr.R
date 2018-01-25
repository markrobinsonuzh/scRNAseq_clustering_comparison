########################
# CIDR
########################


  # source helper files
  source("skript/helper_files/Helper_functions.R")
  # source method CIDR
  source("skript/run_methods/run_functions/run_function_cidr.R")
  
  # source file paths: fileterd , raw etc.
  source("FILES.R")

  # load data sets
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters
  # define number of clusters. default is NULL and cidr uses calinski index to find k.
  # default mode:
  par.k <-  list(
    kumar2015 = NULL,
    trapnell2014 = NULL,
    zhengmix2016=NULL,
    koh2016 = NULL,
    simDataKumar=NULL
  )
  # filtered , unfiltered
  par.k <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 10,
    simDataKumar=3
  )

  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.k) == names(data) )
  # default
  datatype <-  "default"
  
  # RUN cidr
  run_function_cidr( data, labels, par.k, datatype )





### Appendix