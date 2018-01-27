########################
# CIDR
########################

analyze_cidr <- function(dataype, dataset) {
  # source helper files
  source("skript/helper_files/Helper_functions.R")
  # source method CIDR
  source("skript/run_methods/run_functions/run_function_cidr.R")
  
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
  
  # load data sets
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters
  # define number of clusters. default is NULL and cidr uses calinski index to find k.
  # default mode:
  par.k1 <-  list(
    kumar2015 = NULL,
    trapnell2014 = NULL,
    zhengmix2016=NULL,
    koh2016 = NULL,
    simDataKumar=NULL,
    simDataKumar2=NULL
    
  )
  par.nPC1 <-  list(
    kumar2015 = 4,
    trapnell2014 = 4,
    zhengmix2016=4,
    koh2016 = 4,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  # unfiltered
  par.k2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 10,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  par.nPC2 <-  list(
    kumar2015 = 5,
    trapnell2014 = 10,
    zhengmix2016=8,
    koh2016 = 8,
    simDataKumar=3,
    simDataKumar2=3
    
  )
  
  # filtered 
  par.k2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 9,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  par.nPC2 <-  list(
    kumar2015 = 5,
    trapnell2014 = 10,
    zhengmix2016=8,
    koh2016 = 8,
    simDataKumar=3,
    simDataKumar2=3
    
  )
  # optimalk
  par.k3 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=5,
    koh2016 = 13,
    simDataKumar=4,
    simDataKumar2=4
    
  )

  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered") | (datatype=="smooth")) { par.k <- par.k2 ; par.nPC <- par.nPC2}
  else if ((datatype == "default")) { par.k <- par.k1 ; par.nPC <- par.nPC1 }
  else if ((datatype == "optimalk")) { par.k <- par.k3 ; par.nPC <- par.nPC2 }
  else {print("datatype not defined") }
  
  print(par.k)
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.k) == names(data) )
  # RUN cidr
  run_function_cidr( data[paste0(dataset)], labels[paste0(dataset)], par.k[paste0(dataset)],par.nPC2[paste0(dataset)], datatype )
}



### Appendix