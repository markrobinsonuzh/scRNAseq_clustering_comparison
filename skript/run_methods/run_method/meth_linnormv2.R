
#####################
# Linnorm
#####################
# to do : include spikeinns



analyze_linnorm <- function(dataype, dataset) {
  # source method linnomr
  source("skript/run_methods/run_functions/run_function_linnorm.R")
  # source helper
  source("skript/helper_files/Helper_functions.R")
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")
  } else {
    if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
    if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
    
    else {print("datatype not defined") }
  }
  
  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
  # num_center is number of k in kmeans, default is k = c(1:20).
  # Set a lower cutoff for the zhengmix data
  par.minNonZeroPortion <- list(
    kumar2015 = 0.75,
    trapnell2014 =  0.75,
    zhengmix2016 =  0.1,
    koh2016 =  0.75,
    simDataKumar=0.75
  )

  # default:
  par.num_center1 <- list(
    kumar2015 = c(1:20),
    trapnell2014 =  c(1:20),
    zhengmix2016 = c(1:20),
    koh2016 = c(1:20),
    simDataKumar=c(1:20)
  )
  par.BE_strength1 <- list(
    kumar2015 = 0.5,
    trapnell2014 =  0.5,
    zhengmix2016 = 0.5,
    koh2016 = 0.5,
    simDataKumar=0.5
  )

  # filtered , unfiltered:
  par.num_center2 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 =9,
    simDataKumar=4
  )
  par.BE_strength2 <- list(
    kumar2015 = 0.75,
    trapnell2014 =  0.75,
    zhengmix2016 = 0.75,
    koh2016 = 0.75,
    simDataKumar=0.1
  )
  # optimalk
  par.num_center3 <-  list(
    kumar2015 = 3,
    trapnell2014 = 2,
    zhengmix2016=5,
    koh2016 = 11,
    simDataKumar=3
  )
  
  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered")| (datatype=="smooth")) {par.num_center <- par.num_center2; par.BE_strength <- par.BE_strength2
  } else {
    if ((datatype == "default")) { par.num_center <- par.num_center1 ; par.BE_strength <- par.BE_strength1}
    if ((datatype == "optimalk")) { par.num_center <- par.num_center3 ; par.BE_strength <- par.BE_strength2 }
    else {print("datatype not defined") }
  }
  print(par.num_center)
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.num_center) == names(data) )

  
  # RUN cidr
  run_function_linnorm(data[paste0(dataset)], labels[paste0(dataset)],  par.minNonZeroPortion[paste0(dataset)], par.num_center[paste0(dataset)],par.BE_strength[paste0(dataset)], datatype)
  
}
