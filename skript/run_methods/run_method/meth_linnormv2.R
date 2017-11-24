
#####################
# Linnorm
#####################
# to do : include spikeinns



analyze_linnorm <- function(dataype) {
  # source method linnomr
  source("skript/run_methods/run_functions/run_function_linnorm.R")
  # source helper
  source("skript/helper_files/Helper_functions.R")
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
  
  # define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
  # num_center is number of k in kmeans, default is k = c(1:20).
  # Set a lower cutoff for the zhengmix data
  par.minNonZeroPortion <- list(
    kumar2015 = 0.75,
    trapnell2014 =  0.75,
    zhengmix2016 =  0.25,
    koh2016 =  0.75,
    simDataKumar=0.75
  )
  # default:
  par.num_center1 <- list(
    kumar2015 = 1:5,
    trapnell2014 =  1:5,
    zhengmix2016 = 1:5,
    koh2016 = 1:12,
    simDataKumar=1:5
  )
  # filtered , unfiltered:
  par.num_center2 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 10,
    simDataKumar=4
  )
  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered")) {par.num_center <- par.num_center2
  } else {
    if ((datatype == "default")) { par.num_center <- par.num_center1 }
    else {print("datatype not defined") }
  }
  print(par.num_center)
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.num_center) == names(data) )
  
  
  # RUN cidr
  run_function_linnorm(data, labels, datatype , par.minNonZeroPortion, par.num_center)
  
  
}
