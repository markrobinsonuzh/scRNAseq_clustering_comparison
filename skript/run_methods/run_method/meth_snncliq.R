#####################
# SNN Cliq
#####################



analyze_snncliq <- function( datatype, dataset){ 
  # source method rtsnekmeans and helpers
  source("skript/run_methods/run_functions/run_function_snncliq.R") 
  source("skript/helper_files/Helper_functions.R")
  
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
  
  # load data sets
  data <- load_data(files, DATA_DIR)
  # load cell labels
  labels <- load_labels(data) 
  # default mode
  par.k <- list(
    kumar2015 =ncol(data[["kumar2015"]])* c(0.1),
    trapnell2014 = ncol(data[["trapnell2014"]])*c(0.1),
    zhengmix2016 = ncol(data[["zhengmix2016"]])*c(0.1),
    koh2016 = ncol(data[["koh2016" ]])*c(0.1),
    simDataKumar=ncol(data[["simDataKumar"]])*c(0.1),
    simDataKumar2=ncol(data[["simDataKumar2"]])*c(0.1)
    
  )
  par.m <- list(
    kumar2015 = 0.5,
    trapnell2014 =  0.5,
    zhengmix2016 =  0.5,
    koh2016=  0.5,
    simDataKumar= 0.5,
    simDataKumar2=0.5
  )
  par.r <- list(
    kumar2015 = 0.7,
    trapnell2014 = 0.7,
    zhengmix2016 = 0.7,
    koh2016= 0.7,
    simDataKumar=0.7,
    simDataKumar2=0.7
  )

  print(par.k)

  
  # check data files and parameters
  stopifnot( names(files) == names(data) )

  # Run tSNE and kmeans
  run_function_snncliq( data[paste0(dataset)],labels[paste0(dataset)],par.k[paste0(dataset)],par.m[paste0(dataset)],par.r[paste0(dataset)],distan="euclidean", datatype )
  
}
