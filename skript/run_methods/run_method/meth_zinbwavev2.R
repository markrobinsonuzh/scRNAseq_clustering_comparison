#####################
# ZINB-WaVE
#####################
analyze_zinbwave <- function( datatype, dataset ) {
  # uses reduced raw count  matrix. 
  source("skript/helper_files/Helper_functions.R")
  # source method zinbwave
  source("skript/run_methods/run_functions/run_function_zinbwave.R")
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
  
  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters: number of clusters k in kmeans, number of genes n.genes for ZINBWaVE
  #for all runmodes
  n.genes <- list(
    kumar2015 = 1000,
    trapnell2014 = 1000,
    zhengmix2016 = 200,
    koh2016 = 1000,
    simDataKumar = 1000,
    simDataKumar2 = 1000
    
    
  )
  # default:

  
  par.k1 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4,
    simDataKumar2 = 4
    
    
  )
  # filtered


  
  par.k2 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4,
    simDataKumar2 = 4
    
    
  )
  #unfiltered

  
  par.k3 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 10,
    simDataKumar = 4,
    simDataKumar2 = 4
    
    
  )
  # opt k
  
  par.k4 <- list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016 = 4,
    koh2016 = 9,
    simDataKumar = 4,
    simDataKumar2 = 4
    
    
  )
  # which parameter set
  if ( (datatype=="filtered") | (datatype=="smooth")) { n.genes  <-n.genes  ; par.k <-  par.k2}
  else if ((datatype == "unfiltered")) { n.genes  <- n.genes  ; par.k <-  par.k3 }
  else if ((datatype == "default")) { n.genes  <- n.genes ; par.k <-  par.k1 }
  else if ((datatype == "optimalk")) { n.genes  <- n.genes ; par.k <-  par.k4 }
  else {print("datatype not defined") }
  
  print( n.genes  )
  print( par.k )
  
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.k) == names(data) )
  stopifnot( names(n.genes) == names(data) )
  
  # RUN ZINB-WaVE
  run_function_zinbwave( data[paste0(dataset)], labels[paste0(dataset)], par.k[paste0(dataset)],n.genes[paste0(dataset)], datatype ) 

}

