#####################
# RACEID
#####################
# RaceID is an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. 
# The method is based on transcript counts obtained with unique molecular identifies.
analyze_raceid <- function(dataype, dataset) {
  
  source("skript/helper_files/Helper_functions.R")
  # source method raceid
  source("skript/run_methods/run_functions/run_function_raceid.R")
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
  # RUN RaceID
  # set the maxexpr parameter for UMI count data. discarding genes with at least maxexpr transcripts in at least a single cell after normalization or downsampling.
  #the cln parameter is the number of cluster k in kmeans. Default is do.gap true and cln to 0, then gap stat is used to find k.
  par.mintotal <- list(
    kumar2015 = 3000,
    trapnell2014 = 3000,
    zhengmix2016 = 200,
    koh2016 = 3000,
    simDataKumar=3000
  )
  par.maxexpr <-  list(
    kumar2015 = Inf,
    trapnell2014 = Inf,
    zhengmix2016 = Inf,
    koh2016 = Inf,
    simDataKumar=Inf
  )
  
  
  #default:

  do.gap1 <-  list(
    kumar2015 = TRUE,
    trapnell2014 = TRUE,
    zhengmix2016 = TRUE,
    koh2016 = TRUE,
    simDataKumar=TRUE
  )
  par.cln1 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0
  )
  #filtered, unfilterd;
  do.gap2 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 = FALSE,
    koh2016 = FALSE,
    simDataKumar=FALSE
  )
  par.cln2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 9,
    simDataKumar=4
  )
  # which parameter set
  if ((datatype == "unfiltered" ) | (datatype=="filtered")) {do.gap <- do.gap2;  par.cln <- par.cln2
  } else {
    if ((datatype == "default")) { do.gap <- do.gap1;  par.cln <- par.cln1}
    else {print("datatype not defined") }
  }
  print(do.gap )
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(do.gap ) == names(data) )
  
  # run RACEID
  run_function_raceid( data=data[paste0(dataset)], labels=labels[paste0(dataset)], datatype=datatype ,par.mintotal=par.mintotal[paste0(dataset)], par.maxexpr=par.maxexpr[paste0(dataset)],do.gap=do.gap[paste0(dataset)],cln=par.cln[paste0(dataset)]) 
}

### Appendix



