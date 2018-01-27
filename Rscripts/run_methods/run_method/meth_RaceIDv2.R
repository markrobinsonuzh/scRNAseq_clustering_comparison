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
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
  # load data sets
  data <- load_data(files, DATA_DIR)
  # load cell labels
  labels <- load_labels(data) 
  # RUN RaceID
  # set the maxexpr parameter for UMI count data. discarding genes with at least maxexpr transcripts in at least a single cell after normalization or downsampling.
  #the cln parameter is the number of cluster k in kmeans. Default is do.gap true and cln to 0, then gap stat is used to find k.
  
  #default:

  do.gap1 <-  list(
    kumar2015 = TRUE,
    trapnell2014 = TRUE,
    zhengmix2016 = TRUE,
    koh2016 = TRUE,
    simDataKumar=TRUE,
    simDataKumar2=TRUE
    
  )
  par.cln1 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0,
    simDataKumar2=0
    
  )
  par.mintotal1 <- list(
    kumar2015 = 3000,
    trapnell2014 = 3000,
    zhengmix2016 = 3000,
    koh2016 = 3000,
    simDataKumar=3000,
    simDataKumar2=3000
    
  )
  par.maxexpr1 <-  list(
    kumar2015 = Inf,
    trapnell2014 = Inf,
    zhengmix2016 =500,
    koh2016 = Inf,
    simDataKumar=Inf,
    simDataKumar2=Inf
  )
  par.minexprs1 <-  list(
    kumar2015 = 5,
    trapnell2014 = 5,
    zhengmix2016 = 5,
    koh2016 = 5,
    simDataKumar=5,
    simDataKumar2=5
  )
  par.minnumber1 <-  list(
    kumar2015 = 1,
    trapnell2014 = 1,
    zhengmix2016 = 1,
    koh2016 = 1,
    simDataKumar=1,
    simDataKumar2=1
  )
  # unfilterd;
  do.gap2 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 = FALSE,
    koh2016 = FALSE,
    simDataKumar=FALSE,
    simDataKumar2=FALSE
    
  )
  par.cln2 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 10,
    simDataKumar=4,
    simDataKumar2=4
    
  )

  par.mintotal2 <- list(
    kumar2015 = 1000,
    trapnell2014 = 500,
    zhengmix2016 = 420,
    koh2016 = 400,
    simDataKumar=1200,
    simDataKumar2=1200
    
  )
  par.maxexpr2 <-  list(
    kumar2015 = Inf,
    trapnell2014 = Inf,
    zhengmix2016 =500,
    koh2016 = Inf,
    simDataKumar=Inf,
    simDataKumar2=Inf
  )
  par.minexprs2 <-  list(
    kumar2015 = 2,
    trapnell2014 = 2,
    zhengmix2016 = 5,
    koh2016 = 2,
    simDataKumar=5,
    simDataKumar2=2
  )
  par.minnumber2 <-  list(
    kumar2015 = 2,
    trapnell2014 = 2,
    zhengmix2016 = 2,
    koh2016 = 2,
    simDataKumar=2,
    simDataKumar2=2
  )
  # filered and  smoothed data
  do.gap3 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 = FALSE,
    koh2016 = FALSE,
    simDataKumar=FALSE,
    simDataKumar2=FALSE
    
  )
  par.cln3 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=4,
    koh2016 = 9,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  
  par.mintotal3 <- list(
    kumar2015 = 1,
    trapnell2014 = 1,
    zhengmix2016 = 1,
    koh2016 = 1,
    simDataKumar=1,
    simDataKumar2=1
    
  )
  par.maxexpr3 <-  list(
    kumar2015 = Inf,
    trapnell2014 = Inf,
    zhengmix2016 =Inf,
    koh2016 = Inf,
    simDataKumar=Inf,
    simDataKumar2=Inf
  )
  par.minexprs3 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0,
    simDataKumar2=0
  )
  par.minnumber3 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0,
    simDataKumar2=0
  )
  # optimalk
  do.gap4 <-  list(
    kumar2015 = FALSE,
    trapnell2014 = FALSE,
    zhengmix2016 = FALSE,
    koh2016 = FALSE,
    simDataKumar=FALSE,
    simDataKumar2=FALSE
    
  )
  par.cln4 <-  list(
    kumar2015 = 3,
    trapnell2014 = 3,
    zhengmix2016=5,
    koh2016 = 13,
    simDataKumar=3,
    simDataKumar2=4
    
  )
  par.mintotal4 <- list(
    kumar2015 = 1,
    trapnell2014 = 1,
    zhengmix2016 = 1,
    koh2016 = 1,
    simDataKumar=1,
    simDataKumar2=1
    
  )
  par.maxexpr4 <-  list(
    kumar2015 = Inf,
    trapnell2014 = Inf,
    zhengmix2016 =Inf,
    koh2016 = Inf,
    simDataKumar=Inf,
    simDataKumar2=Inf
  )
  par.minexprs4 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0,
    simDataKumar2=0
  )
  par.minnumber4 <-  list(
    kumar2015 = 0,
    trapnell2014 = 0,
    zhengmix2016 = 0,
    koh2016 = 0,
    simDataKumar=0,
    simDataKumar2=0
  )
  # which parameter set
  if ((datatype == "default")) { do.gap <- do.gap1;  par.cln <- par.cln1; par.mintotal <- par.mintotal1; par.maxexpr <- par.maxexpr1; par.minexprs <- par.minexprs1;  par.minnumber <-  par.minnumber1  }
  else if ((datatype == "unfiltered" ) ) {do.gap <- do.gap2;  par.cln <- par.cln2; par.mintotal <- par.mintotal2; par.maxexpr <- par.maxexpr2; par.minexprs <- par.minexprs2;  par.minnumber <-  par.minnumber2}
  else if ( (datatype=="filtered") | (datatype=="smooth")) {do.gap <- do.gap3;  par.cln <- par.cln3; par.mintotal <- par.mintotal3; par.maxexpr <- par.maxexpr3; par.minexprs <- par.minexprs3;  par.minnumber <-  par.minnumber3}
  else if ((datatype == "optimalk")) { do.gap <- do.gap4;  par.cln <- par.cln4; par.mintotal <- par.mintotal4; par.maxexpr <- par.maxexpr4; par.minexprs <- par.minexprs4;  par.minnumber <-  par.minnumber4 }
  else {print("datatype not defined!") }
  
  print(do.gap )
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(do.gap ) == names(data) )
  
  # run RACEID
  run_function_raceid( data=data[paste0(dataset)], labels=labels[paste0(dataset)], datatype=datatype ,par.mintotal=par.mintotal[paste0(dataset)], par.maxexpr=par.maxexpr[paste0(dataset)],
                       par.minnumber= par.minnumber[paste0(dataset)], par.minexprs= par.minexprs[paste0(dataset)],
                       do.gap=do.gap[paste0(dataset)],cln=par.cln[paste0(dataset)]) 
}
### Appendix



