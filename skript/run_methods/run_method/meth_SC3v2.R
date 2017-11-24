########################################
# SC3
#########################################



analyze_sc3 <- function( datatype ){ 
  
  source("skript/helper_files/Helper_functions.R")
# source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")) { source("FILES.R"); print("filtered files")
  } else {
    if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
    else {print("datatype not defined") }
  }
  # source method SC3
  source("skript/run_methods/run_functions/run_function_sc3.R")
  # load data sets

  data <- load_data(files, DATA_DIR)

  # load cell labels
  labels <- load_labels(data) 
  # run the analysis
  # default
par.ks1 <- list(
  kumar2015 = 2:10,
  trapnell2014=2:10,
  zhengmix2016 =2:10,
  koh2016 = 2:10,
  simDataKumar=2:10
)
par.k_estimator1 <- list(
  kumar2015 = TRUE,
  trapnell2014=TRUE,
  zhengmix2016 =TRUE,
  koh2016 =TRUE,
  simDataKumar=TRUE
)
par.k1 <- list(
  kumar2015 = 3,
  trapnell2014=3,
  zhengmix2016 =4,
  koh2016 = 10,
  simDataKumar=4
)
# filtered , unfilterd
par.ks2 <- list(
  kumar2015 = 2:6,
  trapnell2014=2:6,
  zhengmix2016 =2:6,
  koh2016 = 2:12,
  simDataKumar=2:6
)
par.k_estimator2 <- list(
  kumar2015 = FALSE,
  trapnell2014=FALSE,
  zhengmix2016 =FALSE,
  koh2016 =FALSE,
  simDataKumar=FALSE
)
par.k2 <- list(
  kumar2015 = 3,
  trapnell2014=3,
  zhengmix2016 =4,
  koh2016 = 10,
  simDataKumar=4
)

# which parameter set
if ((datatype == "unfiltered" ) | (datatype=="filtered")) { par.ks  <- par.ks2  ;par.k_estimator <- par.k_estimator2; par.k <- par.k2
} else {
  if ((datatype == "default")) { par.ks  <- par.ks1 ;par.k_estimator <- par.k_estimator1 ; par.k <- par.k1 }
  else {print("datatype not defined") }
}
print(par.ks )
print(par.k_estimator)
print(par.k)

# check data files and parameters
stopifnot( names(files) == names(data) )
stopifnot( names(par.ks) == names(data) )
# RUN SC3
run_function_sc3( data, labels, par.ks, par.k_estimator ,par.k,datatype )
}

### Appendix


