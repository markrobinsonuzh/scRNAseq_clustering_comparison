########################################
# SC3
#########################################



analyze_sc3 <- function( datatype , dataset){ 
  
  source("skript/helper_files/Helper_functions.R")
# source file paths: fileterd , raw etc.
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
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
  koh2016 = 2:15,
  simDataKumar=2:10,
  simDataKumar2=2:10
  
)
par.k_estimator1 <- list(
  kumar2015 = TRUE,
  trapnell2014=TRUE,
  zhengmix2016 =TRUE,
  koh2016 =TRUE,
  simDataKumar=TRUE,
  simDataKumar2=TRUE
  
)
par.k1 <- list(
  kumar2015 = NULL,
  trapnell2014=NULL,
  zhengmix2016 =NULL,
  koh2016 = NULL,
  simDataKumar=NULL,
  simDataKumar2=NULL
  
)
pct_dropout_max1<- list(
  kumar2015 = 90,
  trapnell2014=90,
  zhengmix2016 =90,
  koh2016 = 90,
  simDataKumar=90,
  simDataKumar2=90
  
)
# unfiltered
par.ks2 <- list(
  kumar2015 = NULL,
  trapnell2014= NULL,
  zhengmix2016 = NULL,
  koh2016 =  NULL,
  simDataKumar= NULL,
  simDataKumar2= NULL
  
)
par.k_estimator2 <- list(
  kumar2015 = FALSE,
  trapnell2014=FALSE,
  zhengmix2016 =FALSE,
  koh2016 =FALSE,
  simDataKumar=FALSE,
  simDataKumar2=FALSE
  
)
par.k2 <- list(
  kumar2015 = 3,
  trapnell2014=3,
  zhengmix2016 =4,
  koh2016 = 10,
  simDataKumar=4,
  simDataKumar2=4
  
)
# filtered 
par.ks2 <- list(
  kumar2015 = NULL,
  trapnell2014= NULL,
  zhengmix2016 = NULL,
  koh2016 =  NULL,
  simDataKumar= NULL,
  simDataKumar2= NULL
  
)
par.k_estimator2 <- list(
  kumar2015 = FALSE,
  trapnell2014=FALSE,
  zhengmix2016 =FALSE,
  koh2016 =FALSE,
  simDataKumar=FALSE,
  simDataKumar2=FALSE
  
)
par.k3 <- list(
  kumar2015 = 3,
  trapnell2014=3,
  zhengmix2016 =4,
  koh2016 = 9,
  simDataKumar=4,
  simDataKumar2=4
  
)
pct_dropout_max2<- list(
  kumar2015 = 90,
  trapnell2014=90,
  zhengmix2016 =99,
  koh2016 = 99,
  simDataKumar=90,
  simDataKumar2=90
  
)
# optimalk
par.k3 <-  list(
  kumar2015 = 3,
  trapnell2014 = 2,
  zhengmix2016=4,
  koh2016 = 11,
  simDataKumar=4,
  simDataKumar2=4
  
)

# which parameter set
if ( (datatype=="filtered")| (datatype=="smooth")) { par.ks  <- par.ks2  ;par.k_estimator <- par.k_estimator2; par.k <- par.k3; pct_dropout_max <- pct_dropout_max2}
else if ((datatype == "unfiltered" )) { par.ks  <- par.ks2  ;par.k_estimator <- par.k_estimator2; par.k <- par.k2; pct_dropout_max <- pct_dropout_max2}
else if ((datatype == "default")) { par.ks  <- par.ks1 ;par.k_estimator <- par.k_estimator1 ; par.k <- par.k1 ; pct_dropout_max <- pct_dropout_max1 }
else if ((datatype == "optimalk")) { par.ks  <- par.ks2  ;par.k_estimator <- par.k_estimator2; par.k <- par.k3; pct_dropout_max <- pct_dropout_max2 }
else {print("datatype not defined") }

print(par.ks )
print(par.k_estimator)
print(par.k)

# check data files and parameters
stopifnot( names(files) == names(data) )
stopifnot( names(par.ks) == names(data) )
# RUN SC3
run_function_sc3( data[paste0(dataset)], labels[paste0(dataset)], par.ks[paste0(dataset)], par.k_estimator[paste0(dataset)] ,par.k[paste0(dataset)], pct_dropout_max[paste0(dataset)] ,datatype )
}

### Appendix


