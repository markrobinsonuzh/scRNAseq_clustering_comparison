
#####################
# TSCAN
#####################
# TSCAN first preprocess the data by filter out rarley expressed cells and non DE features.
# The TSCAN::preprocess function additionally takes the log plus a pseudocount. function will rule out genes which
# have expression values of less than 1 in at least half of all cells. Genes with covariance smaller than one for the expression values are as well filterd out.
# Next pseudotime analysis is done, using dim reduction with PCA and model based clustering. A  pseudo temporal ordering score (POS) and travelling sales- man problem algorithm (TSP) algorithm combined with 
# Due to included preprocessing we use raw counts. Here the default transformation with log base2 , a pseudocount of one and cutoff for min expr of 2 
# is chosen. The method has an addiotionally criteria for the covariances for features, default is 1. We keep only high expressed genes, here we keep at least 5 pwercent of all genes for highly expressed cells.



analyze_tscan <- function( datatype, dataset){ 
  
source("skript/helper_files/Helper_functions.R")

  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")
  } else {
    if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
    if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
    
    else {print("datatype not defined") }
  }
  
# source method tscan
source("skript/run_methods/run_functions/run_function_tscan.R")

# load data sets
data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# Set a lower cutoff for the zhengmix data
# default
par.minexpr_percent1 <- list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  zhengmix2016 = 0.1,
  koh2016 = 0.5,
  simDataKumar = 0.5
)

par.clusternum1 <- list(
  kumar2015 = 2:10,
  trapnell2014 = 2:10,
  zhengmix2016 =2:10,
  koh2016 =2:10,
  simDataKumar = 2:10
)
# filtered , unfiltered
par.minexpr_percent2 <- list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  zhengmix2016 = 0.1,
  koh2016 = 0.5,
  simDataKumar = 0.5
)

par.clusternum2 <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016= 10,
  simDataKumar=4
)
# optimalk
par.minexpr_percent3 <- list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  zhengmix2016 = 0.1,
  koh2016 = 0.5,
  simDataKumar = 0.1
)

par.clusternum3 <- list(
  kumar2015 = 3,
  trapnell2014 = 2,
  zhengmix2016 = 3,
  koh2016= 10,
  simDataKumar=4
)

# which parameter set
if ((datatype == "unfiltered" ) | (datatype=="filtered")| (datatype=="smooth")) { par.minexpr_percent  <- par.minexpr_percent2  ; par.clusternum <-  par.clusternum2
} else {
  if ((datatype == "default")) { par.minexpr_percent  <- par.minexpr_percent1  ; par.clusternum <-  par.clusternum1 }
  if ((datatype == "optimalk")) { par.minexpr_percent  <- par.minexpr_percent3  ; par.clusternum <-  par.clusternum3 }
  else {print("datatype not defined") }
}
print(par.minexpr_percent )
print(par.clusternum)


# check data files and parameters
stopifnot( names(files) == names(data) )
stopifnot( names(par.clusternum) == names(data) )

# RUN TSCAN
run_function_tscan( data[paste0(dataset)], labels[paste0(dataset)],  par.minexpr_percent[paste0(dataset)]  , par.clusternum = par.clusternum[paste0(dataset)]  , datatype )
}
