#####################
# tSNE + kmeans
#####################

# The following method uses the Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding for dimensionality reduction and kmeans foe clustering.
# Set Random seed for reproducible results.  Rtsne uses by  PCA for dimensionality reduction, by default 50 dimension are retained.
# In Rtsne the perplexity parameter has to be defined, typical values are between 5 and 50. Perplexity parameter is loosly speaking a guess about the number of close neighbors each point has. 
# Perplexity has to at least to be smaller than the number of points. 
# Other hyperparameters as the learning rate epsilon and the number of iteration which can give different results for different value ranges. 
# Variables:  datatype: which analysis; default , filtered, unfiltered


analyze_rtsnekmeans <- function( datatype, dataype){ 
  # set random seed
  rand.seed = 1234
  
  # source method rtsnekmeans and helpers
  source("skript/run_methods/run_functions/run_function_rtsnekmeans.R") 
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
  # define the number of cluster for kmeans clustering 
  par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016= 9,
  simDataKumar=4
  )
  # define the perplexity parameter  and the number of dimension from PCA for tSNE, 
  # default
  par.perp.1 <- list(
    kumar2015 = 30,
    trapnell2014 = 30,
    zhengmix2016 = 30,
    koh2016 = 30,
    simDataKumar=30
  )
  par.initial_dims.1 <- list(
    kumar2015 = 30,
    trapnell2014 = 30,
    zhengmix2016 = 30,
    koh2016= 30,
    simDataKumar=30
    
  )
  # filtered , unfiltered
  par.perp.2<- list(
    kumar2015 = 30,
    trapnell2014 = 30,
    zhengmix2016 = 30,
    koh2016 = 30,
    simDataKumar=30
  )
  par.initial_dims.2 <- list(
    kumar2015 = 20,
    trapnell2014 =20,
    zhengmix2016 = 20,
    koh2016= 20,
    simDataKumar=20
    
  )
  # which parameter set
if ((datatype == "unfiltered" ) | (datatype=="filtered")) { par.perp <- par.perp.2 ; par.initial_dims <-  par.initial_dims.2
} else {
  if ((datatype == "default")) { par.perp <- par.perp.1 ; par.initial_dims <-  par.initial_dims.1 }
  else {print("datatype not defined") }
}
 print(par.perp)
 print(par.initial_dims)
 

  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.perp) == names(data) )

  
  # Run tSNE and kmeans
  run_function_rtsnekmeans(  data[paste0(dataset)], labels[paste0(dataset)], par.k[paste0(dataset)], par.perp[paste0(dataset)], par.initial_dims[paste0(dataset)],datatype )
  
}