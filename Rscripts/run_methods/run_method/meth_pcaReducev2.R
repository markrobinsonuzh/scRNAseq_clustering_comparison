###################
# pcaReduce #######
###################

# pcaReduce uses a PCA and hierarchical clustering to find the number of clusters in the reduced dimension given by PCA. 
# The method expects that large classes of cells ar contained in low dimension PC representation
# and more refined (subsets) of these cells types are contained in higher dimensional PC representations.
# On the latent space (nxq) a kmeans clustering with q+1 clusters is performed. And the PC with the lowest variance is deleted. this process is iteratively repeated until only one single cluster remains.
# The resulting matrix has dimension nxq with q+1 clusters.
# Parameters to define are the number of times the method should be repeated, as pcaReduce is stochastic. Sampling without replacement; we choose nbt = 100, if the number samples is bigger than 100.
# The number of starting principal components q, we choose 50 as the default.
# And the number clusters n which is given by the "ground truth". the stepwise merging of the clusters can be done using sampling based merging (S) or
# merging based on largest probability (M)




analyze_pcareduce <- function( datatype, dataset){ 
  
  source("skript/helper_files/Helper_functions.R")
  # source method pcaReduce
  source("skript/run_methods/run_functions/run_function_pcareduce.R")
  # source file paths: fileterd , raw etc.
  if ((datatype == "default") | (datatype=="filtered")| (datatype=="optimalk")) { source("FILES.R"); print("filtered files")}
  else if ((datatype == "unfiltered")) { source("FILESraw.R"); print("raw files") }
  else if ((datatype == "smooth")) { source("FILESsmooth.R"); print("smooth files") }
  else {print("datatype not defined") }
  
  # load data sets
  
  data <- load_data(files, DATA_DIR)
  
  # load cell labels
  labels <- load_labels(data) 
  
  # parameters
  # set parameters, nbt is number of times to repeat pcareduce; q is number of starting dimensions, n cluster the number of clusters
  # for all
  par.nbt <- list(
    kumar2015 = 100,
    trapnell2014 = 100,
    zhengmix2016 = 100,
    koh2016 = 100,
    simDataKumar=100,
    simDataKumar2=100
    
  )
  par.q <- list(
    kumar2015 = 30,
    trapnell2014 = 30,
    zhengmix2016 = 30,
    koh2016 = 30,
    simDataKumar=30,
    simDataKumar2=30
    
  )

  # for filtered
  n.cluster1 <- list(
    kumar2015=3,
    trapnell2014=3,
    zhengmix2016=4,
    koh2016 = 9,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  # for unfiltered
  n.cluster2 <- list(
    kumar2015=3,
    trapnell2014=3,
    zhengmix2016=4,
    koh2016 = 10,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  # optimalk
  n.cluster3  <-  list(
    kumar2015 = 5,
    trapnell2014 = 3,
    zhengmix2016=5,
    koh2016 = 11,
    simDataKumar=4,
    simDataKumar2=4
    
  )
  if ( (datatype=="filtered")| (datatype=="smooth")) { par.nbt <- par.nbt; par.q <- par.q; n.cluster <- n.cluster1 } 
  else if ((datatype == "default")) { par.nbt <- par.nbt; par.q <- par.q; n.cluster <- n.cluster1 }
  else if ((datatype == "unfiltered")) { par.nbt <- par.nbt; par.q <- par.q; n.cluster <- n.cluster2 }
  else if ((datatype == "optimalk")) {  par.nbt <- par.nbt; par.q <- par.q; n.cluster <- n.cluster3 }
  else {print("datatype not defined") }
  
  
  print(par.nbt)
  print(par.q)
  print( n.cluster)
  # check data files and parameters
  stopifnot( names(files) == names(data) )
  stopifnot( names(par.nbt) == names(data) )

  # RUN pcaReduce
  run_function_pcareduce( data[paste0(dataset)], labels[paste0(dataset)], par.nbt[paste0(dataset)], par.q[paste0(dataset)] , n.cluster[paste0(dataset)] , datatype )

}

