#####################
# tSNE + kmeans
#####################

# The following method uses the Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding for dimensionality reduction and kmeans foe clustering.
# Set Random seed for reproducible results.  Rtsne uses by  PCA for dimensionality reduction, by default 50 dimension are retained.
# In Rtsne the perplexity parameter has to be defined, typical values are between 5 and 50. Perplexity parameter is loosly speaking a guess about the number of close neighbors each point has. 
# Perplexity has to at least to be smaller than the number of points. 
# Other hyperparameters as the learning rate epsilon and the number of iteration which can give different results for different value ranges. 



source("skript/helper_files/Helper_functions.R")
# set random seed
rand.seed <- 1234
# source file paths: fileterd , raw etc.
source("FILES.R") 
# source method CIDR
source("skript/run_methods/run_functions/run_function_rtsnekmeans.R")
# load data sets
data <- load_data(files, DATA_DIR)
# load cell labels
labels <- load_labels(data) 
# define the perplexity parameter for tSNE
par.perp <- list(
  kumar2015 = 30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016 = 30,
  simDataKumar=30
)
# define the number of cluster for kmeans clustering 
par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  zhengmix2016 = 4,
  koh2016= 10,
  simDataKumar=3
  
)

par.initial_dims <- list(
  kumar2015 = 30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016= 30,
  simDataKumar=30
  
)
# define datatype
datatype <- "default"
# check if files, parameters and data are the same:
names(files)==names(data) 

# Run tSNE and kmeans
run_function_rtsnekmeans(  data, labels, par.k, par.perp, par.initial_dims,datatype )

