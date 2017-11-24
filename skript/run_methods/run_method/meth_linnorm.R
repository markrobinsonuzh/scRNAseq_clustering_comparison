
#####################
# Linnorm
#####################
# to do : include spikeinns


source("skript/helper_files/Helper_functions.R")


# source file paths: fileterd , raw etc.

source("FILES.R")
# source method linnomr
source("skript/run_methods/run_functions/run_function_linnorm.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# num_center is number of k in kmeans, default is k = c(1:20).
# Set a lower cutoff for the zhengmix data
par.minNonZeroPortion <- list(
  kumar2015 = 0.75,
  trapnell2014 =  0.75,
  zhengmix2016 =  0.25,
  koh2016 =  0.75,
  simDataKumar=0.75
)
par.num_center <- list(
  kumar2015 = 1:20,
  trapnell2014 =  1:20,
  zhengmix2016 = 1:20,
  koh2016 = 1:20,
  simDataKumar=1:20
)
# define datatype
datatype <- "default"
#check if files, parameters and data are the same:
names(files)==names(data) 

# RUN cidr
run_function_linnorm(data, labels, datatype , par.minNonZeroPortion, par.num_center)
