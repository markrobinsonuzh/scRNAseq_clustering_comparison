#####################
# Seurat
#####################


source("skript/helper_files/Helper_functions.R")


# source file paths: fileterd , raw etc.

source("FILES.R")
# source method CIDR
source("skript/run_methods/run_functions/run_function_seurat.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# parameters
# Resolution parameter resolution, higher number gives more cluster, lower less cluster. we define as the standart by 0.8
# k.param is the number of neirest neighbors
# the number of PC dim to use was determined by an elbow plot and by the jackstraw function
k.param <- list(
  kumar2015 =30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016 = 30,
  simDataKumar=30
)

k.param <- lapply(k.param,round,0)

par.dims.use <-  list(
  kumar2015 = NULL,
  trapnell2014 = NULL,
  zhengmix2016 = NULL,
  koh2016 = NULL,
  simDataKumar=NULL
)
# define datatype
datatype <- "default"
# check if files, parameters and data are the same:
names(files)==names(data) 

# run Seurat
run_function_seurat(  data, labels, k.param , par.dims.use,  datatype )

