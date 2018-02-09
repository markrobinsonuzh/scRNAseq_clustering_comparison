########################################
# SC3
#########################################




source("skript/helper_files/Helper_functions.R")

# source file paths: fileterd , raw etc.

source("FILES.R")

# source method SC3
source("skript/run_methods/run_functions/run_function_sc3.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 
# run the analysis
par.ks <- list(
  kumar2015 = 2:6,
  trapnell2014=2:6,
  zhengmix2016 =2:6,
  koh2016 = 2:12,
  simDataKumar=2:6
)

par.k_estimator <- list(
  kumar2015 = TRUE,
  trapnell2014=TRUE,
  zhengmix2016 =TRUE,
  koh2016 =TRUE,
  simDataKumar=TRUE
)

par.k <- list(
  kumar2015 = 3,
  trapnell2014=3,
  zhengmix2016 =4,
  koh2016 = 10,
  simDataKumar=4
)
# define datatype
datatype <- "default"
# check if files, parameters and data are the same:
names(files)==names(data) 

# RUN SC3
run_function_sc3( data, labels, par.ks, par.k_estimator ,par.k,datatype )

### Appendix


