#####################
# Seurat
#####################

analyze_seurat <- function(dataype, dataset){
  
source("skript/helper_files/Helper_functions.R")
# source method CIDR
source("skript/run_methods/run_functions/run_function_seurat.R")
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

# parameters
# Resolution parameter resolution, higher number gives more cluster, lower less cluster. we define as the standart by 0.8
# k.param is the number of neirest neighbors
# the number of PC dim to use was determined by an elbow plot and by the jackstraw function
# default:
k.param1 <- list(
  kumar2015 =30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016 = 30,
  simDataKumar=30
)


par.dims.use1 <-  list(
  kumar2015 = NULL,
  trapnell2014 = NULL,
  zhengmix2016 = NULL,
  koh2016 = NULL,
  simDataKumar=NULL
)

#filtered, unfilterd;
k.param2 <- list(
  kumar2015 =ncol(data[["kumar2015"]])* c(0.1),
  trapnell2014 = ncol(data[["trapnell2014"]])*c(0.1),
  zhengmix2016 = ncol(data[["zhengmix2016"]])*c(0.1),
  koh2016 = ncol(data[["koh2016" ]])*c(0.1),
  simDataKumar=ncol(data[["simDataKumar"]])*c(0.1)
)

k.param2<- lapply(k.param2,round,0)

par.dims.use2 <-  list(
  kumar2015 = 1:9,
  trapnell2014 = 1:12,
  zhengmix2016 = 1:10,
  koh2016 = 1:15,
  simDataKumar=1:10
)
# which parameter set
if ((datatype == "unfiltered" ) | (datatype=="filtered")) {k.param <- k.param2;  par.dims.use <- par.dims.use2
} else {
  if ((datatype == "default")) { k.param <- k.param1;  par.dims.use <- par.dims.use1 }
  else {print("datatype not defined") }
}
print(k.param)
# run Seurat
run_function_seurat(  data[paste0(dataset)], labels[paste0(dataset)], k.param[paste0(dataset)] , par.dims.use[paste0(dataset)],  datatype )

}
