#####################
# Seurat
#####################

analyze_seurat <- function(dataype, dataset){
  
source("skript/helper_files/Helper_functions.R")
# source method CIDR
source("skript/run_methods/run_functions/run_function_seurat.R")
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
# Resolution parameter resolution, higher number gives more cluster, lower less cluster. we define as the standart by 0.8
# k.param is the number of neirest neighbors
# the number of PC dim to use was determined by an elbow plot and by the jackstraw function


# fixed parameters for all run modes
k.param<- list(
  kumar2015 =ncol(data[["kumar2015"]])* c(0.1),
  trapnell2014 = ncol(data[["trapnell2014"]])*c(0.1),
  zhengmix2016 = ncol(data[["zhengmix2016"]])*c(0.1),
  koh2016 = ncol(data[["koh2016" ]])*c(0.1),
  simDataKumar=ncol(data[["simDataKumar"]])*c(0.1),
  simDataKumar2=ncol(data[["simDataKumar2"]])*c(0.1)
  
)
k.param<- lapply(k.param,round,0)
# resolution parameter
# optimalresolution based on range from 0.6 to 1.2
par.resolution <-  list(
  kumar2015 = 0.6,
  trapnell2014 = 0.6,
  zhengmix2016=0.6,
  koh2016 = 0.7,
  simDataKumar=0.6,
  simDataKumar2=0.6
  
)
# parrameters for no cell filtering
par.mincells1 <- list(
  kumar2015 =0,
  trapnell2014 = 0,
  zhengmix2016 = 0,
  koh2016 = 0,
  simDataKumar=0,
  simDataKumar2=0
  
)
par.mingenes1 <- list(
  kumar2015 =0,
  trapnell2014 = 0,
  zhengmix2016 = 0,
  koh2016 = 0,
  simDataKumar=0,
  simDataKumar2=0
  
)
# default run mode:
k.param1 <- list(
  kumar2015 =30,
  trapnell2014 = 30,
  zhengmix2016 = 30,
  koh2016 = 30,
  simDataKumar=30,
  simDataKumar2=30
  
)
par.dims.use1 <-  list(
  kumar2015 = NULL,
  trapnell2014 = NULL,
  zhengmix2016 = NULL,
  koh2016 = NULL,
  simDataKumar=NULL,
  simDataKumar2=NULL
  
)
#unfiltered dataset:
par.mincells2 <- list(
  kumar2015 =2,
  trapnell2014 = 2,
  zhengmix2016 = 2,
  koh2016 = 2,
  simDataKumar=2,
  simDataKumar2=2
  
)
par.mingenes2 <- list(
  kumar2015 =0,
  trapnell2014 = 0,
  zhengmix2016 = 0,
  koh2016 = 0,
  simDataKumar=0,
  simDataKumar2=0
  
)

#filtered mode:

par.dims.use <-  list(
  kumar2015 = 1:9,
  trapnell2014 = 1:12,
  zhengmix2016 = 1:10,
  koh2016 = 1:15,
  simDataKumar=1:10,
  simDataKumar2=1:10
  
)
# optimalk:

par.dims.use <-  list(
  kumar2015 = 1:9,
  trapnell2014 = 1:12,
  zhengmix2016 = 1:10,
  koh2016 = 1:15,
  simDataKumar=1:10,
  simDataKumar2=1:10
  
)


# which parameter set:
if ( (datatype=="filtered") | (datatype=="smooth")) {par.resolution <- par.resolution; k.param <- k.param;  par.dims.use <- par.dims.use; par.mingenes <-par.mingenes1 ; par.mincells <- par.mincells1 }
else if ((datatype == "unfiltered" ) ) {par.resolution <- par.resolution; k.param <- k.param;  par.dims.use <- par.dims.use; par.mingenes <-par.mingenes2 ; par.mincells <- par.mincells2}
else if ((datatype == "default")) { par.resolution <- par.resolution; k.param <- k.param1;  par.dims.use <- par.dims.use; par.mingenes <-par.mingenes1 ; par.mincells <- par.mincells1 }
else if ((datatype == "optimalk")) { par.resolution <- par.resolution ; k.param <- k.param;  par.dims.use <- par.dims.use; par.mingenes <-par.mingenes1 ; par.mincells <- par.mincells1 }
else {print("datatype not defined") }

# run Seurat
run_function_seurat(  data[paste0(dataset)], labels[paste0(dataset)],par.mincells[paste0(dataset)] ,par.mingenes[paste0(dataset)],par.resolution[paste0(dataset)], k.param[paste0(dataset)] , par.dims.use[paste0(dataset)],  datatype )

}


