##########################
# FILES Data directory
#########################

#### This script defines the files and their respective file paths directories which will be ued in the cluster analysis and evaluations
#### source it before running all the anlysis scripts

# file paths

DATA_DIR <- "data"
files <- list(
  
  kumar2015 = file.path(DATA_DIR, "sceset_red_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_red_GSE52529-GPL16791.rda"),
  zhengmix2016 = file.path(DATA_DIR, "sceset_red_zhengmix.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_red_SRP073808.rda"),
  simDataKumar = file.path(DATA_DIR,"sceset_red_simDataKumar.rda")
)
