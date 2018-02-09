##########################
# FILES original Data directory
#########################

#### This script defines the files and their respective file paths directories which will be ued in the cluster analysis and evaluations
#### source it before running all the anlysis scripts

# file paths

DATA_DIR <- "data"
files <- list(
  
  kumar2015 = file.path(DATA_DIR, "sceset_org_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_org_GSE52529-GPL16791.rda"),
  zhengmix2016 = file.path(DATA_DIR, "sceset_org_zhengmix.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_org_SRP073808.rda"),
  simDataKumar = file.path(DATA_DIR,"sceset_org_simDataKumar.rda"),
  simDataKumar2 = file.path(DATA_DIR,"sceset_org_simDataKumar2.rda")
)
