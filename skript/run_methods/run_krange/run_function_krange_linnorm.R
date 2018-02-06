
#####################
# Linnorm
#####################
# to do : include spikeinns
## try http:// if https:// URLs are not supported
library("Linnorm")
library(scater)
require(dplyr)
set.seed(1234)
source("skript/helper_files/Helper_functions.R")


# file paths

source("FILES.R")
# load data sets

data <- res.cluster <- vector("list", length(files))

names(data) <- names(res.cluster) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# load cell labels
labels <- load_labels(data) 


# RUN linnorm

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# Set a lower cutoff for the zhengmix data
par.minNonZeroPortion <- list(
  kumar2015 = 0.75,
  trapnell2014 =  0.75,
  zhengmix2016 =  0.1,
  koh2016 =  0.75,
  simDataKumar=0.75,
  simDataKumar2=0.75
  

)
par.num_center <- list(
  kumar2015 = c(2:10),
  trapnell2014 = c(2:10),
  koh2016 = c(2:15),
  zhengmix2016=c(2:10),
  simDataKumar=c(2:10),
  simDataKumar2=c(2:10)
  

)

# function 
run_lindworm <- function(data, par.num_center, par.minNonZeroPortion){
  require("Linnorm")
  # create store vectors
  list<- vector("list", length(files))
  names(list) <- names(files)
  list->sys.time->transformedExp->res.cluster 
  # run function
  for (i in names(data)){
    print(i)
    transformedExp[[i]] <- Linnorm(counts(data[[i]]), spikein=NULL, minNonZeroPortion = par.minNonZeroPortion[[i]])
    df.clus <- matrix( nrow = ncol( transformedExp[[i]] ), ncol=length(par.num_center[[i]]) ) 
    for ( j in seq_len(length(par.num_center[[i]])) ){
      print(paste0("k=",j+1))
      df.clus[,j]<- Linnorm.tSNE(transformedExp[[i]], input = "Linnorm",num_center=par.num_center[[i]][j] )$k_means$cluster
    }
    colnames(df.clus) <-  c( paste0(par.num_center[[i]]) )
    res.cluster[[i]] <-  df.clus
  }
  return(res.cluster)
}

# run function

res.cluster <- run_lindworm( data, par.num_center, par.minNonZeroPortion )


# save clusters
dir_cluster <- paste0("results/filtered/linnorm/linnorm_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save experiment labels
labels <- labels[1:5]
file_names <-  paste0("results/filtered/linnorm/linnorm_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/filtered/linnorm/session_info_krange_linnorm.txt")
sessionInfo()
sink()

### Appendix

