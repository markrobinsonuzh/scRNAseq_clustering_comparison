
#####################
# TSCAN
#####################
# TSCAN first preprocess the data by filter out rarley expressed cells and non DE features.
# The TSCAN::preprocess function additionally takes the log plus a pseudocount. function will rule out genes which
# have expression values of less than 1 in at least half of all cells. Genes with covariance smaller than one for the expression values are as well filterd out.
# Next pseudotime analysis is done, using dim reduction with PCA and model based clustering. A  pseudo temporal ordering score (POS) and travelling sales- man problem algorithm (TSP) algorithm combined with 
# Due to included preprocessing we use raw counts. Here the default transformation with log base2 , a pseudocount of one and cutoff for min expr of 2 
# is chosen. The method has an addiotionally criteria for the covariances for features, default is 1. We keep only high expressed genes, here we keep at least 5 pwercent of all genes for highly expressed cells.

require(TSCAN)
require(dplyr)

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


# function for  TSCAN

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# Set a lower cutoff for the zhengmix data
par.minexpr_percent <- list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  zhengmix2016 = 0.1,
  koh2016 = 0.5
)
run_tscan <-  function(data,par.minexpr_percent){
# create store vectors
list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->res.tscan->res.cluster 

names(data)
for (i in names(data)){
  
    res.tscan[[i]] <- preprocess(counts(data[[i]]), minexpr_percent = par.minexpr_percent[[i]], logbase = 2) # preprocessing
    res.cluster[[i]] <- exprmclust(res.tscan[[i]])$clusterid # clustering

  }
}
# run tscan
res.cluster <-  run_tscan( data, par.minexpr_percent  ) 


# save clusters

dir_cluster <- paste0("results/tscan/tscan_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/tscan/tscan_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/tscan/tscan_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/tscan/session_info_tscan_krange.txt")
sessionInfo()
sink()

### Appendix
