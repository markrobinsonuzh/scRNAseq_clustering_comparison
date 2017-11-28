####################a
# SC3
######################
# change sc 3 prepare
# Input in SC3 is by default the exprs slot in a SCEset object. So the data should be normalized and log transformed. AS an option one can also set the argument gene_filter, if a filtering beforhand should br performed.
# SC3 uses distance measures of the filtered and log transformed expression matrix. PCA or Laplacian graphs are used for  dimension reduction. 
# Kmeans clustering is then performed on the d different dimensions. Using the d different clustering results a consensus matrix is computed. On the distances of this consenus matrix a hierarchical
# clustering step is performed. A range of number of clusters k  can be used by used. Its also possible to estimate the number of clusters.
# 
source("skript/helper_files/Helper_functions.R")

# load libraries

library("scater")
library("DESeq2")
library("SC3")

# file paths


source("FILES.R")


# load data sets

data <- vector("list", length(files))
labels <- data
names(data) <- names(labels) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
labels <- load_labels(data) 



# RUN SC3

# number of clusters k
par.k <- list(
  kumar2015 = c(2:10),
  trapnell2014 = c(2:10),
  koh2016 = c(2:15),
  zhengmix2016=c(2:10),
  simDataKumar=c(2:10)
  
)
# run function
run_sc3 <- function( data, par.k ){
  # list to store results
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <- sys.time<- list
for (i in names(data)){
 
    data[[i]]<- sc3_prepare(data[[i]], ks = par.k[i])  # uses the exprs slot of SCEset ; log2transformed, normalized data, filter data 
    #data[[i]]<- sc3_estimate_k(data[[i]]) # optional estimate the number of clusters
    data[[i]]<- sc3(data[[i]], ks = par.k[[i]], gene_filter=FALSE) # perform sc3 clustering

  # store clusters
  p_data <- colData(data[[i]])
  res.cluster[[i]] <- p_data[ , grep("sc3_", colnames(p_data))]
  colnames(res.cluster[[i]]) <-  c( paste0(par.k[[i]]) )
}

return(res.cluster)
}
# mach dataframe
# run the function

res.cluster <-  run_sc3( data, par.k )

# save clusters

dir_cluster <- paste0("results/filtered/SC3/sc3_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/filtered/SC3/sc3_krange_systime_",names(res.cluster),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/filtered/SC3/sc3_krange_labels_",names(res.cluster), ".txt")
for (i in 1:length(res.cluster)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}
###### Save Session Info
sink(file = "results/filtered/SC3/session_info_sc3_krange.txt")
sessionInfo()
sink()

### Appendix
#

