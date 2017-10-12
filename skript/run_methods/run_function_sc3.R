####################
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


DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_red_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_red_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_red_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_red_SRP073808.rda")
)

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
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}


# RUN SC3


# QC data
# done in scater

# run the analysis
par.ks <- list(
  kumar2015=c(2:10),
  trapnell2014=c(2:10),
  xue2013=c(2:10),
  koh2016 =c(2:10)
  
)

run_sc3 <- function( data, par.ks ){
  # list to store results
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <- sys.time<- list
for (i in names(data)){
 
    data[[i]]<- sc3_prepare(data[[i]], ks = par.ks[i])  # uses the exprs slot of SCEset ; log2transformed, normalized data, filter data 
    #data[[i]]<- sc3_estimate_k(data[[i]]) # optional estimate the number of clusters
    data[[i]]<- sc3(data[[i]], ks = par.ks[[i]]) # perform sc3 clustering

  # store clusters
  p_data <- pData(data[[i]])
  res.cluster[[i]] <- p_data[ , grep("sc3_", colnames(p_data))]
  
}

return(res.cluster)
}
# mach dataframe
# run the function

res.cluster <-  run_sc3( data, par.ks )

# save clusters

dir_cluster <- paste0("results/SC3/sc3_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SC3/sc3_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/SC3/sc3_krange_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}
###### Save Session Info
sink(file = "results/SC3/session_info_sc3_krange.txt")
sessionInfo()
sink()

### Appendix
#

