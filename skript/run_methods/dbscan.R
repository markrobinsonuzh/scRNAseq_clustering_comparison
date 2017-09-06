#####################
# DBscan
#####################

source("~/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/WORKIN_DIR.R")
source(paste0(WORKIN_DIR,"skript/helper_files/Helper_functions.R"))



#load libraries

library(dbscan)

# import data as sceset
# file paths

DATA_DIR <- paste0(WORKIN_DIR,"data")
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
)

# load data sets

data <- labels<- vector("list", length(files))

names(data) <-names(labels) <-  names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}


# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}

# create store files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- res.dbscan <- sys.time<- input_matrix<- pca.red <- list


# extract transposed expression data

for (i in names(input_matrix)){
  input_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}
# RUN dbscan

par.k <- list(
  kumar2015 = 97,
  trapnell2014 = 97,
  xue2013 = 15
)


# run k neirest neighbour distance plot to find epsilon
par(mfrow=c(2,2))
for (i in names(input_matrix)){

  kNNdistplot(input_matrix[[i]], k = par.k[[i]]) 
  
}

# run dbscan
par.eps <- list(
  kumar2015 = 280,
  trapnell2014 = 400,
  xue2013 = 600
)
par.Pts <- list(
  kumar2015 = 5,
  trapnell2014 = 5,
  xue2013 = 5
)

for(i in names(input_matrix)){
res.cluster[[i]] <- dbscan(input_matrix[[i]], eps = par.eps[[i]] ,minPts = par.Pts[[i]])$cluster

}



# save clusters

dir_cluster <- paste0(WORKIN_DIR,"results/dbscan/dbscan_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0(WORKIN_DIR,"results/dbscan/dbscan_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0(WORKIN_DIR,"results/dbscan/dbscan_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = paste0(WORKIN_DIR,"results/dbscan/session_info_dbscan.txt"))
sessionInfo()
sink()

### Appendix
# plot Data

# with Rtsne
#library(Rtsne)
#rtsne <- Rtsne(dd, perplexity =5)
#plot(rtsne$Y)
