########################
# CIDR
########################


#load libraries
source("skript/helper_files/Helper_functions.R")

library(Rtsne)
library(scater)
library(dplyr)
library(cidr)

# file paths

source("FILES.R")


# load data sets

list<- vector("list", length(files))
names(list) <- names(files)

list->data->labels->tinput_matrix->sys.time->res.cluster 

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
labels <- load_labels(data) 

# extract transposed expression data

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <- t((assay(data[[i]], "normcounts"))) # use count scaled length scaled tpms, normalized and log2 transformed
}

# RUN cidr
sData <- vector("list", length(files))
names(sData) <-  names(files)
# define number of clusters.
par.k <-  list(
  kumar2015 = 3,
  trapnell2014 = 3,
  koh2016 = 10,
  simDataKumar=3
)

for  (i in names(sData)) {

sys.time[[i]] <- system.time({
sData[[i]] <- scDataConstructor(t(tinput_matrix[[i]])) # creates a scData object with slots for the expression table, lib size, dropout candidates etc...
sData[[i]] <- determineDropoutCandidates(sData[[i]]) # determines the dropout candidates
sData[[i]] <- wThreshold(sData[[i]], plotTornado = FALSE) # sets the imputation weighting threshold
sData[[i]] <- scDissim(sData[[i]]) # computes the dissimilarity matrix for the slot dissim
sData[[i]] <- scPCA(sData[[i]]) # performs PCA on the dissimilarity matrix
sData[[i]] <- nPC(sData[[i]]) # deterimines the optimal number of PC to be used in clustering, populates nPC

# nCluster(sData) # different methods todefine the number of clusters, optional
sData[[i]] <- scCluster(object=sData[[i]], nCluster = par.k[[i]], cMethod = "ward.D2")# hierarchical clustering on PC , nCluster if user defines cluster number, nPC is the number of PC used for clustering (default is 4), cMethod is hierarchical clustering method default is ward.D2

})

res.cluster[[i]] <- sData[[i]]@clusters

}


# save clusters

dir_cluster <- paste0( "results/cidr/cidr_clus_", names(res.cluster), ".txt" )


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0( "results/cidr/cidr_systime_",names(sys.time),".txt" )

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0( "results/cidr/cidr_labels_", names(labels), ".txt" )
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info

sink(file = "results/cidr/session_info_cidr.txt")
sessionInfo()
sink()

### Appendix

par(mfrow=c(2,2))

for (i in names(sData)){

  plot(sData[[i]]@PC[,c(1,2)], 
     pch=sData[[i]]@clusters, main=paste0(names(sData[i])), xlab="PC1", ylab="PC2")
}

