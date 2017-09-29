#####################
# tSNE + kmeans
#####################

# The following method uses the Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding for dimensionality reduction and kmeans foe clustering.
# Set Random seed for reproducible results.  Rtsne uses by  PCA for dimensionality reduction, by default 50 dimension are retained.
# In Rtsne the perplexity parameter has to be defined, typical values are between 5 and 50. Perplexity parameter is loosly speaking a guess about the number of close neighbors each point has. 
# Perplexity has to at least to be smaller than the number of points. 
# Other hyperparameters as the learning rate epsilon and the number of iteration which can give different results for different value ranges. 

#load libraries
source("skript/helper_files/Helper_functions.R")
library(Rtsne)
library(scater)
library(dplyr)

# file paths

DATA_DIR <- "data"
files <- list(
  
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_SRP073808.rda")
  
)

# load data sets

list<- vector("list", length(files))
names(list) <- names(files)

list->data->labels->tinput_matrix->sys.time->res.rtsne->res.cluster 

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}
# extract transposed expression data

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}

# RUN tSNE and kmeans
# set random seed
rand.seed <- 1234

# define the perplexity parameter for tSNE
par.perp <- list(
  kumar2015 = 20,
  trapnell2014 = 20,
  xue2013 = 5,
  koh2016 = 20
)
# define the number of cluster for kmeans clustering 
par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  xue2013 = 8,
  koh2016= 10
)


# Run tSNE and kmeans
# run tSNE
res.rtsne[[1]] <- Rtsne(X = tinput_matrix[[1]] , perplexity = par.perp[[1]] , pca = TRUE)
#      

res.cluster[[i]] <- as.character(kmeans(res.rtsne[[i]]$Y, centers = par.k[[i]])$cluster)






# save clusters

dir_cluster <- paste0("results/RtSNEkmeans/RtSNEkmeans_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/RtSNEkmeans/RtSNEkmeans_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/RtSNEkmeans/RtSNEkmeans_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/RtSNEkmeans/session_info_kmeans.txt")
sessionInfo()
sink()

### Appendix

