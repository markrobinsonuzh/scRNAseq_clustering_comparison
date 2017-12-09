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

source("FILESraw.R")

# load data sets

data <- load_data(files, DATA_DIR)
# load cell labels
labels <- load_labels(data) 

# extract transposed expression data
list<- vector("list", length(files))
names(list) <- names(files)

list->data->labels->tinput_matrix->sys.time->res.rtsne->res.cluster 

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <-  t(assay(data[[i]], "normcounts"))# use count scaled length scaled tpms, normalized and log2 transformed
}


# RUN tSNE and kmeans
# set random seed
rand.seed <- 1234
names(tinput_matrix)
# define the perplexity parameter for tSNE
par.perp <- list(
  raw.kumar2015 = 20,
  raw.trapnell2014 = 20,
  raw.zhengmix2016 = 20,
  raw.koh2016 = 20,
  raw.simDataKumar=20
)
# define the number of cluster for kmeans clustering 
par.k <- list(
  raw.kumar2015 = 3,
  raw.trapnell2014 = 3,
  raw.zhengmix2016 = 4,
  raw.koh2016= 10,
  raw.simDataKumar=3
  
)
# Run tSNE and kmeans
for (i in names(data)){
  
  sys.time[[i]] <- system.time({
    res.rtsne[[i]] <- Rtsne(X= tinput_matrix[[i]] ,perplexity=par.perp[[i]] , pca = TRUE)
    res.cluster[[i]] <- as.character(kmeans(res.rtsne[[i]]$Y, centers = par.k[[i]])$cluster)
  })

}


# save clusters

dir_cluster <- paste0("results/RtSNEkmeans/RtSNEkmeans_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/RtSNEkmeans/RtSNEkmeans_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/RtSNEkmeans/RtSNEkmeans_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  lab_i <- as.data.frame(labels[[i]])
  write.table(lab_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/RtSNEkmeans/session_info_kmeans.txt")
sessionInfo()
sink()

### Appendix

