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

source("FILES.R")

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
  zhengmix2016 = 20,
  koh2016 = 20
)
# define the number of cluster for kmeans clustering 
par.k <- list(
  kumar2015 = c(2:10),
  trapnell2014 = c(2:10),
  zhengmix2016 = c(2:10),
  koh2016= c(2:10)
)


#tSNE and kmeans function

run_RtSNEkmeans <- function(tinput_matrix,par.k, par.perp) {
for (i in names(data)){
  
 
    res.rtsne[[i]] <- Rtsne(X= tinput_matrix[[i]] ,perplexity=par.perp[[i]] , pca = TRUE)
    
    df.clus <- matrix( nrow = nrow(tinput_matrix[[i]]), ncol = length(par.k[[i]]) ) 
    
    for ( j in seq_len(length(par.k[[i]])) ) {
      df.clus[,j]<- as.integer( kmeans(res.rtsne[[i]]$Y, centers = par.k[[i]][j])$cluster )
    }
    colnames(df.clus) <-  c( paste0(par.k[[i]]) )
    res.cluster[[i]] <- df.clus
}
  return(res.cluster)
}

### run function
res.cluster <- run_RtSNEkmeans(tinput_matrix,par.k, par.perp)


# save clusters

dir_cluster <- paste0("results/RtSNEkmeans/RtSNEkmeans_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/RtSNEkmeans/RtSNEkmeans_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/RtSNEkmeans/RtSNEkmeans_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/RtSNEkmeans/session_info_kmeans_krange.txt")
sessionInfo()
sink()

### Appendix


