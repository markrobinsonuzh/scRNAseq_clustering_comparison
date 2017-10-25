#####################
# tSNE + kmeans optimal number of cluster
#####################

# The following method uses the Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding for dimensionality reduction and kmeans foe clustering.
# Set Random seed for reproducible results.  Rtsne uses by  PCA for dimensionality reduction, by default 50 dimension are retained.
# In Rtsne the perplexity parameter has to be defined, typical values are between 5 and 50. Perplexity parameter is loosly speaking a guess about the number of close neighbors each point has. 
# Perplexity has to at least to be smaller than the number of points. 
# Other hyperparameters as the learning rate epsilon and the number of iteration which can give different results for different value ranges. 
pdf("results/plots/optimalk_wss_tsnekmeans.pdf")
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
  xue2013 = 5,
  koh2016 = 20
)
# define the number of cluster for kmeans clustering 
## maximal number of Clusters
par.nk <- list(
  kumar2015=10,
  trapnell2014=15,
  xue2013=15,
  koh2016=15
)


# Run tSNE and kmeans
# run tSNE
for (i in names(tinput_matrix)){
res.rtsne[[i]] <- Rtsne(X = tinput_matrix[[i]] , perplexity = par.perp[[i]] , pca = TRUE, dim = 10)
}

# compute the wss
wss <- vector("list", length(files))
names(wss) <- names(files)

for (i in names(wss)){
  wss[[i]] <- (nrow(res.rtsne[[i]]$Y )-1)*sum(apply(res.rtsne[[i]]$Y,2,var))
}

# compute wss for different k
for (i in names(wss)) {{
  for (j in 2:par.nk[[i]] ) wss[[i]][j] <- sum(kmeans(res.rtsne[[i]]$Y,
                                                      centers=j)$withinss)
  
}}
# save the data
save(wss,file="results/number_k/rescluswsstsne.rda")

 
# plot
par(mfrow=c(2,2))
for ( i in 1:length(files)) {
  plot(1:15, wss[[i]][1:15], type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares", main=paste0(names(files)[i]))
}

dev.off()