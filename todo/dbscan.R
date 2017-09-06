#####################
# DBscan
#####################
source("~/Desktop/masterthesis/skript/Helper_functions.R")


#load libraries

library(pcaReduce)

# import data as sceset
# file paths

DATA_DIR <- "~/Desktop/masterthesis/data"
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
res.cluster <- sys.time<- input_matrix<- pca.red <- list


# extract expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}
dim(input_matrix[[1]])
# RUN dbscan
kNNdistplot(input_matrix[[1]], k = 97)
dbscan(input_matrix[[1]], eps = 260 ,minPts = 5)









# number of clusters in kmeans
rand.seed <- 1234
par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  xue2013 = 8
)

par.perp <- list(
  kumar2015 = 20,
  trapnell2014 = 20,
  xue2013 = 5
)
# Run tSNE and kmeans
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- list


for (i in names(data)){
  
  sys.time [[i]] <- system.time({
    data[[i]] <- plotTSNE(data[[i]], exprs_values= "counts",rand_seed = rand.seed, perplexity= par.perp[[i]],return_SCESet = TRUE, draw_plot= FALSE) # use Rtsne? function plotTSNE uses Rtsne anyway
    pData(data[[i]])$tSNE_kmeans <- as.character(kmeans(data[[i]]@reducedDimension, centers = par.k[[i]])$clust)
  })
  res.cluster[[i]] <- pData(data[[i]])$tSNE_kmean
  
}


# save clusters

dir_cluster <- paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "~/Desktop/masterthesis/results/tSNEkmeans/session_info_kmeans.txt")
sessionInfo()
sink()

### Appendix
# plot Data

# with Rtsne
#library(Rtsne)
#rtsne <- Rtsne(dd, perplexity =5)
#plot(rtsne$Y)
