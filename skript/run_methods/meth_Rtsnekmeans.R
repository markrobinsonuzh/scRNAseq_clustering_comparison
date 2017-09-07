#####################
# tSNE + kmeans
#####################

#load libraries
source("skript/helper_files/Helper_functions.R")

library(Rtsne)
library(scater)
# file paths

DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
  
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

# define perplexitz parameter for tsne an number of clusters in kmeans
rand.seed <- 1234
par.perp <- list(
  kumar2015 = 20,
  trapnell2014 = 20,
  xue2013 = 5
)

par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 3,
  xue2013 = 8
)


# Run tSNE and kmeans


for (i in names(data)){
  
  sys.time [[i]] <- system.time({
    res.rtsne[[i]] <- Rtsne(X= tinput_matrix[[i]] ,perplexity=par.perp[[i]] , pca = TRUE, verbose = TRUE)
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
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/RtSNEkmeans/session_info_kmeans.txt")
sessionInfo()
sink()

### Appendix

