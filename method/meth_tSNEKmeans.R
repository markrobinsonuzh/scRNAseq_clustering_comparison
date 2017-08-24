#####################
# tSNE + kMeans
#####################

#load libraries

library(scater)
library(Rtsne)
set.seed(1234567)

# file paths

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rds")
  
)

# load data sets

data <- labels<- vector("list", length(files))

names(data) <-names(labels) <-  names(files)

for (i in 1:length(data)){
    f <- files[[i]]
    load(f)
    data[[i]] <- res
    
}

# load cell labels
for(i in 1:length(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$source_name_ch1)
}


# RUN tSNE and kmeans

# number of clusters in kmeans
rand.seed <- 1234
par.k <- list(
  kumar2015 <- 3,
  trapnell2014 <- 12,
  xue2013 <- 8
)

par.perp <- list(
  kumar2015 <- 20,
  trapnell2014 <- 20,
  xue <- 5
)
# Run tSNE and kmeans
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- list



for (i in 1:length(data)){
  
  sys.time [[i]] <- system.time({
  data[[i]] <- plotTSNE(data[[i]], exprs_values= "counts",rand_seed = rand.seed, perplexity= par.perp[[i]],return_SCESet = TRUE, draw_plot= FALSE) # use Rtsne? function plotTSNE uses Rtsne anyway
  pData(data[[i]])$tSNE_kmeans <- as.character(kmeans(data[[i]]@reducedDimension, centers = par.k[[i]])$clust)
  })
  res.cluster[[i]] <- pData(data[[i]])$tSNE_kmean
  
}


# save clusters

file_names <- paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_clus_", 
                       names(res.cluster), ".txt")

for (i in 1:length(file_names)) {
  res_i <- res.cluster[[i]]
  write.table(res_i, file = file_names[i], row.names = FALSE, quote = FALSE, sep = "\t")
}

# save systemtime

file_names <-  paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_systime_",names(sys.time))
for (i in 1:length(sys.time)){
  sys_i <- sys.time[[i]]["elapsed"]
  write.table(sys_i, file=file_names[i], sep="\t")
  
}


# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/tSNEkmeans/tSNEkmeans_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
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


