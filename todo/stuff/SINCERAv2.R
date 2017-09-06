#####################
# SINCERA
#####################

library(SINCERA)
library(SC3)
library(pheatmap)
library(scRNA.seq.funcs)

# import data as sceset
# file paths
DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds")
)

# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# create storage files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- input_matrix<- pca.red <- ztrans.matrix <- dist.matrix <-hc <-    list



# extract expression data
for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]])
}

##### SINCERA analysis


# transform data
# perform gene-by-gene per-sample z-score transformation
for(i in 1:length(input_matrix)){
  ztrans.matrix[[i]] <- apply(input_matrix[[i]], 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
}

# hierarchical clustering
# create dissimilarity matrix
for(i in 1:length(input_matrix)){
dist.matrix[[i]] <- as.dist((1 - cor(t(ztrans.matrix[[i]]), method = "pearson"))/2)
}
# clustering

for(i in 1:length(dist.matrix)){
hc[[i]] <- hclust(dist.matrix[[i]], method = "average")
}

# stop
# !whats happenng here; slice clusters so that no singleton exists
num.singleton <- 0
kk <- 1
for (i in 2:dim(ztrans.matrix[[1]])[2]) {
  clusters <- cutree(hc[[1]], k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)

# store variable
pData(data[[1]])$SINCERA <- as.character(cutree(hc[[1]], k = i))
plotPCA(data[[1]], colour_by = "SINCERA")
#
# Let's now visualize the SINCERA results as a heatmap:
pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = 14,
  kmeans_k = 100,
  show_rownames = FALSE
)


