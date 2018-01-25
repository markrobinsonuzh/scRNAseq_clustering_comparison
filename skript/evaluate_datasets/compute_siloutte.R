########################################### 
## Compute the sihlhouette statistics   ##  
########################################### 
# libs
require(cluster)
require(dplyr)
source("skript/helper_files/Helper_functions.R")


#----------------------------------------------------------------------------------------------
### load the datasets
# source original files
source("FILES.R")
# load data sets
data <- load_data(files, DATA_DIR)
# load cell labels, as integer
labels <- load_labels(data) 
lab <- lapply(labels, function(x) {  as.factor(x)%>%as.integer(x)  })
# extract counts
cts <- lapply(data, counts)
#transpose
t.cts <- lapply(cts,t)
# dictances Euclidesn

d <- lapply(t.cts, dist)
save(d, file = "results/distances/d.RData")
#compute distance matrix using Euclidean distances
s <- res.s <- vector("list", length(d))
names(res.s) <- names(s) <- names(d)
for (i in names(d)){

  print(i)
  s[[i]] <- silhouette(x=lab[[i]],dist=d[[i]])
  res.s[[i]] <- summary(s[[i]])
}
save(res.s, file = "results/distances/res.s.RData")

