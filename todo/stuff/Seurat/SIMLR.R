# SIMLR
library("SIMLR")
#browseVignettes("SIMLR")
library(igraph)
library(scater)

load("/Users/angeloduo/Desktop/masterthesis/data/sceset_GSE60749_GPL13112.rds")
load("/Users/angeloduo/Desktop/masterthesis/data/sceset_GSE52529-GPL16791.rds")

set.seed(11111)
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

### with SCEset
output = SIMLR(X = res, c = 12, cores.ratio = 0)
plot(output$ydata,
     xlab = "SIMLR component 1",
     ylab = "SIMLR component 2",
     pch = 20,
     col= output$y$cluster,
     main="SIMILR 2D visualization for GSE60749-GPL13112")


