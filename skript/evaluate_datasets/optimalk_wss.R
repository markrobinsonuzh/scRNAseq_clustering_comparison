#####################################
## Within sum of squares for clusters      
#####################################
pdf("~/Desktop/masterthesis/results/plots/optimalk_wss.pdf")
### load libraries
source("~/Desktop/masterthesis/skript/helper_functions/Helper_functions.R")
library(cluster)
library(dplyr)
library(fpc)
library(scater)

### set seed

set.seed(1234)

# Directories

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
  
)


# load data sets

data <- vector("list", length(files))
tinput_matrix<- data
names(data) <- names(tinput_matrix) <-  names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# extract transposed expression data

for (i in 1:(length(tinput_matrix))){
  tinput_matrix[[i]] <- t(exprs(data[[i]])) # use count scaled length scaled tpms, normalized and log2 transformed
}



#####################################
## Within sum of squares for clusters     
#####################################

## maximal number of Clusters
par.nk <- list(
  kumar2015=15,
  trapnell2014=15,
  xue2013=15
  
)
# 
wss <- vector("list", length(files))
names(wss) <- names(files)
for (i in names(wss)){
wss[[i]] <- (nrow(tinput_matrix[[i]] )-1)*sum(apply(tinput_matrix[[i]],2,var))
}




# compute wss for different k
for (i in 1:3) {{
  for (j in 2:par.nk[[i]] ) wss[[i]][j] <- sum(kmeans(tinput_matrix[[i]],
                                       centers=j)$withinss)

}}

# plot
par(mfrow=c(2,2))
for ( i in 1:3) {
plot(1:15, wss[[i]][1:15], type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main=paste0(names(files)[i]))
}

dev.off()
