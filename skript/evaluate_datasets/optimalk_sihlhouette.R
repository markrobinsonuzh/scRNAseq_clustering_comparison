#############################################
# sihlhuette plot
############################################

#####################################
## Within sum of squares for clusters      
#####################################

### load libraries
source("skript/helper_files/WORKIN_DIR.R")

source("skript/helper_files/Helper_functions.R")
library(cluster)


pdf("results/plots/optimalk_wss.pdf")
### set seed

set.seed(1234)

# Directories

DATA_DIR <- paste0(WORKIN_DIR,"data")

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

#kmeans
df <- pam(tinput_matrix[[1]], k=3)
si <- silhouette(df)
plot(si)
