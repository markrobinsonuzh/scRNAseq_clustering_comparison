#####################################
## Gap statistic for datasets      ##  
#####################################
# Workin directory
source("~/Desktop/masterthesis/scRNAseq_clustering_comparison/skript/helper_files/WORKIN_DIR.R")

### load libraries
source(paste0(WORKIN_DIR,"skript/helper_files/Helper_functions.R"))

library(cluster)
library(dplyr)
library(fpc)
library(scater)

pdf(paste0(WORKIN_DIR,"results/plots/plot_optimalk_gap.pdf"))

### set seed

set.seed(1234)

# load data

DATA_DIR <- paste0(WORKIN_DIR,"data")
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
  
)


# load data sets

data <- vector("list", length(files))
input_matrix<- data
names(data) <- names(input_matrix) <-  names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# extract transposed expression data

for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <- exprs(data[[i]]) # use count scaled length scaled tpms, normalized and log2 transformed
}
input_matrix[[i]] <- t(input_matrix[[i]])
################################################
##### Compute gap statistic, number of Monte Carlo samples = 100 (increase it to 1000) , maximum of clusters K.max=10
################################################

# define number of Monte Carlo samples par.B and number of clusters par.K.max
res.clusgap <- vector("list", length(files))
names(res.clusgap)  <-  names(input_matrix)

par.B <- list(
  kumar2015=1,
  trapnell2014=1,
  xue2013=1
  
)

par.K.max <- list(
  kumar2015=2,
  trapnell2014=2,
  xue2013=2
)


for (i in names(input_matrix))
res.clusgap[[i]] <- clusGap(input_matrix[[i]], kmeans, K.max=par.K.max[[i]], B = par.B[[i]] , verbose = interactive(), spaceH0 = "scaledPCA")

save(res.clusgap,file=paste0(WORKIN_DIR,"results/number_k/resclusgap.rda"))
#### store results for further use



#### plot the the gap statistic
plots_gap <- names(files) <- vector("list", length(files))
plot_gap_statistic <- function(gaps, stddevs, num_clusters) {
  qplot(num_clusters, gaps, xlab = "number of clusters", ylab = "Gap", geom = "line", main = "Estimating the number of clusters via the gap statistic")+geom_errorbar(aes(num_clusters, ymin = gaps - stddevs, ymax = gaps + stddevs), size = 0.3, width = 0.2, colour = "darkblue")
}

for (i in names(res.clusgap)){
plots_gap[[i]] <- plot_gap_statistic(gap = res.clusgap[[i]]$Tab[,3], stddevs = res.clusgap[[i]]$Tab[,4], num_clusters = seq_len(par.K.max[[i]]) )
print(plots_gap[[i]])
}
dev.off()