#####################################
## Gap statistic for datasets      ##  
#####################################

### load libraries and helper files
source("skript/helper_files/Helper_functions.R")

library(cluster)
library(dplyr)
library(fpc)
library(scater)

# define plot directories
pdf("results/plots/plot_optimalk_gap.pdf")

### set seed

set.seed(1234)

# load data
source("FILES.R")


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
  kumar2015=100,
  trapnell2014=100,
  zhengmix2016=100,
  koh2016=100
  
)

par.K.max <- list(
  kumar2015=12,
  trapnell2014=12,
  zhengmix2016=12,
  koh2016=12
)


for (i in names(input_matrix))
res.clusgap[[i]] <- clusGap(input_matrix[[i]], kmeans, K.max=par.K.max[[i]], B = par.B[[i]] , verbose = interactive(), spaceH0 = "scaledPCA")

save(res.clusgap,file="results/number_k/resclusgap.rda")
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