#######################################################
# Stability analysis of methods using subsampling
# for one dataset only
#######################################################

# load librareis
library(dplyr)
library(mclust)
library(purrr)
# load helper files
source("skript/helper_files/Helper_functions.R")
#load methods for stability anlysis
source("skript/stability_analysis/interface_clusterboot.R")

# which dataset
dataset <- "kumar2015"

#-------------------------------------------------------------

# load a data set and the labels thereof
source("FILES.R")
data <- load_data(files, DATA_DIR)
labels <- load_labels(data)[[paste0(dataset)]] %>% as.factor %>% as.integer
data.normcts <- exprs(data[[paste0(dataset)]]) %>% rbind(labels=labels)
data.raw <-  assay(data[[paste0(dataset)]], "counts") %>% rbind(labels=labels)

#-------------------------------------------------------------

# subsample the dataset 30 times
# function

sample_rep <- function(data) {
  data[, sample(1:ncol(data),100, replace = FALSE)]
}

#bootstrap
sample.norm <- NULL
for ( i in 1:30) {
  sample.norm[[i]] <-  (sample_rep(data.normcts ) )
}
sample.raw <- NULL
for ( i in 1:30) {
  sample.raw[[i]] <-  (sample_rep(data.raw ) )
}

#-------------------------------------------------------------
# wrapper function for computing ARI from subsamples

# fun.method: clustering method
# sample.full: full dataset, as matrix
# sub.sample: subsample
# ... : parameters for clustering method
stability.analysis <- function(fun.method,sample.full,subsample, ... ) {
  # merge function
  merge.tbl <- function(res.sub,res.org){
    merge(res.sub,res.org, by="id" )
  }
  # clustering
  res.org <- fun.method(sample.full, ... )$cluster%>%as.data.frame
  res.sub <- lapply( subsample, fun.method, ... )%>%lapply("[[", 6)
  #merge list
  x <- lapply(res.sub, merge.tbl,res.org )
  #ari
  res <- lapply(x, function(x)adjustedRandIndex(x$res.cluster.x ,x$res.cluster.y ) )
  return(unlist(res))
 
}
# seurat
seurat.stability.analysis <- function(fun.method,sample.full,subsample, ... ) {
  # merge function
  merge.tbl <- function(res.sub,res.org){
    merge(res.sub,res.org, by="id" )
  }
  # clustering
  res.org <-seuratCBI(data.raw, par.resolution=0.6, k.param = 25 , par.dims.use=1:9 )$cluster%>%as.data.frame
  res.sub <- lapply( sample.raw, seuratCBI, par.resolution=0.6, k.param = 25 , par.dims.use=1:9 )%>%lapply("[[", 6)
  #merge list
  x <- lapply(res.sub, merge.tbl,res.org )
  #ari
  res <- lapply(x, function(x)adjustedRandIndex(x$res.cluster.x ,x$res.cluster.y ) )
  return(unlist(res))
  
}

#-------------------------------------------------------------
# run methods on subsamples, comparing between partitions in subsamples and full dataset

res.cidr <- stability.analysis(fun.method=cidrCBI, sample.full=data.normcts, subsample=sample.norm,par.k=3, par.nPC=4  )
res.tscan <-  stability.analysis(fun.method=tbscanCBI,sample.full=data.normcts, subsample=sample.norm, par.minexpr_percent=0.5  ,par.clusternum=3)
res.linnorm <- stability.analysis( fun.method=linnormCBI, sample.full=data.raw, subsample=sample.raw,par.minNonZeroPortion=0.75, par.num_center=3 ,par.BE_strength=0.5 )

res.simlr <- stability.analysis(fun.method=simlrCBI, sample.full=data.normcts, subsample=sample.norm,par.c=3, normalize=TRUE )
res.raceid <-  stability.analysis(fun.method=raceidCBI,sample.full=data.normcts, subsample=sample.norm, par.mintotal=3000, par.maxexpr=Inf ,do.gap=FALSE,cln=3)
res.rtsnekmeans <-  stability.analysis( fun.method=rtsnekmeansCBI,sample.full=data.normcts, subsample=sample.norm ,k=3, perplexity=30)
res.seurat <-  seurat.stability.analysis( fun.method=seuratCBI, sample.full=data.normcts, subsample=sample.norm ,par.resolution=0.6, k.param = 25 , par.dims.use=1:9)

res.zinbwave <-  stability.analysis(fun.method=zinbwaveCBI,sample.full=data.raw, subsample=sample.raw, par.k=3, n.genes=1000 )
res.pcareduce <-  stability.analysis( fun.method=pcareduceCBI2, sample.full=data.normcts, subsample=sample.norm,par.nbt=100, par.q=30 ,n.cluster=3)


#-------------------------------------------------------------
# save ari results

write.table(res.cidr , file = "results/stability_analysis/res.cidr.txt", sep="\t")
write.table(res.linnorm , file = "results/stability_analysis/res.linnorm.txt", sep="\t")
write.table(res.tscan  , file = "results/stability_analysis/res.tscan.txt", sep="\t")
write.table(res.simlr , file = "results/stability_analysis/res.simlr.txt", sep="\t")
write.table(res.raceid , file = "results/stability_analysis/res.raceid.txt", sep="\t")
write.table(res.rtsnekmeans , file = "results/stability_analysis/res.rtsnekmeans.txt", sep="\t")
write.table(res.seurat  , file = "results/stability_analysis/res.seurat.txt", sep="\t")
write.table(res.zinbwave  , file = "results/stability_analysis/res.zinbwave.txt", sep="\t")
write.table(res.pcareduce  , file = "results/stability_analysis/res.pcareduce.txt", sep="\t")

### Appendix


