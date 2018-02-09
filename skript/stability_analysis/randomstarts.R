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

# copy sample 30 times

sample_rep <- function(data) {
  data[]
}

#bootstrap
sample.norm <- NULL
for ( i in 1:30) {
  sample.norm[[i]] <-  data.normcts  
}
sample.raw <- NULL
for ( i in 1:30) {
  sample.raw[[i]] <-  data.raw  
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

# adaption of function for methods with internal random seed parameters
r.stability.analysis <- function(fun.method,sample.full,subsample, ... ) {
  # merge function
  merge.tbl <- function(res.sub,res.org){
    merge(res.sub,res.org, by="id" )
  }
  #vector with randomstarts
  r.seed <- sample(1:100,30,replace=F)
  # clustering
  res.org <- fun.method(sample.full, r.seed=r.seed[1], ... )$cluster%>%as.data.frame
  res.sub <- mapply( fun.method, subsample, MoreArgs =list(r.seed= r.seed), ..., SIMPLIFY = FALSE)%>%lapply("[[", 6)
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
res.raceid <-  r.stability.analysis(fun.method=raceidCBI,sample.full=data.normcts, subsample=sample.norm, par.mintotal=3000, par.maxexpr=Inf ,do.gap=FALSE,cln=3  )
res.rtsnekmeans <-  stability.analysis( fun.method=rtsnekmeansCBI,sample.full=data.normcts, subsample=sample.norm ,k=3, perplexity=30)
res.seurat <-  r.stability.analysis( fun.method=seuratCBI, sample.full=data.normcts, subsample=sample.norm ,par.resolution=0.6, k.param = 25 , par.dims.use=1:9)

#res.sc3 <-  stability.analysis(fun.method=sc3CBI,sample.full=data.normcts, subsample=sample.norm, par.ks = NULL, par.k_estimator=FALSE , par.k =3, pct_dropout_max =90, rand_seed=NULL )
res.zinbwave <-  stability.analysis(fun.method=zinbwaveCBI,sample.full=data.normcts, subsample=sample.norm, par.k=3, n.genes=1000 )
res.pcareduce <-  stability.analysis( fun.method=pcareduceCBI, sample.full=data.normcts, subsample=sample.norm,par.nbt=100, par.q=30 ,n.cluster=3)


#-------------------------------------------------------------
# save ari results

write.table(res.cidr , file = "results/rstart_analysis/res.cidr.txt", sep="\t")
write.table(res.linnorm , file = "results/rstart_analysis/res.linnorm.txt", sep="\t")
write.table(res.tscan  , file = "results/rstart_analysis/res.tscan.txt", sep="\t")
write.table(res.simlr , file = "results/rstart_analysis/res.simlr.txt", sep="\t")
write.table(res.raceid , file = "results/rstart_analysis/res.raceid.txt", sep="\t")
write.table(res.rtsnekmeans , file = "results/rstart_analysis/res.rtsnekmeans.txt", sep="\t")
write.table(res.seurat  , file = "results/rstart_analysis/res.seurat.txt", sep="\t")
write.table(res.zinbwave  , file = "results/rstart_analysis/res.zinbwave.txt", sep="\t")
write.table(res.pcareduce  , file = "results/rstart_analysis/res.pcareduce.txt", sep="\t")
#-------------------------------------------------------------

### Appendix

# seurat
stability.analysis.seurat <- function(fun.method,sample.full,subsample, ... ) {
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
  res <- lapply(x, function(x) adjustedRandIndex(x$res.cluster.x ,x$res.cluster.y ) )
  return(unlist(res))
  
}
