######################################
# Interface functions for bootclust  #
######################################
#-----------------------------------------------------
# RtsneKmeans
rtsnekmeansCBI <- function(  data, k, perplexity ) {
  require(Rtsne)
  data <- t(data)
  res.rtsne <- Rtsne(X= data , perplexity=perplexity , pca = TRUE, check_duplicates = FALSE)
  res.cluster <- as.integer( kmeans(res.rtsne$Y, centers = k )$cluster )
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "rtsnekmeansCBI",
    cluster=cbind(cluster, id=rownames(data) )
  )
  return(list)

}

#-----------------------------------------------------
# CIDR
cidrCBI <- function ( data, par.k,par.nPC) {
  require(cidr)
  # RUN cidr
  sData <- vector("list", length(data))
  names(sData) <-  names(data)
  # define number of clusters.
      sData <- scDataConstructor(data) # creates a scData object with slots for the expression table, lib size, dropout candidates etc...
      sData <- determineDropoutCandidates(sData) # determines the dropout candidates
      sData <- wThreshold(sData, plotTornado = FALSE) # sets the imputation weighting threshold
      sData <- scDissim(sData) # computes the dissimilarity matrix for the slot dissim
      sData <- scPCA(sData) # performs PCA on the dissimilarity matrix
      sData <- nPC(sData) # deterimines the optimal number of PC to be used in clustering, populates nPC
      # nCluster(sData) # different methods todefine the number of clusters, optional
      sData <- scCluster(object = sData, nCluster = par.k, cMethod = "ward.D2", nPC = par.nPC)# hierarchical clustering on PC , nCluster if user defines cluster number, nPC is the number of PC used for clustering (default is 4), cMethod is hierarchical clustering method default is ward.D2
      res.cluster <- sData@clusters
      # out
      clusterlist <- vector("list", length( unique(res.cluster) ))
      clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
      cluster <- as.data.frame(res.cluster)
      list <- list(
        result=res.cluster,
        nc= length( unique(res.cluster) ),
        clusterlist=clusterlist,
        partition=res.cluster,
        clustermethod= "cidrCBI",
        cluster=cbind(cluster, id=colnames(data) )
      )
        return(list)
        
    
  }
#-----------------------------------------------------
# Linnorm
linnormCBI <- function ( data, par.minNonZeroPortion, par.num_center ,par.BE_strength) {
  require(Linnorm)
  transformedExp <- Linnorm(data, spikein=NULL, minNonZeroPortion = par.minNonZeroPortion, BE_strength=par.BE_strength)
  res.cluster <- Linnorm.tSNE(transformedExp, input = "Linnorm",num_center=par.num_center )$k_means$cluster
  # out
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "linnormCBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
  
}
#-----------------------------------------------------
# SIMLR

simlrCBI <- function (data,  par.c, normalize) {
  require(SIMLR)
  require(igraph)
  require(scater)
  
  # RUN SIMLR
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.SIMLR <- res.cluster <-  list
  # run SIMLR
  res.SIMLR = SIMLR(X =data, c = par.c, no.dim = NA,k=10, if.impute=FALSE, normalize = normalize ) # use exprs slot of SCeset; log2, normalized count_lstpm
  res.cluster <- res.SIMLR$y$cluster  
  # out
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "simlrCBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
}

#-----------------------------------------------------
# SC3

sc3CBI <- function ( data, par.ks, par.k_estimator ,par.k, pct_dropout_max, rand_seed ) {
  require("DESeq2")
  require("SC3")
  sceset <- SingleCellExperiment( assays = list(counts = data))
  assay(sceset, "exprs") <- data # tpm is populated with the gene level tpms
  counts(sceset)
  
  names(assays(sceset) ) 
  
  # list to store results
  list <- vector("list", length(data))
  names(list) <- names(data)
  res.cluster <- list

  data <- SingleCellExperiment(data)
  data<- sc3_prepare(sceset, ks=par.ks)# uses the exprs slot of SCEset ; log2transformed, non normalized data, filter data 
  
  data<- sc3_estimate_k(sceset) # optional estimate the number of clusters
  # estimated k (default ) or user supplied k (filtered, unfiltered)
  ifelse(par.k_estimator==TRUE, par.k <- metadata(data)$sc3$k_estimation, par.k <- par.k)
  # run sc3
  data<- sc3(data, ks = par.k , pct_dropout_max=pct_dropout_max,
             gene_filter = TRUE ,rand_seed=NULL, n_cores=2, rand_seed=rand_seed) # perform sc3 clustering
  
  # store clusters
  p_data <- colData(data)
  res.cluster <- p_data[ , grep("sc3_", colnames(p_data))]
  # out
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "sc3CBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
  
}
#-----------------------------------------------------
# pcaReduce
pcareduceCBI <- function ( data, par.nbt, par.q ,n.cluster) {
  require(pcaReduce)
  # extract k dimension 
  
  par.k <- function(i){
    (par.q+2)-(n.cluster)
  }
  
  # run pce Reduce, vary q
      pca.red <- PCAreduce(t(data), nbt = par.nbt, q = par.q, method = 'S')[[1]]
      res.cluster  <- as.character(pca.red[ ,par.k(i)])
      # out
      clusterlist <- vector("list", length( unique(res.cluster) ))
      clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
      cluster <- as.data.frame(res.cluster)
      list <- list(
        result=res.cluster,
        nc= length( unique(res.cluster) ),
        clusterlist=clusterlist,
        partition=res.cluster,
        clustermethod= "pcareduceCBI",
        cluster=cbind(cluster, id=colnames(data) )
      )
      return(list)
  
}

# pcaReduce with consensus
pcareduceCBI2<- function ( data, par.nbt, par.q ,n.cluster) {
  require(pcaReduce)
  # extract k dimension 
  
  par.k <- function(i){
    (par.q+2)-(n.cluster)
  }
  
  # run pce Reduce, vary q
  pca.red <- PCAreduce(t(data), nbt = par.nbt, q = par.q, method = 'S')
  part <- lapply(pca.red, function(x)x[,par.k(i)])
  part<- lapply(part, function(x) as.cl_partition(x))
  part <- as.cl_ensemble(part)
  cl_cons <- cl_consensus(  part , method = "SE",  control = list(nruns=50))
  res.cluster<- as.character( cl_class_ids( cl_cons ) )
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "pcareduceCBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
  
}
#-----------------------------------------------------
# Seurat
seuratCBI <- function( data, par.resolution, k.param , par.dims.use , r.seed) {
require(Seurat)
res.cluster <-  vector("list", length(data))
names(res.cluster) <- names(data)

### Run SEurat
# normalize data
  # create Seurat object
  data <- CreateSeuratObject(raw.data = data, min.cells = 0, min.genes = 0, project = "whatever") # use raw count_lstpm
  ## Normalizing the data. After removing unwanted cells from the dataset, 
  # the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" 
  #that normalizes the gene expression measurements for each cell by the total expression, 
  #multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
  data <- NormalizeData(object = (data), normalization.method = "LogNormalize", scale.factor = 1e4)
  # Detection of variable genes across the single cells
  data <- FindVariableGenes(object =(data), mean.function = ExpMean, dispersion.function = LogVMR)
  ### Scaling the data and removing unwanted sources of variation
  data <- ScaleData(object = data)
  


### Perform lin dimension reduction and cluster the cells

  ### Perform linear dimensional reduction
  
    data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = FALSE)
    #data <- JackStraw(object = data, num.replicate = 100, do.print = FALSE)
    ### Cluster the cells
    # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
    # but with a different resolution value (see docs for full details)
    data <- FindClusters(object = data, reduction.type = "pca", dims.use = par.dims.use, k.param = k.param ,
                              resolution = par.resolution, print.output = 0, save.SNN = FALSE, random.seed = r.seed)
 
   res.cluster <-  data@ident
   # out
   clusterlist <- vector("list", length( unique(res.cluster) ))
   clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x-1 } )# here minus one to corecct for zero in cluster
   cluster <- as.data.frame(res.cluster)
   list <- list(
     result=res.cluster,
     nc= length( unique(res.cluster) ),
     clusterlist=clusterlist,
     partition=res.cluster,
     clustermethod= "seuratCBI",
     cluster=cbind("res.cluster"=data@ident,"id"=names(data@ident) )# use colnames here, as it is alread named
     
   )
   return(list)
}


#-----------------------------------------------------
# TSCAN
tbscanCBI <- function ( data, par.minexpr_percent  ,par.clusternum) {
  
  require(TSCAN)

  # create store vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  list->res.tscan->res.cluster 
  res.tscan <- preprocess(data, minexpr_percent = par.minexpr_percent, logbase = 2) # preprocessing
  res.cluster <- exprmclust(res.tscan, clusternum = 
                                   par.clusternum )$clusterid # clustering
  # out
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "tscanCBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
}
#-----------------------------------------------------
# raceid


raceidCBI <- function ( data, par.mintotal, par.maxexpr,do.gap,cln, r.seed) {
  
source("method_resources/RaceID/RaceID_class.R")
source("skript/helper_files/Helper_functions.R")
require(tsne)
require(pheatmap)
require(MASS)
require(cluster)
require(mclust)
require(flexmix)
require(lattice)
require(fpc)
require(amap)
require(RColorBrewer) 
require(locfit)
  

list<- vector("list", length(data))
names(list) <- names(data)
list->sc->res.cluster 
  
    sc <- SCseq(as.data.frame(data)) # extract the expression data) # extract the expression data
    sc <- filterdata(sc, mintotal = par.mintotal, minexpr = 5, 
                          minnumber = 1, maxexpr = par.maxexpr , downsample = FALSE, dsn = 1, rseed = r.seed)

    res.cluster<- clustexp(sc, metric="pearson", cln=cln, do.gap=do.gap, clustnr=20, B.gap=50, SE.method="Tibs2001SEmax", 
                                SE.factor=.25, bootnr=50, rseed=r.seed)@kmeans$kpart
    # out
    clusterlist <- vector("list", length( unique(res.cluster) ))
    clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
    cluster <- as.data.frame(res.cluster)
    list <- list(
      result=res.cluster,
      nc= length( unique(res.cluster) ),
      clusterlist=clusterlist,
      partition=res.cluster,
      clustermethod= "raceidCBI",
      cluster=cbind(cluster, id=colnames(data) )
    )
    return(list)
}
#-----------------------------------------------------
# zibwave

zinbwaveCBI <- function( data, par.k  , n.genes) {
  
  require(zinbwave)
  require(Rtsne)
  filter_hvg <- function(data, n.genes){
    filter <- rowSums( ( data  )>5)>5 
    data <- data[filter,]# filter out genes with at least five counts
    data  %>% log1p %>% rowVars -> vars
    names(vars) <- rownames( data )
    vars <- sort(vars, decreasing = TRUE)
    data <- data[names(vars)[1:n.genes],]
    return(data)
  }

  # create vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  list->transformedExp-> res.zinb -> d -> tsne_data -> res.cluster
 
    # filter for selction of hvg genes
    data <- filter_hvg( data, n.genes )
    # run ZINBWaVE

      res.zinb <- zinbFit( round(data,0) ,
                                K=2, epsilon=n.genes, verbose=TRUE,
                                nb.repeat.initialize = 2, 
                                maxiter.optimize = 25,
                                stop.epsilon.optimize = 1e-04)  # round data as it assumes whole counts
      
      d<- dist(getW( res.zinb ))
      #tsne_data <- Rtsne(d, is_distance = TRUE, pca = FALSE, perplexity=10, max_iter=5000)
      res.cluster <- kmeans(d, centers=par.k )$cluster
  
  # out
  clusterlist <- vector("list", length( unique(res.cluster) ))
  clusterlist <- lapply(seq_len(length( unique(res.cluster))), function(x){ res.cluster == x } )
  cluster <- as.data.frame(res.cluster)
  list <- list(
    result=res.cluster,
    nc= length( unique(res.cluster) ),
    clusterlist=clusterlist,
    partition=res.cluster,
    clustermethod= "zinbwaveCBI",
    cluster=cbind(cluster, id=colnames(data) )
  )
  return(list)
}
####

