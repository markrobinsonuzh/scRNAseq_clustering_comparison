##########################################
# Seurat range of resolution parameter
##########################################
set.seed(1234)
source( "skript/helper_files/Helper_functions.R" )


#load libraries

library(Seurat)

# file paths
source("FILES.R")


# load data sets

data <- load_data(files, DATA_DIR )


# load cell labels
labels <- load_labels(data)

# Seurat
# Resolution parameter resolution, higher number gives more cluster, lower less cluster. we define as the standart by 0.8
# k.param is the number of neirest neighbors
# the number of PC dim to use was determined by an elbow plot and by the jackstraw function

k.param <- list(
  kumar2015 =ncol(data[["kumar2015"]])* c(0.1),
  trapnell2014 = ncol(data[["trapnell2014"]])*c(0.1),
  zhengmix2016 = ncol(data[["zhengmix2016"]])*c(0.1),
  koh2016 = ncol(data[["koh2016"]])*c(0.1),
  simDataKumar =ncol(data[["simDataKumar"]])*c(0.1),
  simDataKumar2 =ncol(data[["simDataKumar2"]])*c(0.1)
  
)
par.dims.use <-  list(
  kumar2015 = 1:9,
  trapnell2014 = 1:12,
  zhengmix2016 = 1:10,
  koh2016 = 1:15,
  simDataKumar=1:10,
  simDataKumar2=1:10
)
par.resolution <- list(
  kumar2015 = c(0.6,0.7,0.8,0.9,1.0,1.1,1.2),
  trapnell2014 =c(0.6,0.7,0.8,0.9,1.0,1.1,1.2),
  zhengmix2016 = c(0.6,0.7,0.8,0.9,1.0,1.1,1.2),
  koh2016 = c(0.6,0.7,0.8,0.9,1.0,1.1,1.2),
  simDataKumar =c(0.6,0.7,0.8,0.9,1.0,1.1,1.2),
  simDataKumar2 =c(0.6,0.7,0.8,0.9,1.0,1.1,1.2)
  
  
)
#Seurat function
run_seurat <- function( data,  k.param, par.resolution, par.dims.use ) {
  
  require(Seurat)
  res.cluster <-  sys.time <-  vector("list", length(files))
  names(res.cluster) <- names(sys.time) <- names(files)
  ### Run SEurat
  for (i in names(data)) {
    # create Seurat object
    data[[i]] <- CreateSeuratObject(raw.data = counts( data[[i]] ), min.cells = 0, min.genes = 0, project = "krange") # use raw count_lstpm
    ## Normalizing the data. After removing unwanted cells from the dataset, 
    # the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" 
    #that normalizes the gene expression measurements for each cell by the total expression, 
    #multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
    data[[i]] <- NormalizeData(object = (data[[i]]), normalization.method = "LogNormalize", scale.factor = 1e4)
    # Detection of variable genes across the single cells
    data[[i]] <- FindVariableGenes(object =(data[[i]]), mean.function = ExpMean, dispersion.function = LogVMR, do.plot =FALSE)
    ### Scaling the data and removing unwanted sources of variation
    data[[i]] <- ScaleData(object = data[[i]])
    
    ### Perform lin dimension reduction and cluster the cells
    ### Perform linear dimensional reduction
    
    
    data[[i]] <- RunPCA(object = data[[i]], pc.genes = data[[i]]@var.genes, do.print = FALSE)
    data[[i]] <- JackStraw(object = data[[i]], num.replicate = 100, do.print = FALSE)
    JackStrawPlot(data[[i]] , PCs=1:15)
    PCElbowPlot(object =   data[[i]])
    
    ### Cluster the cells
    # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
    # but with a different resolution value (see docs for full details)
    df.clus <- matrix( nrow = ncol(data[[i]]@raw.data), ncol = length( par.resolution[[i]]  ))
    for ( j in seq_len( length(par.resolution) ) ){
      df.clus[,j] <- FindClusters(object = data[[i]], reduction.type = "pca", k.param = k.param[[i]],
                                  dims.use = par.dims.use[[i]] , resolution = par.resolution[[i]][j] , 
                                  print.output = 0, save.SNN = TRUE  )@ident
      
    } 
    
    colnames(df.clus) <- c( paste0(  par.resolution[[i]]) )
    res.cluster[[i]] <-  df.clus
  }
  return( res.cluster )
}
# Run seurat

res.cluster <-  run_seurat( data, k.param, par.resolution, par.dims.use )

# save clusters

dir_cluster <- paste0("results/filtered/Seurat/seurat_krange_clus_", names(res.cluster), ".txt")
save_clusters(res.cluster,dir_cluster)



# save experiment labels

file_names <-  paste0("results/filtered/Seurat/seurat_krange_labels_",names(res.cluster), ".txt")
for (i in 1:length(res.cluster)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
  
}


###### Save Session Info
sink(file = "results/filtered/Seurat/session_info_Seurat_krange.txt")
sessionInfo()
sink()

##### Apendix


