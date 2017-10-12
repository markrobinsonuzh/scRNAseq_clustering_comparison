#####################
# Seurat
#####################

source( "skript/helper_files/Helper_functions.R" )


#load libraries

library(Seurat)

# file paths

DATA_DIR <- "data"
files <- list(
  
  kumar2015 = file.path(DATA_DIR, "sceset_red_GSE60749-GPL13112.rda" ),
  trapnell2014 = file.path(DATA_DIR, "sceset_red_GSE52529-GPL16791.rda" ),
  xue2013 = file.path(DATA_DIR, "sceset_red_GSE44183-GPL11154.rda" ),
  koh2016 = file.path(DATA_DIR,"sceset_red_SRP073808.rda" )
  
)

# load data sets

data <- labels <-res.cluster <-  sys.time <-  vector("list", length(files))

names(data) <- names(res.cluster) <- names(labels) <- names(sys.time) <- names(files)


for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}

# Seurat
# Resolution parameter, higher number gives more cluster, lower less cluster.
par.resolution <- list(
    kumar2015 = round(seq(0.1, 2,length.out = 10),1),
    trapnell2014 = round(seq(0.1, 2,length.out = 10),1),
    xue2013 = rep(0.6,10),
    koh2016 = round(seq(0.1, 2,length.out = 10),1)
  )
par.dims.use <-  list(
  kumar2015 = 1:10,
  trapnell2014 = 1:10,
  xue2013 = 1:10,
  koh2016 = 1:10
)

#Seurat function

run_seurat <- function( data, par.resolution, par.dims.use ) {
  require(Seurat)
  res.cluster <-  sys.time <-  vector("list", length(files))
  names(res.cluster) <- names(sys.time) <- names(files)
### Run SEurat
for (i in names(data)) {
  # create Seurat object
  data[[i]] <- CreateSeuratObject(raw.data = counts(data[[i]]), min.cells = -Inf, min.genes = -Inf, project = "test") # use raw count_lstpm
  ## Normalizing the data. After removing unwanted cells from the dataset, 
  # the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" 
  #that normalizes the gene expression measurements for each cell by the total expression, 
  #multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
  data[[i]] <- NormalizeData(object = (data[[i]]), normalization.method = "LogNormalize", scale.factor = 1e4)
  # Detection of variable genes across the single cells
  data[[i]] <- FindVariableGenes(object =(data[[i]]), mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
  ### Scaling the data and removing unwanted sources of variation
  data[[i]] <- ScaleData(object = data[[i]])

### Perform lin dimension reduction and cluster the cells
  ### Perform linear dimensional reduction
  sys.time[[i]] <- system.time({
    
    data[[i]] <- RunPCA(object = data[[i]], pc.genes = data[[i]]@var.genes, do.print = FALSE)
    #data[[i]] <- JackStraw(object = data[[i]], num.replicate = 30, do.print = FALSE)
    ### Cluster the cells
    # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
    # but with a different resolution value (see docs for full details)
    df.clus <- matrix( nrow = ncol(data[[i]]@raw.data), ncol = length( par.resolution[[i]]  ))
    for ( j in seq_len( length(par.resolution) ) ){

    df.clus[,j] <- FindClusters(object = data[[i]], reduction.type = "pca", 
                            dims.use = par.dims.use[[i]] , resolution = par.resolution[[i]][j] , 
                            print.output = 0, save.SNN = TRUE  )@ident
                           
    } 
  })
  colnames(df.clus) <- c( paste0("par.res", par.resolution[[i]]) )
  res.cluster[[i]] <-  df.clus
}
  return( res.cluster )
}
# Run seurat

res.cluster <-  run_seurat( data, par.resolution, par.dims.use )

# save clusters

dir_cluster <- paste0("results/Seurat/seurat_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/Seurat/seurat_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("results/Seurat/seurat_krange_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
  
}


###### Save Session Info
sink(file = "results/Seurat/session_info_Seurat_krange.txt")
sessionInfo()
sink()

##### Apendix




