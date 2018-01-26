#####################
# Seurat
#####################


run_function_seurat <-  function(  data, labels, par.mincells ,par.mingenes ,par.resolution , k.param , par.dims.use,  datatype ){
  
  require(Seurat)
  
  res.cluster <-  sys.time <-  vector("list", length(data))
  names(res.cluster)  <- names(sys.time) <- names(data)
  
### Run SEurat
# normalize data
for (i in names(data)) {
  print(i)
  # create Seurat object
  data[[i]] <- CreateSeuratObject(raw.data = counts( data[[i]] ), min.cells = par.mincells[[i]], min.genes = par.mingenes[[i]], project = "scRNAseq") # use raw count_lstpm
  ## Normalizing the data. After removing unwanted cells from the dataset, 
  # the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" 
  #that normalizes the gene expression measurements for each cell by the total expression, 
  #multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
  data[[i]] <- NormalizeData(object = (data[[i]]), normalization.method = "LogNormalize", scale.factor = 1e4)
  # Detection of variable genes across the single cells
  data[[i]] <- FindVariableGenes(object =(data[[i]]), mean.function = ExpMean, dispersion.function = LogVMR)
  ### Scaling the data and removing unwanted sources of variation
  data[[i]] <- ScaleData(object = data[[i]])
  
}

### Perform lin dimension reduction and cluster the cells
for (i in names(data)){
  ### Perform linear dimensional reduction
sys.time[[i]] <- system.time({
data[[i]] <- RunPCA(object = data[[i]], pc.genes = data[[i]]@var.genes, do.print = FALSE)
#data[[i]] <- JackStraw(object = data[[i]], num.replicate = 100, do.print = FALSE)
PCElbowPlot(object =   data[[i]])
### Cluster the cells
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
# but with a different resolution value (see docs for full details)
data[[i]] <- FindClusters(object = data[[i]], reduction.type = "pca", dims.use = par.dims.use[[i]], k.param = k.param[[i]] ,
                          resolution = par.resolution[i], print.output = 0, save.SNN = TRUE)
})
res.cluster[[i]] <-  data[[i]]@ident
}

# save clusters

dir_cluster <- paste0("results/",datatype,"/Seurat/seurat_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/Seurat/seurat_systime_",names(res.cluster),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

dir_labels <-  paste0("results/",datatype,"/Seurat/seurat_labels_", names(res.cluster), ".txt")
save_labels(labels, dir_labels )


###### Save Session Info
sink(file = paste0("results/",datatype,"/Seurat/session_info_Seurat.txt") )
sessionInfo()
sink()
}


