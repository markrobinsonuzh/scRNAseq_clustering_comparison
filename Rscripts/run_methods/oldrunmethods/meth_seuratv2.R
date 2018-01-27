#####################
# Seurat
#####################


#load libraries
library(Seurat)
library(scater)
source("skript/helper_files/Helper_functions.R")

# file paths

source("FILES.R")
# load data sets

data <- labels <-res.cluster <-  sys.time <-  vector("list", length(files))

names(data) <- names(res.cluster) <- names(labels) <- names(sys.time) <- names(files)


for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
labels <- load_labels(data) 

# Seurat
# Resolution parameter resolution, higher number gives more cluster, lower less cluster. we define as the standart by 0.8
# k.param is the number of neirest neighbors
# the number of PC dim to use was determined by an elbow plot and by the jackstraw function
k.param <- list(
  kumar2015 =ncol(data[[1]])* c(0.1),
  trapnell2014 = ncol(data[[2]])*c(0.1),
  zhengmix2016 = ncol(data[[3]])*c(0.1),
  koh2016 = ncol(data[[4]])*c(0.1),
  simDataKumar=ncol(data[[5]])*c(0.1)
)

k.param <- lapply(k.param,round,0)
par.dims.use <-  list(
  kumar2015 = 1:9,
  trapnell2014 = 1:12,
  zhengmix2016 = 1:10,
  koh2016 = 1:15,
  simDataKumar=1:10
)
### Run SEurat
# normalize data
for (i in names(data)) {
  # create Seurat object
  data[[i]] <- CreateSeuratObject(raw.data = counts(data[[i]]), min.cells = 0, min.genes = 0, project = "test") # use raw count_lstpm
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
                          resolution = 0.8, print.output = 0, save.SNN = TRUE)
})
res.cluster[[i]] <-  data[[i]]@ident
}

# save clusters

dir_cluster <- paste0("results/Seurat/seurat_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/Seurat/seurat_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("results/Seurat/seurat_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
  
}


###### Save Session Info
sink(file = "results/Seurat/session_info_Seurat.txt")
sessionInfo()
sink()

##### Apendix
data[[1]] <- CreateSeuratObject(raw.data = counts(data[[1]]), min.cells = 0, min.genes = 0, project = "test") # use raw count_lstpm




