#________________________________
#Compute silhouette widths for datasets
#_________________________________
library(cluster)
library(scater)

# Overview Table for datasets

datasets <- c("Kumar", "Trapnell" ,"Koh",
              "SimKumar4easy", "SimKumar4hard", "SimKumar8hard", 
              "KohTCC", "TrapnellTCC", "KumarTCC", 
              "Zhengmix4eq", "Zhengmix4uneq" ,"Zhengmix8eq")
filterings <-  c("filteredExpr10","filteredHVG10", "filteredM3Drop10")

dir <-"/Volumes/Shared/data/seq/scRNAseq_clustering_comparison/data/"

# load full data
list.full <-as.list( paste0(dir,"sce_full/","sce_full_", datasets, ".rds") )
names(list.full) <- datasets
full_data <- lapply(list.full, function(x) readRDS(x)  )
# compute silhouette widths

# Eucl. distances from transposed count matrix
s <- lapply(full_data,function(x){
  d <- dist( t(scater::exprs(x))) 
  s <- cluster::silhouette(
    as.integer(
      as.factor(
        colData(x)$phenoid)),
    d)
  ssum <- summary(s)
  return(ssum)
})

#save somewhere
saveRDS(s, file = "ouput/res_silhouette.rds")
