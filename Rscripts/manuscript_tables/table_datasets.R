# Overview Table for datasets
library(cluster)
library(scater)

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

# number of cells, number of features per dataset, median counts per cell, median gnees per cell, number of suppopulations
ncells <- sapply(full_data, function(x) ncol(x)) # ncell
ngenes <- sapply(full_data, function(x) nrow(x)) # nfeatures
mcounts <-sapply(full_data, function(x) median(x$total_counts, na.rm = TRUE)/1e6) # median counts per cell (Mio)
mgenes <- sapply(full_data, function(x) median(x$total_features, na.rm = TRUE))# median features per cell, genes with zero expression excluded
npop <- sapply(full_data, function(x) length( unique(SummarizedExperiment::colData(x)$phenoid) ) )# n populations per dataset

# avg.silhoutte width
readRDS(file = "res_silhouette.rds")
avg.s <- sapply(s,function(x)x$avg.width)

# table
tbl <- cbind(ncells, mcounts, mgenes,ngenes, npop, avg.s )

print( xtable::xtable(tbl) )
write.csv(tbl, file="tbldata.csv")

