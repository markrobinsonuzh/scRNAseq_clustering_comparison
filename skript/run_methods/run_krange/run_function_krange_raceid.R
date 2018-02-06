
#####################
# RACEID
#####################
# RaceID is an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. 
# The method is based on transcript counts obtained with unique molecular identifies.
# load libraries
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
source("method_resources/RaceID/RaceID_class.R")
source("skript/helper_files/Helper_functions.R")


# file paths

source("FILES.R")

# load data sets

data <- res.cluster <- vector("list", length(files))

names(data) <- names(res.cluster) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# load cell labels
labels <- load_labels(data) 


# RUN RaceID

par.maxexpr <-  list(
  kumar2015 = Inf,
  trapnell2014 = Inf,
  zhengmix2016 = Inf,
  koh2016 = Inf,
  simDataKumar=Inf,
  simDataKumar2=Inf
  
)


# define the number of cluster for kmeans clustering 
par.k <- list(
  kumar2015 = c(2:10),
  trapnell2014 = c(2:10),
  zhengmix2016=c(2:10),
  koh2016 = c(2:15),
  simDataKumar=c(2:10),
  simDataKumar2=c(2:10)
  
)


run_raceid <-  function(data,par.k,par.maxexpr,par.minexpr_percent){
# create store vectors
list<- vector("list", length(files))
names(list) <- names(files)
list->sc->res.cluster 
for (i in names(data)){
  
    sc[[i]] <- SCseq(as.data.frame(counts(data[[i]]))) # extract the expression data
    sc[[i]] <- filterdata(sc[[i]], mintotal = 1, minexpr = 0, 
                          minnumber = 0, maxexpr = par.maxexpr[[i]] , downsample = FALSE, dsn = 1, rseed = 1234)
    df.clus <- matrix( nrow = ncol(data[[i]]), ncol = length(par.k[[i]]) ) 
    
    for ( j in seq_len(length(par.k[[i]])) ) {
      df.clus[,j]<- clustexp(sc[[i]], metric="pearson", cln=par.k[[i]][j], do.gap=FALSE, clustnr=20, B.gap=50, SE.method="Tibs2001SEmax", 
                                SE.factor=.25, bootnr=50, rseed=1234)@kmeans$kpart
    }
    colnames(df.clus) <-  c( paste0(par.k[[i]]) )
    res.cluster[[i]] <- df.clus
}
return(res.cluster)

}

# run fuction

res.cluster <-  run_raceid(data,par.k, par.maxexpr, par.minexpr_percent)

# save clusters

dir_cluster <- paste0("results/filtered/raceid/raceid_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)


# save experiment labels

file_names <-  paste0("results/filtered/raceid/raceid_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/filtered/raceid/session_info_raceid_krange.txt")
sessionInfo()
sink()

### Appendix
