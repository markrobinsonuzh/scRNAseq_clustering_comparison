
#####################
# RACEID
#####################
# RaceID is an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. 
# The method is based on transcript counts obtained with unique molecular identifies.

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
# set the maxexpr parameter for UMI count data. discarding genes with at least maxexpr transcripts in at least a single cell after normalization or downsampling.
par.maxexpr <-  list(
  kumar2015 = Inf,
  trapnell2014 = Inf,
  zhengmix2016 = 500,
  koh2016 = Inf
)

# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
# Set a lower cutoff for the zhengmix data
par.minexpr_percent <- list(
  kumar2015 = 0.75,
  trapnell2014 = 0.75,
  zhengmix2016 = 0.1,
  koh2016 = 0.75
)

run_raceid <-  function(data,par.maxexpr,par.minexpr_percent){
# create store vectors
list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->sc->res.cluster 
for (i in names(data)){
  
    sc[[i]] <- SCseq(as.data.frame(counts(data[[i]]))) # extract the expression data
    sc[[i]] <- filterdata(sc[[i]], mintotal = 3000, minexpr = 5, 
                          minnumber = 1, maxexpr = par.maxexpr[[i]] , downsample = FALSE, dsn = 1, rseed = 1234)
    
    res.cluster[[i]]<- clustexp(sc[[i]], metric="pearson", cln=0, do.gap=TRUE, clustnr=20, B.gap=50, SE.method="Tibs2001SEmax", 
                                SE.factor=.25, bootnr=50, rseed=1234)@kmeans$kpart

  
}
}

# run fuction

res.cluster <-  run_raceid(data, par.maxexpr, par.minexpr_percent)
# save clusters

dir_cluster <- paste0("results/raceid/raceid_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/raceid/raceid_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/raceid/raceid_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/raceid/session_info_raceid_krange.txt")
sessionInfo()
sink()

### Appendix
