
#####################
# TSCAN
#####################
# TSCAN first preprocess the data by filter out rarley expressed cells and non DE features.
# The TSCAN::preprocess function additionally takes the log plus a pseudocount. function will rule out genes which
# have expression values of less than 1 in at least half of all cells. Genes with covariance smaller than one for the expression values are as well filterd out.
# Next pseudotime analysis is done, using dim reduction with PCA and model based clustering. A  pseudo temporal ordering score (POS) and travelling sales- man problem algorithm (TSP) algorithm combined with 
# Due to included preprocessing we use raw counts. Here the default transformation with log base2 , a pseudocount of one and cutoff for min expr of 2 
# is chosen. The method has an addiotionally criteria for the covariances for features, default is 1. We keep only high expressed genes, here we keep at least 5 pwercent of all genes for highly expressed cells.
set.seed(1234)
require(TSCAN)
require(dplyr)

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


# function for  TSCAN


#  filtered data
# define the minimum percentage of highly expressed cells (expression value bigger than minexpr_value) for the genes/features to be retained.
par.minexpr_percent <- list(
  kumar2015 = 0,
  trapnell2014 = 0,
  zhengmix2016 = 0,
  koh2016 = 0,
  simDataKumar = 0,
  simDataKumar2 = 0
  
)
# define k
par.k <- list(
  kumar2015 = 2:3,
  trapnell2014 = 2:6,
  zhengmix2016 = 2:10,
  koh2016= 2:15,
  simDataKumar=2:10,
  simDataKumar2=2:10
)
run_tscan <-  function(data, par.minexpr_percent, par.k ){
# create store vectors
list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->res.tscan->res.cluster 

for (i in names(data)){
  
    res.tscan[[i]] <- preprocess(counts(data[[i]]), minexpr_percent = par.minexpr_percent[[i]], logbase = 2) # preprocessing
    
    df.clus <- matrix( nrow = ncol(data[[i]]), ncol = length(par.k[[i]]) ) 
   
     for ( j in seq_len(length(par.k[[i]])) ) {
      print(i)
       print(paste0("k=", j))
      df.clus[,j] <- exprmclust( res.tscan[[i]], clusternum = par.k[[i]][j] )$clusterid # clustering
    }
    colnames(df.clus) <-  c( paste0(par.k[[i]]) )
    res.cluster[[i]] <- df.clus
    }
return(res.cluster)
}


# run tscan
res.cluster <-  run_tscan( data, par.minexpr_percent,par.k ) 


# save clusters

dir_cluster <- paste0("results/tscan/tscan_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)



# save experiment labels

file_names <-  paste0("results/tscan/tscan_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/tscan/session_info_tscan_krange.txt")
sessionInfo()
sink()

### Appendix
