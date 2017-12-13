
#####################
# TSCAN
#####################
# TSCAN first preprocess the data by filter out rarley expressed cells and non DE features.
# The TSCAN::preprocess function additionally takes the log plus a pseudocount. function will rule out genes which
# have expression values of less than 1 in at least half of all cells. Genes with covariance smaller than one for the expression values are as well filterd out.
# Next pseudotime analysis is done, using dim reduction with PCA and model based clustering. A  pseudo temporal ordering score (POS) and travelling sales- man problem algorithm (TSP) algorithm combined with 
# Due to included preprocessing we use raw counts. Here the default transformation with log base2 , a pseudocount of one and cutoff for min expr of 2 
# is chosen. The method has an addiotionally criteria for the covariances for features, default is 1. We keep only high expressed genes, here we keep at least 5 percent of all genes for highly expressed cells.




# RUN TSCAN
run_function_tscan <-  function( data, labels ,par.minexpr_percent  ,par.clusternum , datatype ){
  
require(TSCAN)
require(dplyr)

# create store vectors
list<- vector("list", length(data))
names(list) <- names(data)
list->sys.time->res.tscan->res.cluster 
for (i in names(data)){
  print(i)
  sys.time [[i]] <- system.time({
    res.tscan[[i]] <- preprocess(counts(data[[i]]), minexpr_percent = par.minexpr_percent[[i]], logbase = 2) # preprocessing
    res.cluster[[i]] <- exprmclust(res.tscan[[i]], clusternum = 
                                     par.clusternum[[i]] )$clusterid # clustering
  })

  
}

# save clusters

dir_cluster <- paste0("results/",datatype,"/tscan/tscan_clus_", names(res.cluster), ".txt")
save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/tscan/tscan_systime_",names(sys.time),".txt")
save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels<-  paste0("results/",datatype,"/tscan/tscan_labels_",names(labels), ".txt")
save_labels(labels, dir_labels )


###### Save Session Info
sink(file = "results/",datatype,"/tscan/session_info_tscan.txt")
sessionInfo()
sink()
}
### appendix

