#--------------------------------------
# compute consensus, clue and clusterExperiment package
#_______________________________________
suppressPackageStartupMessages({
require(plyr)
require(dplyr)
require(tidyr)
require(clusterExperiment)
require(clue)
require(multidplyr)
require(ggplot2)
})
res <- readRDS(file="../../clustering_summary.rds")

# consensus function 
# with clue
#
cons.clue <- function(cluster, run, cell, dataset){
  require(dplyr)
  require(clue)
    m <- data.frame(cbind("cluster"=as.integer(cluster), "run"=as.integer(run), "id"=c(rep( 1:length(unique(cell)) , 5 ) ))) %>%
    tidyr::spread( key=run, value=cluster) %>%
    dplyr::select(c( -id)) %>%
    as.matrix
    re <- plyr::alply(m,2, clue::as.cl_partition )
    re <- clue::as.cl_ensemble(re)
    re <- clue::cl_consensus(re,method = "SE",  control = list(nruns=50))
    re <- clue::cl_class_ids(re) 
    re <- rep(re,5)%>%as.character
  return(re) 
}
# consensus function
# with clusterExperiment
cons.many <- function(cluster, run, cell){
  require(clusterExperiment)
  require(dplyr)
  m <-  as.data.frame(cbind("cluster"=as.integer(cluster), "run"=as.integer(run),  "id"=c(rep( 1:length(unique(cell)) , 5 ) )) )%>% # remove cbind as takes time
    tidyr::spread( run, cluster )%>% 
    dplyr::select( -id )%>%
    as.matrix
  re <- combineMany(m, proportion= 0.7, clusterFunction = "hierarchical01", propUnassigned = 0.5, minSize = 1)$cluster%>%
    rep(5)
  return(re)
}
#_____________________________________________
# Creating 2-core cluster 
cluster <- create_cluster(2)

# Register function and variable
cluster_assign_value(cluster, "cons.many", cons.many)

# Check registered items
cluster_ls(cluster)
# return items
cluster_get(cluster, "cons.many")

# compute consensus,  RaceID and SIMLlarge empty, SC3 with NA in cluster, Seurat; group by resolution
system.time(
  res_consmany <- res %>% filter(!is.na(cluster), !method%in% c( "Seurat", "SIMLRlargescale", "SC3svm", "SC3","RaceID") ) %>% 
    partition(dataset, method, k, cluster=cluster)%>%
    mutate( consensus.many=cons.many(cluster, run, cell))%>% collect()
)

saveRDS(res_consmany, file="consensus_clusterexp.rds")
# Unregister function 
cluster_rm(cluster, c("cons.many"))

#_____________________________________________
# Creating 2-core cluster 
cluster <- create_cluster(2)

# Register function and variable
cluster_assign_value(cluster, "cons.clue", cons.clue)

# Check registered items
cluster_ls(cluster)
# return items
cluster_get(cluster, "cons.clue")

# compute consensus, RaceID and SIMLlarge empty, SC3 with NA in cluster, Seurat; group by resolution

  res_consclue <- res2 %>% filter(!is.na(cluster), !method %in% c( "Seurat", "SIMLRlargescale", "SC3svm", "SC3","RaceID") ) %>% 
    partition( dataset, method, k, cluster=cluster) %>%
    mutate( consensus.clue=cons.clue(cluster, run, cell, dataset))%>% collect()


saveRDS(res_consclue, file="consensus_clue.rds")
# Unregister function 
cluster_rm(cluster, c("cons.clue"))
#_____________________________________________

