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
res <- readRDS(file="output/clustering_summary/clustering_summary.rds")
# consensus function 
# with clue package
#
cons.clue <- function(cluster, run, cell){
  require(dplyr)
  require(clue)
  
  tryCatch(m <-data.frame(cbind("cluster"=as.integer(cluster), "run"=as.integer(run), "id"=c(rep( 1:length(unique(cell)) , 5 ) ))) %>%
             tidyr::spread( key=run, value=cluster) %>%
             dplyr::select(c( -id)) %>%as.matrix, error=function(e) NA)
  ifelse( all(is.na(m)), 
          re <- rep(NA, 5*nrow(m) ) , 
          { re <- plyr::alply(m,2, clue::as.cl_partition )
          re <- lapply(re, function(x) { x <- as.matrix(x$.Data) # some Nas, replace Nas by zeros
          x[is.na(x)] <- 0
          return(x)}
          ) 
          re <- clue::as.cl_ensemble(re)
          re <- clue::cl_consensus(re,method = "SE",  control = list(nruns=50))
          re <- clue::cl_class_ids(re) 
          re <- rep(re,5)%>%as.character  }
  )
  return(re) 
}

# consensus function
# with clusterExperiment
cons.many <- function(cluster, run, cell){
  require(clusterExperiment)
  require(dplyr)
  m <-  tryCatch(as.data.frame(cbind("cluster"=as.integer(cluster), "run"=as.integer(run),  "id"=c(rep( 1:length(unique(cell)) , 5 ) )) )%>% # remove cbind as takes time
    tidyr::spread( run, cluster )%>% 
    dplyr::select( -id )%>%
    as.matrix, error=function(e)NA)
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
# are there any NAs?
all.na <- res%>%group_by(dataset, method) %>% dplyr::summarise(all.na = all(is.na(cluster) ) )%>% filter(all.na == TRUE) # remove RaceID, complete missing
somek.na <- res%>%group_by(dataset, method, k) %>% dplyr::summarise(all.na = all(is.na(cluster) ) )%>% filter(all.na == TRUE) # "FlowSOM" "RaceID"  "TSCAN"   "CIDR"  have some k with complete missing
na.cell <- res%>%group_by(dataset, method, k, run) %>% dplyr::summarise(na.cell = anyNA(cluster) ) %>% filter(na.cell == TRUE)  #  "SC3svm"  "SC3"     have some missing cells

# compute consensus,  RaceID and SIMLlarge empty, SC3 with NA in cluster, Seurat; group by resolution

  res_consmany <- res %>% filter(!is.na(cluster), !method%in% c( "Seurat", "SIMLRlargescale", "SC3svm", "SC3","RaceID") ) %>% 
    partition(dataset, method, k, cluster=cluster)%>%
    mutate( consensus.many=cons.many(cluster, run, cell))%>% collect()


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
# compute consensus,  Seurat; group by resolution

  res_consclue <- res %>% filter( !method %in% c( "Seurat")) %>%
   partition( dataset, method, k, cluster=cluster) %>%
   do( dplyr::mutate(.,consensus.clue = cons.clue(cluster, run, cell) )) %>%collect()
  
  # consensus for Seurat by resolution
  res_consclue.seurat <- res %>% filter( method %in% c( "Seurat")) %>%
    partition( dataset, method, resolution, cluster=cluster) %>%
    do( dplyr::mutate(.,consensus.clue = cons.clue(cluster, run, cell) )) %>%collect()
  
  
  res_consclue2 <- bind_rows(res_consclue2,res_consclue.seurat  )
  

saveRDS(res_consclue2, file="output/consensus/consensus_clue.rds")
# Unregister function 
cluster_rm(cluster, c("cons.clue"))
#_____________________________________________

############
cons.clue <- function(cluster, run, cell){
  require(dplyr)
  require(clue)
  
  tryCatch(m <-data.frame(cbind("cluster"=as.integer(cluster), "run"=as.integer(run), "id"=c(rep( 1:length(unique(cell)) , 5 ) ))) %>%
    tidyr::spread( key=run, value=cluster) %>%
    dplyr::select(c( -id)) %>%as.matrix, error=function(e) NA)
  ifelse( all(is.na(m)), 
        re <- rep(NA, 5*nrow(m) ) , 
        { re <- plyr::alply(m,2, clue::as.cl_partition )
        re <- lapply(re, function(x) { x <- as.matrix(x$.Data) # some Nas, replace Nas by zeros
        x[is.na(x)] <- 0
        return(x)}
        ) 
         re <- clue::as.cl_ensemble(re)
         re <- clue::cl_consensus(re,method = "SE",  control = list(nruns=50))
         re <- clue::cl_class_ids(re) 
         re <- rep(re,5)%>%as.character  }
        )
  return(re) 
}


