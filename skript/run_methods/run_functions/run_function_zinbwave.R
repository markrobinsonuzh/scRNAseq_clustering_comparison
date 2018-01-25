#####################
# ZINB-WaVE
#####################
# uses count  matrix. 
filter_hvg <- function(data, n.genes){
  filter <- rowSums( ( assay(data, "normcounts")  )>5)>5 
  data <- data[filter,]# filter out genes with at least five counts
  counts( data ) %>% log1p %>% rowVars -> vars
  names(vars) <- rownames( data )
  vars <- sort(vars, decreasing = TRUE)
  data <- data[names(vars)[1:n.genes],]
  return(data)
}

# run ZINBWave
run_function_zinbwave <-  function( data, labels, par.k,n.genes,datatype ){
  require(zinbwave)
  require(scater)
  require(Rtsne)
  require(dplyr)
  # create vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  list->sys.time->transformedExp-> res.zinb -> d -> tsne_data -> res.cluster
  for (i in names(data)){
    # filter for selction of hvg genes
    data[[i]] <- filter_hvg( data[[i]], n.genes[[i]] )
    # run ZINBWaVE
    sys.time [[i]] <- system.time({
    res.zinb[[i]] <- zinbFit( round(assay(data[[i]], "counts"),0) ,
                                K=2, epsilon=n.genes[[i]], verbose=TRUE,
                              nb.repeat.initialize = 2, 
                              maxiter.optimize = 25,
                              stop.epsilon.optimize = 1e-04)  # round data as it assumes whole counts

    d[[i]]<- dist(getW( res.zinb[[i]] ))
    res.cluster[[i]] <- kmeans(d[[i]], centers=par.k[[i]] )$cluster
    })
  }

  # save clusters
  
  dir_cluster <- paste0("results/",datatype,"/zinbwave/zinbwave_clus_", names(res.cluster), ".txt")
  
  
  save_clusters(res.cluster,dir_cluster)
  
  # save systemtime
  
  dir_systime <-  paste0("results/",datatype,"/zinbwave/zinbwave_systime_",names(sys.time),".txt")
  
  save_systemtime(sys.time, dir_systime)
  
  # save experiment labels
  
  dir_labels <-  paste0("results/",datatype,"/zinbwave/zinbwave_labels_",names(labels), ".txt")
  save_labels(labels, dir_labels )
  
  
  ###### Save Session Info
  sink(file = "results/",datatype,"/zinbwave/session_info_zinbwave.txt")
  sessionInfo()
  sink()
  

}
### Appendix

#commondispersion = TRUE, verbose = TRUE, 
