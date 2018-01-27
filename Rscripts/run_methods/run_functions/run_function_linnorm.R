
#####################
# Linnorm
#####################
# to do : include spikeinns

run_function_linnorm <-  function( data, labels, par.minNonZeroPortion, par.num_center ,par.BE_strength, datatype ){
  require(Linnorm)
  
  # create store vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  list->sys.time->transformedExp->res.cluster 
  
  for (i in names(data)){
    print(i)
    sys.time [[i]] <- system.time({
      transformedExp[[i]] <- Linnorm(counts(data[[i]]), spikein=NULL, minNonZeroPortion = par.minNonZeroPortion[[i]], BE_strength=par.BE_strength[[i]])
      res.cluster[[i]] <- Linnorm.tSNE(transformedExp[[i]], input = "Linnorm",num_center=par.num_center[[i]] )$k_means$cluster
    })
    
    
  }
  # save clusters
  
  dir_cluster <- paste0("results/",datatype,"/linnorm/linnorm_clus_", names(res.cluster), ".txt")
  
  
  save_clusters(res.cluster,dir_cluster)
  
  # save systemtime
  
  dir_systime <-  paste0("results/",datatype,"/linnorm/linnorm_systime_",names(sys.time),".txt")
  
  save_systemtime(sys.time, dir_systime)
  
  # save experiment labels
  
  dir_labels<-  paste0("results/",datatype,"/linnorm/linnorm_labels_",names(labels), ".txt")
  save_labels(labels, dir_labels )
  
  ###### Save Session Info
  sink(file = paste0("results/", datatype ,"/linnorm/session_info_linnorm.txt") )
  sessionInfo()
  sink()
  
  
}

