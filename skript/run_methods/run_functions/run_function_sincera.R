#####################
# Sincera
#####################




run_function_sincera <- function(data, labels, filter.cell,datatype){
  require(SINCERA)
  # creae container vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  
  list->tinput_matrix->sys.time->res.rtsne->res.cluster 
  # extract transposed expression data
  for (i in names(data) ){
    tinput_matrix[[i]] <-  t(assay(data[[i]], "normcounts"))# use count scaled length scaled tpms, normalized and log2 transformed
  }
  # Run sincera
  for (i in names(data)){
    print(i)
    # construct object
    sc <- construct(exprmatrix=tinput_matrix[[i]],
                    samplevector=data[[i]]$phenoid)
    ifelse( filter.cell==TRUE, sc <- prefilterGenes(sc, pergroup=TRUE, min.expression=5, min.cells=2, min.samples=2), NULL)
    # Use expr.minimum function to set a minimum expression value
    sc <- expr.minimum(sc, value=0.01)
    # Perform per-sample z-score transformation
    sc <- normalization.zscore(sc, pergroup=FALSE)
    
    
  }
  
  
  # save clusters
  
  dir_cluster <- paste0("results/",datatype,"/RtSNEkmeans/RtSNEkmeans_clus_", names(res.cluster), ".txt")
  save_clusters(res.cluster,dir_cluster)
  
  # save systemtime
  
  dir_systime <-  paste0("results/",datatype,"/RtSNEkmeans/RtSNEkmeans_systime_",names(sys.time),".txt")
  save_systemtime(sys.time, dir_systime)
  
  # save experiment labels
  
  dir_labels <-  paste0("results/",datatype,"/RtSNEkmeans/RtSNEkmeans_labels_",names(labels), ".txt")
  save_labels(labels, dir_labels )
  
  
  ###### Save Session Info
  sink(file = "results/",datatype,"/RtSNEkmeans/session_info_kmeans.txt")
  sessionInfo()
  sink()
}

### Appendix

