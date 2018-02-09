#####################
# tSNE + kmeans
#####################

# The following method uses the Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding for dimensionality reduction and kmeans foe clustering.
# Set Random seed for reproducible results.  Rtsne uses by  PCA for dimensionality reduction, by default 50 dimension are retained.
# In Rtsne the perplexity parameter has to be defined, typical values are between 5 and 50. Perplexity parameter is loosly speaking a guess about the number of close neighbors each point has. 
# Perplexity has to at least to be smaller than the number of points. 
# Other hyperparameters as the learning rate epsilon and the number of iteration which can give different results for different value ranges. 
# RUN tSNE and kmeans
# set random seed


run_function_snnclique <- function(data, labels,par.k,par.m,distan ,datatype){

  # create container vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  
  list->tinput_matrix->sys.time->res.rtsne->res.cluster 
  # extract transposed expression data
  for (i in names(data) ){
    tinput_matrix[[i]] <-  t(assay(data[[i]], "normcounts"))#use count scaled length scaled tpms, normalized and log2 transformed
  }
  
  # Run tSNE and kmeans
  for (i in names(tinput_matrix)){
    # construct a graph 
    scRNA.seq.funcs::SNN(
      data = tinput_matrix[[i]],
      outfile = "snn-cliq.txt",
      k = par.k[[i]],
      distance = distan
    )
    
    # find clusters in the graph
    
    sys.time[[i]] <- system.time({
      snn.res <- 
        system(
          paste0(
            "python method_resources/snn-cliq/Cliq.py ", 
            "-i snn-cliq.txt ",
            "-o res-snn-cliq.txt ",
            "-r ", par.r[[i]],
            " -m ", par.m[[i]]
          ),
          intern = TRUE
        )
    })
    
    #
    cat(paste(snn.res, collapse = "\n"))
    snn.res <- read.table("res-snn-cliq.txt")
    
    res.cluster<- as.character(snn.res[,1])
    # remove files that were created during the analysis
    system("rm snn-cliq.txt res-snn-cliq.txt")
  }
  # save clusters
  
  dir_cluster <- paste0("results/",datatype,"/snncliq/snncliq_clus_", names(res.cluster), ".txt")
  save_clusters(res.cluster,dir_cluster)
  
  # save systemtime
  
  dir_systime <-  paste0("results/",datatype,"/snncliq/snncliq_systime_",names(sys.time),".txt")
  save_systemtime(sys.time, dir_systime)
  
  # save experiment labels
  
  dir_labels <-  paste0("results/",datatype,"/snncliq/snncliq_labels_",names(labels), ".txt")
  save_labels(labels, dir_labels )
  
  
  ###### Save Session Info
  sink(file = "results/",datatype,"/snncliq/session_info_kmeans.txt")
  sessionInfo()
  sink()
}

