#####################
# SNN-cliq
#####################



run_function_snncliq <- function(data, labels, par.k, par.m, par.r, distan , datatype){

  # create container vectors
  list<- vector("list", length(data))
  names(list) <- names(data)
  
  list->tinput_matrix->sys.time->snn.res->res.cluster 
  # extract transposed expression data
  for (i in names(data) ){
    tinput_matrix[[i]] <-  t(assay(data[[i]], "normcounts"))#use count scaled length scaled tpms, normalized and log2 transformed
  }
  
  # Run tSNE and kmeans
  for (i in names(tinput_matrix)){
    # construct a graph 
    scRNA.seq.funcs::SNN(
      data = tinput_matrix[[i]],
      outfile = "method_resources/snn-cliq/snn-cliq.txt",
      k = par.k[[i]],
      distance = distan
    )
    
    # find clusters in the graph
    
    sys.time[[i]] <- system.time({
      snn.res[[i]] <- 
        system(
          paste0(
            "python method_resources/snn-cliq/Cliq.py ", 
            "-i method_resources/snn-cliq/snn-cliq.txt ",
            "-o method_resources/snn-cliq/res-snn-cliq.txt ",
            "-r ", par.r[[i]],
            " -m ", par.m[[i]]
          ),
          intern = TRUE
        )
    })
    
    #
    cat(paste(snn.res[[i]], collapse = "\n"))
    snn.res[[i]] <- read.table("method_resources/snn-cliq/res-snn-cliq.txt")
    
    res.cluster[[i]]<- as.character(snn.res[[i]][,1])
    # remove files that were created during the analysis
    system("rm method_resources/snn-cliq/snn-cliq.txt method_resources/snn-cliq/res-snn-cliq.txt")
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

