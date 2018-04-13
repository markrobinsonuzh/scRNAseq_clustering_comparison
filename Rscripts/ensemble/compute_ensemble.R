#--------------------------------------
# Ensemble between methods, per run 
#_______________________________________
suppressPackageStartupMessages({
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(clusterExperiment)
  require(clue)
  require(multidplyr)
  require(ggplot2)
  require(viridis)
  require(ggthemes)
  require(pheatmap)
  require(reshape2)
  require(mclust)
  require(RColorBrewer)
  
})
###  Helperfunction for computing ensemble clusters
# data: dataframe with variables method, k , dataset, cell, trueclass
# method: character vector with method names for ensemble
# Out: tibble longformat with method cluster,ensemble cluster, trueclass, dataset
#------------------------------------------------------------------

# load files
df <- readRDS(file = "../clustering_summary.rds")
df.sub <- df%>%filter(!dataset %in% grep("Zheng", unique(df$dataset) , value=TRUE) )

         
helper_ensemble <- function(methods, df){
  
  l <- vector("list",length(unique(df$dataset)) )
  names(l) <- unique(df$dataset)
  
  for(i in unique(df$dataset) ){
    print(i)
    
    res <- df %>% filter( dataset %in% i ) %>% filter(k==length(unique(trueclass))) 
    
    # wide format
    res.w <- dcast(res%>%filter(!method %in% c("Seurat"), method %in% methods), trueclass+cell ~ method + run, 
                   value.var = c("cluster"))
    
    res2 <- res.w%>%select(-trueclass)%>% tibble::column_to_rownames('cell')%>%as.matrix # name data.frame
    
    m <- matrix(NA, nrow=nrow(res2),ncol = 5)
    for(j in 1:5){
      run <- grep(j, colnames(res2))
      res3 <- res2[, run  ] 
      re <- res3%>%plyr::alply(2, clue::as.cl_partition ) # partition
      re <- lapply(re, function(x) { x <- as.matrix(x$.Data) # replace Nas by zeros
      x[is.na(x)] <- 0
      return(x)}
      ) 
      #lapply(re, function(x)apply(x$.Data, 2, count))# count cluster, doesnt work if zeros replaced
      #lapply(re, function(x)sum( apply(x$.Data, 2, count) ))# total sum
      re <- clue::as.cl_ensemble(re)
      re <- clue::cl_consensus(re,method = "SE",  control = list(nruns=50)) # no NAS!
      m[,j]<- clue::cl_class_ids(re)%>%as.matrix
      rownames(m) <- row.names(re$.Data)
    }
    colnames(m) <-  paste0("cons_run",c(1:5))
    out <- cbind( dataset=rep(unique( res$dataset), nrow(m) ) ,m, as.matrix(res.w) ) %>%as.data.frame
    l[[i]] <-out
    
  }
  res.df <-  plyr::rbind.fill(l)
  res.df <- reshape2::melt(data=res.df  , id.vars=c(1,7,8), measure.vars=c(2:6),value.name=c("cons_cluster") )

  return( res.df )
}
# which ensemble combinations

comb.ensmbl <- combn( unique(df$method),2 , simplify = FALSE)

names(comb.ensmbl) <- sapply(comb.ensmbl, function(x) paste0(x, collapse = ""))

ensembles <- list( RtsneKmeans_CIDR= c("CIDR","RtsneKmeans"),
                   SIMLR_CIDR= c("SIMLR", "RtsneKmeans"))


out <- lapply(ensembles, helper_ensemble, df=df.sub  )

# saVE

saveRDS(out, file= paste0("../clustering_ensemble.rds" ) )
