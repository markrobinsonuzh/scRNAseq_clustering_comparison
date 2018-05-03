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
# data: dataframe with variables: method, k , dataset, cell, trueclass
# method: character vector with method names for subsequent ensemble clustering
# Out: tibble longformat with method cluster,ensemble cluster, trueclass, dataset
# 
#------------------------------------------------------------------

# load files
df <- readRDS(file = "output/clustering_summary/clustering_summary.rds")
df.sub <- df%>%filter(!method %in% c("Seurat") )

helper_ensemble_trueclass <- function(methods, df){
  
  l <- vector("list",length(unique(df$dataset)) )
  names(l) <- unique(df$dataset)
  
  for(i in unique( df$dataset ) ){
    print(i)
    res <- df %>% filter( dataset %in% i ) %>% filter( k==length(unique(trueclass)) ) 
    # remove methods with all NAs
    isNA <- res%>%group_by(dataset, method)%>%summarise(isNA=all(is.na(cluster)) )%>%filter(isNA==TRUE)
    res <- res%>%filter(!method %in% isNA$method)
    # skip dataset if results for one method does not exist
    if ( sum( unique( res$method) %in% methods) != length(methods) )  {next}
    else {
    # wide format
    res.w <- dcast(res%>%filter( method %in% methods), trueclass+cell ~ method + run, 
                   value.var = c("cluster"))
    
    res2 <- res.w%>%select(-trueclass)%>% tibble::column_to_rownames('cell')%>%as.matrix # name data.frame
    # all NAs, replace by zero
     ifelse( all(is.na(res2)==TRUE),
             res2[] <- 0,
          res2 <- res2 )
    # compute conesnsus:
    m <- matrix(NA, nrow=nrow(res2),ncol = 5)
    
    for( j in 1:5 ){
      run <- grep(j, colnames(res2))
      res3 <- res2[, run  ] 
      re <- res3%>%plyr::alply(2, clue::as.cl_partition ) # partition
      re <- lapply(re, function(x) { x <- as.matrix(x$.Data) # some Nas, replace Nas by zeros
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
    colnames(m) <-  paste0(c(1:5))
 
    out <- cbind( dataset=rep(unique( res$dataset), nrow(m) )  ,m, as.matrix(res.w) ,method= rep( paste(methods, collapse = "."), nrow(m) )) %>%as.data.frame
    print(unique(out$method))
    l[[i]] <-out
  
    }
  }
  
  # collapse list
  res.df <-  plyr::rbind.fill(l)
  res.df <- tryCatch(reshape2::melt(data=res.df, id.vars=c("dataset","trueclass","cell","method"), measure.vars=c("1","2","3","4","5"),value.name=c("cons_cluster"), variable.name=c("run") ),
                     error=function(e) NULL)#
  return( res.df )
}



# compute ensembles, 2 methods:
# list of  ensemble combinations
comb.ensmbl2 <- combn( unique(df.sub$method), 2 , simplify = FALSE)
names(comb.ensmbl2) <- sapply(comb.ensmbl2, function(x) paste0(x, collapse = "."))
out2 <- lapply(comb.ensmbl2, helper_ensemble_trueclass, df=df.sub  )
out2 <- plyr::rbind.fill(out2)

saveRDS(out2, file= paste0("output/ensemble/clustering_ensemble_allmethods2.rds" ) )

# compute ensembles, 3 methods:
comb.ensmbl3 <- combn( unique(df.sub$method), 3 , simplify = FALSE)
names(comb.ensmbl3) <- sapply(comb.ensmbl3, function(x) paste0(x, collapse = "."))
out3 <- lapply(comb.ensmbl3, helper_ensemble_trueclass, df=df.sub  )
out3 <- plyr::rbind.fill(out3)
saveRDS(out3, file= paste0("output/ensemble/clustering_ensemble_allmethods3.rds" ) )


