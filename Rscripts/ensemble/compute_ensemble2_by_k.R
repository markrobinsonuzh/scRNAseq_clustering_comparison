#--------------------------------------
# Ensemble between methods, per run , all k
#_______________________________________
suppressPackageStartupMessages({
  
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(purrr)
  require(clue)
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
df <- readRDS(file = "output/clustering_summary/clustering_summary.rds")

# sample seurat 
df.sub <- df%>%filter( !method %in% c("Seurat") ) 
df.seurat <- df%>%filter(method=="Seurat")%>%
  group_by(dataset, method,k)%>%filter( resolution==sample(unique(resolution),1) )
df.sub <- plyr::rbind.fill(df.sub , df.seurat)

# helper for computing ensembles
helper_ensemble <- function(methods, df){
  
  l <- vector("list",length(unique(df$dataset)) )
  names(l) <- unique(df$dataset)
  for( i in unique( df$dataset ) ){
  print(i)
  res <- df %>% filter( dataset %in% i ) 
  # skip dataset if results for one method not exist
  combined_l <- vector("list", length(2:max(res$k))  )
      for ( u in 2:max(res$k) ) {
      res.k <- res %>% filter( k == u )
      # skip dataset if results for one method not exist
      if ( sum( unique( res.k$method) %in% methods) != length(methods)  ) {next}
      else {
      
      # wide format
      res.w <- dcast(res.k%>%filter(method %in% methods), trueclass+cell ~ method + run, 
                   value.var = c("cluster"))
      res2 <- res.w%>%select(-trueclass)%>% tibble::column_to_rownames('cell')%>%as.matrix # name data.frame
      # all NAs
      ifelse( all(is.na(res2)==TRUE),
            res2[] <- 0,
            res2 <- res2 )
    
      m <- matrix(NA, nrow=nrow(res2),ncol = 5)
       for( j in 1:5 ){
       run <- grep(j, colnames(res2))
       res3 <- res2[, run  ] 
       re <- res3%>%plyr::alply(2, clue::as.cl_partition ) # partition
       re <- lapply(re, function(x) { x <- as.matrix(x$.Data) # some Nas, replace Nas by zeros
       x[is.na(x)] <- 0
       return(x)}) 
       #lapply(re, function(x)apply(x$.Data, 2, count))# count cluster, doesnt work if zeros replaced
       #lapply(re, function(x)sum( apply(x$.Data, 2, count) ))# total sum
       re <- clue::as.cl_ensemble(re)
       re <- clue::cl_consensus(re,method = "SE",  control = list(nruns=50)) # no NAS!
       m[,j]<- clue::cl_class_ids(re)%>%as.matrix
       rownames(m) <- row.names(re$.Data)
       }
      colnames(m) <-  paste0(c(1:5))
    
    out <- cbind( dataset=rep(unique( res$dataset), nrow(m) ), m , as.matrix(res.w) ,method= rep( paste(methods, collapse = "."), nrow(m) ), k= rep(u, nrow(m)) ) %>%as.data.frame
    print(unique(out$method))
    combined_l[[u]]<- out
      }
      }
  
  l[[i]] <- plyr::rbind.fill(combined_l)
    
  }
  res.df <-  plyr::rbind.fill(l)
  res.df <- reshape2::melt(data=res.df  , id.vars=c("dataset","trueclass","cell","method", "k"), measure.vars=c("1","2","3","4","5"),value.name=c("cons_cluster"), variable.name=c("run") )
  return( res.df )
}
# which ensemble combinations
#comb.ensmbl <- combn( unique(df.sub$method), 2 , simplify = FALSE)
comb.ensmbl <- list(l1=unique(df.sub$method),l2=unique(df.sub$method) )%>%cross()%>%purrr::map((paste))

names(comb.ensmbl) <- sapply(comb.ensmbl, function(x) paste0(x, collapse = ".") )
# remove identical ensembles
comb.ensmbl <- comb.ensmbl %>% discard(function(x) x[1] == x[2] )
# run 
out <- lapply( comb.ensmbl , helper_ensemble, df=df.sub  )
out <- plyr::rbind.fill(out)
# saVE
saveRDS(out, file= paste0("output/ensemble/clustering_ensemble_allmethods2_by_k.rds" ) )
