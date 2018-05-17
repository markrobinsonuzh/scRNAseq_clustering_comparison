# ------------------------------------------------------------------------
# Stability analysis per method, by computing ARI for each run partition
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(purrr)
  require(reshape2)
})

res <-  readRDS(file="output/clustering_summary/clustering_summary.rds")
# ------------------------------------
# Compute the ARis
# ------------------------------------


res_summary<- res  %>% dplyr::group_by(dataset, method, k) %>% nest() 

# wide format
cast.x <-  function(x){
  d <- reshape2::dcast(x,cell~run, value.var ="cluster")
  return(d)
}

# 
res_nested <- res_summary  %>% mutate(data.wide  =  purrr::map( data, cast.x  )  ) 

# function compute ARI 
ARi_df <- function(x){
  stopifnot(class(x)=="data.frame")
  stopifnot(class(x[,1])=="character")
  
  x <- select(x, -cell)
  columns <- combn(ncol(x),2)
  ari2 <-array(NA,ncol(columns) )
  for (i in 1:10 ){
  ari2[i] <- mclust::adjustedRandIndex(  x[,columns[1,i]], x[,columns[2,i]]  )
  }
  return(ari2)
}
res_arrr <-res_nested  %>%  mutate(aririri  = purrr::map( data.wide, ARi_df   )  ) 


res_arrr$aririri


# ------------------------------------
# PLot the values
# ------------------------------------

# unnest the data

res_stab <- res_arrr %>% select(dataset ,  method ,k,aririri)%>%unnest()
res_stab <- res_stab%>% filter(method=="PCAKmeans", dataset=="sce_filteredExpr10_Kumar")

ggplot(res_stab , aes(k, aririri) ) +
  geom_point( )+
  facet_grid(dataset~method)

