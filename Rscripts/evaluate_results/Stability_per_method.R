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
# load colors
source("Rscripts/Colorscheme.R")  

res <-  readRDS(file="output/clustering_summary/clustering_summary.rds")
# ------------------------------------
# Compute the stability, based on ARI
# ------------------------------------
# nest df
res_summary <- res  %>% dplyr::group_by(dataset, method, k) %>% nest() 

res_summary<- res_summary%>%mutate(truenclust=purrr::map_int(data, function(x){
  y <- length(unique(x$trueclass))
  return(y)
}))
  
# wide format
cast.map <-  function(x){
  d <- reshape2::dcast(x,cell~run, value.var ="cluster")
  return(d)
}

res_nested <- res_summary  %>% mutate(data.wide  =  purrr::map( data, cast.map  )  ) 

# function for computing ARI 
ARi_df <- function(x){
  stopifnot(class(x)=="data.frame")
  stopifnot(class(x[,1])=="character")
  
  x <- select(x, -cell)
  columns <- combn(ncol(x),2)
  ari.nk <-array(NA,ncol(columns) )
  for (i in 1:10 ){
  ari.nk[i] <- mclust::adjustedRandIndex(  x[,columns[1,i]], x[,columns[2,i]]  )
  }
  stab <- as.data.frame( cbind(ari.stab=ari.nk, run1=columns[1,], run2=columns[2,]) )
  return(stab)
}
# compute ARI
res_stab.tmp <-res_nested  %>%  mutate(stability  = purrr::map( data.wide, ARi_df   )  ) 
# unnest
res_stab <- res_stab.tmp %>% select(dataset , method ,k,  stability, truenclust)%>%unnest()  %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) 
res_stab$k <- as.integer(res_stab$k)

# ------------------------------------
# PLot stability by k
# ------------------------------------

pdf("plots/performance/plot_stability.pdf", width=20)
# methods combined
ggplot( res_stab,
        aes(x = k, y = ari.stab, group = method, color = method))+ 
  geom_smooth() + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale+
  facet_grid(filtering~dataset, scales = "free_x" )+
  ylim(NA, 1)+
  labs(y="Stability (ARI)")
# methods seperated
ggplot( res_stab,
        aes(x = k, y = ari.stab, group = method, color = method))+ 
  geom_smooth() + 
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  theme_bw() +
  manual.scale+
  facet_grid(dataset+filtering~method, scales = "free_x" )+
  ylim(NA, 1)+
  labs(y="Stability (ARI)")
#___________________________________
# stability at truenclust
##________________________________
ggplot( res_stab %>% filter(k==truenclust),
        aes(x = method, y = ari.stab, group = method, color = method))+ 
  geom_boxplot() + 
  theme_bw() +
  manual.scale+
  facet_grid(filtering~dataset, scales = "free_x")+
  ylim(NA, 1)+
  labs(y="Stability (ARI)", title="k==truenclust")+
  theme( axis.text.x=element_text(size=10, angle=90))

#___________________________________
# plot heat map on median stability with truenclust
##________________________________
ggplot( res_stab %>% filter(k==truenclust)%>%group_by(filtering, dataset, method, k)%>%
          dplyr::summarise( median.stability = median(ari.stab) ),
          aes(x = reorder( method,median.stability, FUN=median, na.rm=TRUE  ),
              y = reorder(dataset,median.stability, FUN=median, na.rm=TRUE  ),
              fill = median.stability ))+
          geom_tile(color="white", size=0.1)+
          facet_wrap(~ filtering)+
          scale_fill_viridis(name="median stability (ARI)", direction=-1, na.value = "grey")+
          theme_tufte(base_family="Helvetica")+
          labs(x=NULL, y=NULL, title="median stability (ARI), k = truenclust") +
          coord_equal()+
          theme(axis.text.x=element_text(size=10, angle=90))+
          theme(axis.text.y=element_text(size=10))+
          theme(panel.border=element_blank())+
          theme(legend.title=element_text(size=12))+
          theme(legend.title.align=1)+
          theme(legend.text=element_text(size=10))+
          theme(legend.position="right")+
          theme(legend.key.size=unit(2, "cm"))+
          theme(legend.key.width=unit(0.5, "cm"))+
          theme(axis.ticks=element_blank()) 


dev.off()








