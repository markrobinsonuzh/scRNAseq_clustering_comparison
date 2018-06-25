#------------------------------------------------------------------------------
# Figure 2 for Evaluation of clustering methods for single-cell RNA-seq data
#______________________________________________________________________________
#load libs
suppressPackageStartupMessages({
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(purrr)
  require(reshape2)
  library(ggthemes)
  library(viridis) 
})
# load colors
source("Rscripts/Colorscheme.R")  
# load summary results

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
#
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

# plot heat map on median stability with truenclust
#______________________________________________________
p1 <- ggplot( res_stab %>% filter(k==truenclust)%>%group_by(filtering, dataset, method, k)%>%
          dplyr::summarise( median.stability = median(ari.stab) ),
        aes(x = reorder( method,median.stability, FUN=mean, na.rm=TRUE  ),
            y = reorder(dataset,median.stability, FUN=mean, na.rm=TRUE  ),
            fill = median.stability ))+
  geom_tile(color="white", size=0.1)+
  facet_wrap(~ filtering)+
  scale_fill_viridis(name="ARI", direction=-1, na.value = "grey")+
  theme_tufte(base_family="Helvetica")+
  labs(x=NULL, y=NULL, title="") +
  coord_equal()+
  theme(axis.text.x=element_text(size=8, angle=90, hjust = 1, vjust = 1))+
  #theme(axis.text.y=element_text(size=15))+
  theme(panel.border=element_blank())+
  #theme(legend.title=element_text(size=15))+
  theme(legend.title.align=1)+
  #theme(legend.text=element_text(size=15))+
  theme(legend.position="right")+
  #theme(legend.key.size=unit(2, "cm"))+
  #theme(legend.key.width=unit(0.5, "cm"))+
  #theme(strip.text = element_text(size=16))+
  theme(axis.ticks=element_blank())


# ------------------------------------
# compute Entropy 
# ------------------------------------
shanon_entropy <- function(cluster){
  p <-c(table(cluster)) / length(cluster)
  s <- -1*sum(p*log2(p))
  return(s)
}

res_entropy <- res %>% 
  dplyr::group_by( dataset, method, run, k) %>% 
  dplyr::filter(!is.na(cluster)) %>% 
  dplyr::summarize(s = shanon_entropy(cluster),
                   s.true=  shanon_entropy(trueclass),
                   s.true.norm = s.true/log2(unique(k)),
                   s.norm=s/log2(unique(k)),
                   ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k)
  )  %>%
  tidyr::separate( dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
# plot entropy 
# ------------------------------------
p2 <- ggplot(data = res_entropy, 
              aes(x = reorder(method, s.norm, FUN=median, na.rm=TRUE),
                  y = s.norm, group=method, color=method))+       
         geom_boxplot()+  
         #facet_grid(dataset~., scale="free")+
         #geom_hline(aes(yintercept = s.true.norm), linetype = "dashed")+ 
         #geom_point(aes( x=truenclust , y=s.true ), color=1, shape=4)+
         manual.scale+
         theme_bw()+
         theme(axis.text.x = element_text(angle=90, size=8))+
         #theme(axis.text.y = element_text( size=15))+
         theme(legend.position = "none")+
         labs(x="method", y=expression("normalised entropy" *" "* frac(H,H[max]) ) ) 

# ------------------------------------
# Difference in k
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset,method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   timing = median(timing)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
# difference in k to maximum of median ARI (median by method and k)
diff_abs <- res_summary %>% 
  dplyr::group_by(dataset, filtering, method, truenclust, k) %>%
  dplyr::summarize(medARI=median(ARI)) %>%
  dplyr::filter(medARI==max(medARI)) %>%
  dplyr::mutate(k_diff= (k-truenclust))
# plot difference in k to max performance
# ------------------------------------
p3 <- ggplot(diff_abs, aes(x = method, y = k_diff, group = method, color = method)) + 
  geom_boxplot( outlier.color = NA, alpha=0.5) + 
  geom_dotplot( binaxis = "y", stackdir = "center",  dotsize = 0.2, stackratio=1)+
  theme_bw() +
  manual.scale+
  theme(axis.text.x = element_text(size=rel(1),angle = 90, hjust = 1, vjust = 1))+
  facet_grid(filtering~., scales="free")+
  theme(legend.position = "none")

# ---------------------------------------------------
#  normalized by median Rtsnekmeans , time combined
# ---------------------------------------------------

median.tsne <- res_summary%>%select(filtering, dataset, method, k, run, truenclust, timing)%>%
  dplyr::group_by( dataset, filtering , k)%>%
  filter(method=="RtsneKmeans")%>%dplyr::summarise(med.t=median(  timing ))%>%
  ungroup()
res.time <- res_summary %>% 
  group_by(filtering, dataset, method, k) %>% 
  dplyr::summarise(median.timing= median(timing) )%>%
  ungroup()

res.time <-full_join(res.time, median.tsne, by=c("dataset", "filtering", "k" )  ) %>% 
  dplyr::mutate(norm.time=median.timing/med.t)
p4 <- ggplot(res.time, 
             aes(x =reorder(method, norm.time, FUN=median, order=TRUE, na.rm=TRUE) , 
                 y = norm.time, group = method, color = method)) +
        manual.scale+
        geom_boxplot()+
        scale_y_log10()+
        labs(y= "log normalised time ", x="method")+
        theme_bw()+
        #labs(title="Runtime, normalized by the method Rtsnekmeans", size=16)+
        theme(axis.text.x=element_text( size=8, angle=90) )+
        theme(legend.position = "none")+
        #theme(axis.title = element_text(size=15))+
        #theme(legend.text = element_text(size=15))+
        theme(legend.position = "none")
        #theme(strip.text = element_text(size=16))
# ---------------------------------------------------
#  PLot grid Figure 2
# ---------------------------------------------------

row.one <- cowplot::plot_grid(p1, labels = c('a'))
row.two <- cowplot::plot_grid(p2,p3,p4, labels = c('b', 'c', 'd'), nrow = 1)
fig2 <-  cowplot::plot_grid(row.one, row.two , labels=c('', ''), ncol=1)


cowplot::ggsave("plots/figure2.pdf", fig2)




