args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(ggthemes)
  library(viridis) 
  
})

## Read clustering results
res <- readRDS(file="../../output/clustering_summary/clustering_summary.rds")
res2 <- readRDS(file="../../output/clustering_summary/clustering_summary_old.rds")

pdf("../../plots/performance/res_performance.pdf", width=15, height = 8)

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset,method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   timing = median(timing)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
unique(res_summary$filtering)
# --------------------------------------
# ## Calculate performance indices for each method and clustering run
# --------------------------------------
# remove the method Seurat in Zheng, as of high # of k

print(ggplot(res_summary%>%filter(!( method=="Seurat" & dataset %in% c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq") )), aes(x = k, y = ARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_smooth() + 
        theme_bw() +
        scale_color_brewer(palette = "Set3" )  +
        facet_grid(filtering ~ dataset, scales = "free_x")
)

print( ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
               dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust))%>%
                filter(!( method=="Seurat" & dataset %in% c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq") )),
             aes(x = k, y = medianARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_path() + 
        theme_bw() +
        scale_color_brewer(palette = "Set3" ) +
        facet_grid(filtering ~ dataset, scales = "free_x") 
)


# --------------------------------------
# ## Plot timing one boxplot per dataset, over all ks and runs as they are similar
# --------------------------------------

print(ggplot(res_summary, aes(x = method, y = timing, group = method, color = method)) + 
        geom_boxplot() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        scale_y_log10()+
        theme_bw() +
        scale_color_brewer(palette = "Set3")+
        theme(axis.text.x = element_text(size=rel(0.8),angle = 90, hjust = 1, vjust = 1)) )

# --------------------------------------
# ## Plot timing 
# --------------------------------------
# by k
print(ggplot(res_summary, aes(x = k, y = timing, group = method, color = method)) + 
        geom_smooth() + 
        facet_grid(filtering ~ dataset, scales = "free") + 
        scale_color_discrete(name = "") + 
        scale_color_brewer(palette = "Set3")+
        scale_y_log10())

# one line per method, normalized against maximum of Rtsnekmeans, max()
print(ggplot(res_summary %>% dplyr::group_by(dataset, filtering)%>%
             dplyr::mutate(normtime = (timing/ max(res_summary%>%
             filter(method=="RtsneKmeans")%>%select(timing)) ) ), 
             aes(x = k, y = normtime, group = method, color = method)) +
             scale_color_brewer(palette = "Set3")+
             geom_smooth()
)
# --------------------------------------
# ## Heatmap median ARI of truenclust, from https://github.com/hrbrmstr/facetedcountryheatmaps
# --------------------------------------
# on true k , median ARI, ordered by median
print(  res_summary %>% dplyr::filter(k == truenclust) %>%
          dplyr::group_by(dataset, filtering, method, k) %>%
          dplyr::summarize(medianARI = median(ARI)) %>%
          ggplot(aes(x = reorder(method,medianARI, median , na.rm=FALSE), y = dataset, fill = medianARI))+
          geom_tile(color="white", size=0.5, na.rm =FALSE)+
          facet_wrap(~ filtering) +
          scale_fill_viridis(name="medianARI", direction=-1 )+
          theme_tufte(base_family="Helvetica")+
          labs(x=NULL, y=NULL, title="median ARI, k = truenclust") +
          coord_equal() +
          theme(axis.text.x=element_text(size=8, angle=90))+
          theme(axis.text.y=element_text(size=8))+
          theme(legend.title=element_text(size=8))+
          theme(legend.title.align=1)+
          theme(legend.text=element_text(size=6))+
          theme(legend.position="right")+
          theme(legend.key.size=unit(2, "cm"))+
          theme(legend.key.width=unit(0.5, "cm"))+
          theme(axis.ticks=element_blank())
)
# on estimated k, median ARI

print(  res_summary %>% dplyr::filter(k == estnclust) %>%
          dplyr::group_by(dataset, filtering, method, k) %>%
        dplyr::summarize(medianARI = median(ARI)) %>%
        ggplot(aes(x = reorder(method,medianARI, median , na.rm=FALSE), y = dataset, fill = medianARI))+
        geom_tile(color="white", size=0.1)+
        facet_wrap(~ filtering) +
        scale_fill_viridis(name="medianARI", direction=-1, na.value = "grey")+
        theme_tufte(base_family="Helvetica")+
        labs(x=NULL, y=NULL, title="median ARI, k = estnclust") +
        coord_equal() +
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
)

dev.off()

date()
sessionInfo()


