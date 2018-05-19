
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
res <- readRDS(file="output/ensemble/clustering_ensemble_allmethods2.rds")

pdf("plots/ensemble/performance_ensmbl2_truenclust.pdf", width=20, height = 15)

# ------------------------------------
# compute ARI, no of unique clusters
# ------------------------------------
res_summary <- res%>% dplyr::group_by(dataset, method, run) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# --------------------------------------
# ## Calculate performance indices for each method and clustering run
# --------------------------------------
# remove the method Seurat in Zheng, as of high # of k
print( ggplot(res_summary, aes(x = method, y = ARI)) + 
        geom_boxplot() + 
        theme_bw() +
        theme(axis.text.x=element_text(size=6, angle=90, hjust=1)) +
        scale_color_brewer(palette = "Set3" )  +
        facet_grid(filtering ~ dataset, scales = "free_x")
)

print( ggplot(res_summary, aes(x = method, y = ARI)) + 
        geom_point() + 
        theme_bw() +
        theme(axis.text.x=element_text(size=6, angle=90, hjust=1)) +
        scale_color_brewer(palette = "Set3" )  +
        facet_grid(filtering ~ dataset, scales = "free_x")
)
print( ggplot(res_summary, aes(x = method, y = ARI, color=filtering)) + 
         geom_point() + 
         theme_bw() +
         theme(axis.text.x=element_text(size=6, angle=90, hjust=1)) +
         scale_color_brewer(palette = "Set3" )  +
         facet_grid( ~ dataset, scales = "free_x")
)
# --------------------------------------
# ## Heatmap median ARI of truenclust, from https://github.com/hrbrmstr/facetedcountryheatmaps
# --------------------------------------
# on true k , median ARI, ordered by median
print(    res_summary  %>%
          dplyr::group_by(dataset, filtering, method) %>%
          dplyr::summarize(medianARI = median(ARI)) %>%
          ggplot(aes(x = reorder(method,medianARI, median , na.rm=FALSE), y = dataset, fill = medianARI))+
          geom_tile(color="white", size=0.6, na.rm =FALSE)+
          facet_wrap(~ filtering, nrow=4, ncol=1) +
          scale_fill_viridis(name="medianARI", direction=-1 )+
          theme_tufte(base_family="Helvetica")+
          labs(x=NULL, y=NULL, title="median ARI, k = truenclust") +
          coord_equal() +
          theme(axis.text.x=element_text(size=8, angle=90, hjust=1))+
          theme(axis.text.y=element_text(size=8))+
          theme(legend.title=element_text(size=8))+
          theme(legend.title.align=1)+
          theme(legend.text=element_text(size=6))+
          theme(legend.position="right")+
          theme(legend.key.size=unit(2, "cm"))+
          theme(legend.key.width=unit(0.5, "cm"))+
          theme(axis.ticks=element_blank())
)

dev.off()

date()
sessionInfo()


