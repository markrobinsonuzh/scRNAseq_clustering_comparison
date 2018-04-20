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
res <- readRDS(file="output/ensemble/clustering_ensemble_allmethods2_by_k.rds")


pdf("plots/ensembles/res_performance_ensemble2_by_k.pdf", width=15, height = 8)

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
# --------------------------------------
# ## Calculate performance indices for each method and clustering run
# --------------------------------------
print(ggplot(res_summary, aes(x = k, y = ARI, group = method, color = method)) + 
        #geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_smooth() + 
        theme_bw() +
        facet_grid(filtering ~ dataset, scales = "free_x")
)


print( ggplot(res_summary %>% dplyr::group_by(dataset, filtering, method, k) %>%
                dplyr::summarize(medianARI = median(ARI), truenclust = unique(truenclust))%>%
                filter(!( method=="Seurat" & dataset %in% c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq") )),
              aes(x = k, y = medianARI, group = method, color = method)) + 
         #geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
         geom_path() + 
         theme_bw() +
         #scale_color_brewer(palette = "Set3" ) +
         facet_grid(filtering ~ dataset, scales = "free_x") 
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


dev.off()

date()
sessionInfo()



