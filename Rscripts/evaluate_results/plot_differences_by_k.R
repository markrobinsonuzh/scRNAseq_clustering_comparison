# differences between k at max ARI and k == truth
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
res <- readRDS(file="output/clustering_summary/clustering_summary.rds")
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


diff_abs3 <- res_summary %>% dplyr::group_by(dataset,filtering,method, truenclust,k) %>%summarize(medARI=median(ARI)) %>%filter(medARI==max(medARI)) %>%mutate(k_diff= (k-truenclust))


pdf("plots/performance/difference_in_k.pdf", height = 10, width=8)

print(ggplot(diff_abs3, aes(x = method, y = k_diff, group = method, color = method)) + 
        geom_boxplot( outlier.color = NA, alpha=0.5) + 
        geom_dotplot( binaxis = "y", stackdir = "center",  dotsize = 0.6, stackratio=1)+
        theme_bw() +
        scale_color_brewer(palette = "Set3")+
        labs(title = "difference_maxARI")+
        theme(axis.text.x = element_text(size=rel(1),angle = 90, hjust = 1, vjust = 1))+
        #geom_text(aes(label=dataset) ,size=3)+
        facet_wrap(~filtering, ncol=1)
        )
print(ggplot(diff_abs3, aes(x = method, y = k_diff, group = method, color = method)) + 
        geom_dotplot( binaxis = "y", stackdir = "center",  dotsize = 0.6, stackratio=1) + 
        theme_bw() +
        #scale_color_brewer(palette = "Set3")+
        labs(title = "difference_maxARI")+
        theme(axis.text.x = element_text(size=rel(1),angle = 90, hjust = 1, vjust = 1))+
        #geom_text(aes(label=dataset) ,size=2,  color=1,position=position_jitter(width=.2,height=.2) )+
        ggrepel::geom_text_repel(aes(label = dataset),
                        size = 2, color=1) +
        facet_wrap(~filtering, ncol=1)
)

# absolute values or not
diff_med <- res_summary %>% dplyr::group_by(dataset,filtering,method, k) %>%
  dplyr::summarize(medianARI= median(ARI))%>%ungroup()%>%dplyr::group_by(dataset,filtering,method, k)%>%
  dplyr::summarize(max=max(medianARI))


# estnclust
diff_abs <- res_summary %>% dplyr::group_by(dataset,filtering,method, run) 

# maximum ARI not important, look at difference of estnclust and truenclust
diff_estnclust <- res_summary %>%dplyr::group_by(dataset,filtering,method,estnclust,truenclust) %>%
  summarize(medARI=median(ARI)) %>%
  filter(medARI==max(medARI))%>%filter(estnclust != is.na(estnclust))%>%mutate(k_diff= (estnclust-truenclust))

diff_estnclust <- res_summary %>%  filter(estnclust != is.na(estnclust))%>%group_by(method, dataset, filtering, estnclust, truenclust)%>%
  summarise(k_diff= unique(estnclust)- unique(truenclust))
print(ggplot(diff_estnclust, aes(x = method, y = k_diff, group = method, color = method)) + 
        geom_boxplot( outlier.color = NA) + 
        geom_point(position = position_jitter(width = 0.1), alpha = 0.8)+
        theme_bw() +
        labs(title = "difference_estnclust")+
        scale_color_brewer(palette = "Dark2")+
        theme(axis.text.x = element_text(size=rel(1),angle = 90, hjust = 1, vjust = 1))+
        geom_text(aes(label=dataset), nudge_x = 0.0, nudge_y = 0.0, check_overlap = T,  size=3)+
        facet_wrap(~filtering, ncol=1))

dev.off()

