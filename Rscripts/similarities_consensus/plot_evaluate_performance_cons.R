
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
  library(viridis)
  library(ggthemes)
  
})
source("Colorscheme.R")
# ----------------------------------------------
# Plot summary of consensus clusterin, for control
# ----------------------------------------------

## Read clustering results
res <- readRDS(file="output/consensus/consensus_clue.rds")
pdf("plots/performance/res_performance_cons.pdf", width=20, height = 15)
# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset,method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(consensus.clue , trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k))%>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()
# --------------------------------------
# ## Calculate performance indices for each method and clustering run
# --------------------------------------
# remove the method Seurat in Zheng

print(ggplot(res_summary, 
        aes(x = k, y = ARI, group = method, color = method)) + 
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_line(size=1) + 
        theme_bw() +
        scale_color_brewer(palette = "Set3" )  +
        facet_grid(filtering ~ dataset, scales = "free_x")
)

# --------------------------------------
# ## Heatmap  ARI of truenclust, https://github.com/hrbrmstr/facetedcountryheatmaps
# --------------------------------------

# on true k , ARI, ordered by median
print(  p1 <- res_summary %>% dplyr::filter(k == truenclust) %>%
          dplyr::group_by(dataset, filtering, method, k) %>%
          ggplot(aes(x = reorder(method,ARI, median , na.rm=FALSE), y = dataset, fill = ARI))+
          geom_tile(color="white", size=0.1, na.rm =FALSE)+
          facet_wrap(~ filtering) +
          scale_fill_viridis(name="ARI", direction=-1 )+
          theme_tufte(base_family="Helvetica")+
          labs(x=NULL, y=NULL, title="consensus ARI, k = truenclust") +
          coord_equal() +
          theme(axis.text.x=element_text(size=12, angle=90))+
          theme(axis.text.y=element_text(size=12))+
          theme(legend.title=element_text(size=12))+
          theme(legend.title.align=1)+
          theme(legend.text=element_text(size=11))+
          theme(legend.position="right")+
          theme(legend.key.size=unit(2, "cm"))+
          theme(legend.key.width=unit(0.5, "cm"))+
          theme(axis.ticks=element_blank())
)
# on estimated k,  ARI

print(  p2 <- res_summary %>% dplyr::filter(k == estnclust) %>%
          dplyr::group_by(dataset, filtering, method, k) %>%
          ggplot(aes(x = reorder(method,ARI, median , na.rm=FALSE), y = dataset, fill = ARI))+
          geom_tile(color="white", size=0.1)+
          facet_wrap(~ filtering) +
          scale_fill_viridis(name="ARI", direction=-1, na.value = "grey")+
          theme_tufte(base_family="Helvetica")+
          labs(x=NULL, y=NULL, title="Consensus ARI, k = estnclust") +
          coord_equal() +
          theme(axis.text.x=element_text(size=12, angle=90))+
          theme(axis.text.y=element_text(size=12))+
          theme(panel.border=element_blank())+
          theme(legend.title=element_text(size=12))+
          theme(legend.title.align=1)+
          theme(legend.text=element_text(size=11))+
          theme(legend.position="right")+
          theme(legend.key.size=unit(2, "cm"))+
          theme(legend.key.width=unit(0.5, "cm"))+
          theme(axis.ticks=element_blank())
)

dev.off()
date()
sessionInfo()

