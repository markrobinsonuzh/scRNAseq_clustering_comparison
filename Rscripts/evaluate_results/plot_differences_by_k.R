args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(outrds)

# differences between k at max ARI and k == truth

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

# Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = clusteringsummary)

# ------------------------------------
# compute ARI, no of unique clusters, no of estimated k, median time
# ------------------------------------
res_summary <- res %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   elapsed = median(elapsed)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# difference in k to maximum of median ARI (median by method and k)
diff_abs <- res_summary %>% 
  dplyr::group_by(dataset, filtering, method, truenclust, k) %>%
  dplyr::summarize(medARI = median(ARI)) %>%
  dplyr::filter(medARI == max(medARI)) %>%
  dplyr::mutate(k_diff = (k - truenclust))

pdf(gsub("rds$", "pdf", outrds), width = 12)
p1 <- ggplot(diff_abs, aes(x = method, y = k_diff, group = method, color = method)) + 
        geom_boxplot(outlier.color = NA, alpha = 0.5) + 
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2, stackratio = 1) +
        theme_bw() +
        manual.scale +
        labs(title = "Difference to maximum performance") +
        theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, vjust = 1)) +
        facet_grid(~ filtering, scales = "free") +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              strip.text = element_text(size = 16))
print(p1)

p2 <- ggplot(diff_abs, aes(x = method, y = k_diff, group = method, color = method)) + 
  geom_dotplot(binaxis = "y", stackdir = "center",  dotsize = 0.1, stackratio = 1) + 
  theme_bw() +
  manual.scale +
  labs(title = "Difference to maximum performance", xlab = "method", ylab = "Difference in k") +
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, vjust = 1)) +
  ggrepel::geom_text_repel(aes(label = dataset),
                           size = 2, color = 1) +
  facet_wrap(~ filtering, ncol = 1) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16))
print(p2)

#  difference to estimated number of clusters
diff_estnclust <- res_summary %>%  
  dplyr::filter(k == truenclust) %>%
  group_by(method, dataset, filtering, estnclust, truenclust) %>%
  dplyr::summarize(k_diff = unique(estnclust) - unique(truenclust)) %>%
  dplyr::ungroup()
diff_estnclust$method <- factor(diff_estnclust$method)
p3 <- ggplot(na.omit(diff_estnclust), aes(x = method, y = k_diff, color = method)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.1, stackratio = 1) +
  theme_bw() +
  labs(title = "Difference of number of clusters to estimated number of clusters", 
       xlab = "method", ylab = "Difference in k") +
  manual.scale +
  facet_wrap(~filtering, scales="free_x")+
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, vjust = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16))
print(p3)
dev.off()

ggsave(filename = gsub("\\.rds$", "_2.pdf", outrds),
       cowplot::plot_grid(p1, p3, labels = "auto", ncol = 1)
)

saveRDS(NULL, outrds)
date()
sessionInfo()
        

