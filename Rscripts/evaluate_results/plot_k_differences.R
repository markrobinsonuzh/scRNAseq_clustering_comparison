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

## Initialize list to hold plots
plots <- list()

## Compute ARI, true number of clusters, estimated number of clusters, 
## elapsed time
res_summary <- res %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k),
                   elapsed = median(elapsed)) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

## Calculate the difference between the k that gives the maximal ARI and the true k
diff_abs <- res_summary %>% 
  dplyr::group_by(dataset, filtering, method, truenclust, k) %>%
  dplyr::summarize(medARI = median(ARI)) %>%
  dplyr::filter(medARI == max(medARI)) %>%
  dplyr::mutate(k_diff = (k - truenclust))

plots[["diff_kmax_ktrue"]] <- 
  ggplot(diff_abs, aes(x = method, y = k_diff, group = method, color = method)) + 
  geom_boxplot(outlier.color = NA, alpha = 0.5, size = 1.1) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.2, stackratio = 1) +
  theme_bw() +
  manual.scale +
  labs(title = "Difference between k giving maximal performance and true k", x="", y="Difference in k") +
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~ filtering, scales = "free") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16),
        legend.position = "none")

## Calculate the difference between the estimated and true number of clusters
diff_estnclust <- res_summary %>%  
  dplyr::filter(k == truenclust) %>%
  group_by(method, dataset, filtering, estnclust, truenclust) %>%
  dplyr::summarize(k_diff = unique(estnclust) - unique(truenclust)) %>%
  dplyr::ungroup()
diff_estnclust$method <- factor(diff_estnclust$method)

plots[["diff_kest_ktrue"]] <- 
  ggplot(na.omit(diff_estnclust), aes(x = method, y = k_diff, color = method)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.1, stackratio = 1) +
  theme_bw() +
  labs(title = "Difference between estimated and true k", x="", y="Difference in k") +
  manual.scale +
  facet_wrap(~filtering, scales="free_x")+
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, vjust = 1),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16),
        legend.position = "none")

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 7)
print(plots[["diff_kmax_ktrue"]])
print(plots[["diff_kest_ktrue"]])
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()
        

