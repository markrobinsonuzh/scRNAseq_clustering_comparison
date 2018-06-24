args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(outrds)

# ------------------------------------------
# Plots of entropy, by k, Seurat is excluded as f the number of k
# ------------------------------------------
suppressPackageStartupMessages({
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
})

# Load colors
source("Rscripts/Colorscheme.R") 

res <- readRDS(file = clusteringsummary)

shannon_entropy <- function(cluster) {
  p <- c(table(cluster)) / length(cluster)
  s <- -1*sum(p*log2(p))
  return(s)
}

# ------------------------------------
# compute Entropy and ARI
# ------------------------------------
res_entropy <- res %>% 
  dplyr::group_by(dataset, method, run, k) %>% 
  dplyr::filter(!is.na(cluster)) %>% 
  dplyr::summarize(s = shannon_entropy(cluster),
                   s.true = shannon_entropy(trueclass),
                   s.true.norm = s.true/log2(unique(k)),
                   s.norm = s/log2(unique(k)),
                   ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass)),
                   estnclust = unique(est_k)
                   ) %>%
  tidyr::separate(dataset, sep = "_", into = c("sce", "filtering", "dataset")) %>%
  dplyr::select(-sce) %>% dplyr::ungroup()

# plot entropy per k
# ------------------------------------
pdf(gsub("rds$", "pdf", outrds), width = 12, height = 6)

print(ggplot(data = res_entropy %>% filter(!is.na(s)), 
             aes(x = k, y = s, group = method, color = method)) +       
        geom_smooth() +
        facet_grid(filtering ~ dataset, scale = "free_x") +
        geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
        geom_point(aes(x = truenclust, y = s.true), color = 1, shape = 4) +
        manual.scale +
        theme_bw() +
        labs(x = "method", y = "entropy") 
)

print(ggplot(data = res_entropy %>% filter(k == truenclust), 
             aes(x = ARI, y = s, group = method, color = method)) +       
        geom_point() +
        facet_grid(filtering ~ dataset, scale = "free") +
        geom_hline(aes(yintercept = s.true), linetype = "dashed") + 
        manual.scale +
        theme_bw() +
        labs(x = "ARI", y = "entropy") 
)

print(ggplot(data = res_entropy %>% filter(k == truenclust), 
              aes(x = method, y = s.norm, group = method, color = method)) +       
         geom_boxplot() +
         facet_grid(. ~ dataset, scale = "free") +
         geom_hline(aes(yintercept = s.true.norm), linetype = "dashed") +
         manual.scale +
         theme_bw() +
         labs(x = "method", y = "norm.entropy", 
              title = "normalised entropy per dataset, k=truenclust") 
)

print(ggplot(data = res_entropy %>% filter(k == truenclust), 
              aes(x = reorder(method, s.norm, FUN = mean),
                              y = s.norm, group = method, color = method)) +       
         geom_boxplot() +
         manual.scale +
         theme_bw() +
         labs(x = "method", y = "norm.entropy", 
              title = "normalised entropy all datasets, k=truenclust") 
)

print(ggplot(data = res_entropy, 
             aes(x = k,
                 y = s.norm, group = method, color = method)) +
        geom_smooth() +
        manual.scale +
        theme_bw() +
        labs(x = "method", y = "norm.entropy", 
             title = "normalised entropy all datasets, by k") 
)

print(ggplot(data = res_entropy, 
             aes(x = reorder(method, s.norm, FUN = median, na.rm = TRUE),
                 y = s.norm, group = method, color = method)) +
         geom_boxplot() +
         manual.scale +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1)) +
         theme(axis.text.y = element_text(size = 15)) +
         theme(legend.position = "none") +
         labs(x = "method", y = expression("normalised entropy" *" "* frac(H, H[max])),
              title = "normalised entropy all datasets, all k") 
)

dev.off()

date()
sessionInfo()

