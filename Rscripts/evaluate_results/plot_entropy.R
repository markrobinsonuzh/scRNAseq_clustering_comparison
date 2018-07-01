args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(outrds)

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Load colors
source("Rscripts/Colorscheme.R") 

## Read clustering results
res <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

shannon_entropy <- function(cluster) {
  p <- c(table(cluster)) / length(cluster)
  s <- -1*sum(p*log2(p))
  return(s)
}

## Compute entropy and ARI
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
plots[["entropy_byds_byk"]] <- 
  ggplot(data = res_entropy %>% filter(!is.na(s)), 
         aes(x = k, y = s, group = method, color = method)) +       
  geom_smooth() +
  facet_grid(filtering ~ dataset, scale = "free_x") +
  geom_vline(aes(xintercept = truenclust), linetype = "dashed") + 
  geom_point(aes(x = truenclust, y = s.true), color = 1, shape = 4) +
  manual.scale +
  theme_bw() +
  labs(x = "k", y = "entropy") 

plots[["entropy_vs_ari_byds"]] <- 
  ggplot(data = res_entropy %>% filter(k == truenclust), 
             aes(x = ARI, y = s, group = method, color = method)) +       
        geom_point() +
        facet_grid(filtering ~ dataset, scale = "free") +
        geom_hline(aes(yintercept = s.true), linetype = "dashed") + 
        manual.scale +
        theme_bw() +
        labs(x = "ARI", y = "entropy") 

plots[["normentropy_by_ds"]] <- 
  ggplot(data = res_entropy %>% filter(k == truenclust), 
         aes(x = method, y = s.norm, group = method, color = method)) +       
  geom_boxplot() +
  facet_grid(. ~ dataset, scale = "free") +
  geom_hline(aes(yintercept = s.true.norm), linetype = "dashed") +
  manual.scale +
  theme_bw() +
  labs(x = "method", y = "norm.entropy", 
       title = "normalised entropy per dataset, k=truenclust") 

plots[["normentropy_allds_truek"]] <- 
  ggplot(data = res_entropy %>% filter(k == truenclust), 
         aes(x = reorder(method, s.norm, FUN = mean),
             y = s.norm, group = method, color = method)) +       
  geom_boxplot() +
  manual.scale +
  theme_bw() +
  labs(x = "method", y = "norm.entropy", 
       title = "normalised entropy all datasets, k=truenclust") 

plots[["normentropy_allds_allk"]] <- 
  ggplot(data = res_entropy, 
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

# Difference to truth at truenclust
plots[["deltaentropy_at_truth"]] <- 
  ggplot(data = res_entropy %>% filter(k == truenclust) %>%
           mutate(ds = s - s.true), 
         aes(x = method, y = ds, group = method, color = method)) +       
  geom_boxplot() +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  facet_grid(. ~ dataset, scale = "free") +
  manual.scale +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1),
        legend.position = "none") +
  labs(x = "method", y = "entropy", 
       title = "difference entropy at truth, k = truenclust")

# Difference to truth at truenclust, normalized entropy
plots[["deltanormentropy_at_truth"]] <- 
  ggplot(data = res_entropy %>% filter(k == truenclust) %>%
           mutate(ds.norm = s.norm - s.true.norm), 
         aes(x = method, y = ds.norm, group = method, color = method)) +       
  geom_boxplot(size = 1.1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  manual.scale +
  theme_bw() +
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = expression("Difference between normalised entropy" *" "* frac(H, H[max],) *" for clustering and truth")) 

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 6)
print(plots[["entropy_byds_byk"]])
print(plots[["entropy_vs_ari_byds"]])
print(plots[["normentropy_by_ds"]])
print(plots[["normentropy_allds_truek"]])
print(plots[["normentropy_allds_allk"]])
print(plots[["deltaentropy_at_truth"]])
print(plots[["deltanormentropy_at_truth"]])

dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()

