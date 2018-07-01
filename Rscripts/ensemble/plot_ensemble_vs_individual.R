args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(ensemblerds)
print(clusteringsummary)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(ggthemes)
  library(viridis) 
  library(data.tree)
  library(ggtree)
})

## Read clustering results
res_ensmbl2 <- readRDS(file = ensemblerds)
res_summary <- readRDS(file = clusteringsummary)

## Initialize list to hold plots
plots <- list()

# ------------------------------------
# compute ARI, no of unique clusters
# ------------------------------------
sum_ensmbl2 <- res_ensmbl2 %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARIensmbl = mclust::adjustedRandIndex(cons_cluster, trueclass),
                   truenclust = length(unique(trueclass))) %>%
  tidyr::separate(method, sep = "[[:punct:]]", into = c("methone", "methtwo"), remove = FALSE) %>% 
  dplyr::ungroup()

sum_summary <- res_summary %>% dplyr::group_by(dataset, method, run, k) %>%
  dplyr::summarize(ARI = mclust::adjustedRandIndex(cluster, trueclass),
                   truenclust = length(unique(trueclass))) 
sum_ensmbl2$run <- as.integer(sum_ensmbl2$run)

# join datasets
sum_all <- left_join(sum_ensmbl2, sum_summary, 
                      by = c("dataset", "methone" = "method", "run", "k", "truenclust")) %>%
  dplyr::rename(method.ensmbl = method, ARIone = ARI) %>% 
  mutate(ARIdiff1 = ARIensmbl - ARIone)
sum_all2 <- left_join(sum_all, sum_summary, 
                      by = c("dataset", "methtwo" = "method", "run", "k", "truenclust")) %>% 
  dplyr::rename(ARItwo = ARI)

# compute median
res_median <- sum_all %>% dplyr::group_by(dataset, method.ensmbl, 
                                          k, truenclust, methone, methtwo) %>% 
  dplyr::summarize(ARIensmbl = median(ARIensmbl), 
                   ARIone = median(ARIone), ARIdiff1 = median(ARIdiff1))

plots[["ensembl_vs_first_allk"]] <- 
  ggplot(res_median, aes(ARIdiff1)) +
  geom_histogram() +
  facet_grid(methone ~ methtwo, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 2) +
  theme(axis.text.x = element_text(size = 6)) +
  labs(title = "Difference single method vs ensemble, all k",
       x = "Difference ARI")

plots[["ensembl_vs_first_truek"]] <- 
  ggplot(res_median %>% filter(k == truenclust), aes(ARIdiff1)) +
  facet_grid(methone ~ methtwo, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 2) +
  geom_histogram() +
  theme(axis.text.x = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  labs(title = "",
       x = "ARI difference between ensemble and method in row, for true k")

sum_max <- sum_all2 %>% select(dataset, methone, methtwo, k, ARIensmbl, run,
                               ARIone, ARItwo, truenclust) %>%
  group_by(dataset, methone, methtwo, k) %>% 
  summarise(med.ARIone = median(ARIone, na.rm = TRUE),
            med.ARItwo = median(ARItwo, na.rm = TRUE),
            med.ensemble = median(ARIensmbl, na.rm = TRUE),
            truenclust = unique(truenclust)) %>% ungroup()

# best method per ensemble
# diff
sum_max2 <- sum_max %>% dplyr::mutate(best = pmax(med.ARIone, med.ARItwo, na.rm = TRUE),
                                      worst = pmin(med.ARIone, med.ARItwo, na.rm = TRUE),
                                      diff.best = med.ensemble - best,
                                      diff.worst = med.ensemble - worst)
df.long <- reshape2::melt(sum_max2, id = c("dataset", "methone", "methtwo",
                                           "k", "truenclust", "med.ARIone", 
                                           "med.ARItwo", "med.ensemble", "best", "worst"))
df.long2 <- df.long %>% mutate(ARI = case_when(variable == "diff.best" ~ best,
                                               variable == "diff.worst" ~ worst))

# plot histo of differences in ARI, by the two methods 
plots[["ensembl_vs_bestworst_allk"]] <- 
  ggplot(df.long2 %>% filter() , aes(variable, value)) +
  geom_point(position = "jitter", alpha = 0.2, aes(color = ARI)) +
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  viridis::scale_color_viridis(option = "C", direction = -1) +
  labs(color = "ARI best/worst", x = "difference to ensemble", y = "ARI difference", 
       title = "difference between ensemble and best/worst method, for all k") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  scale_x_discrete(labels = c("diff.best" = "ensemble-best", "diff.worst" = "ensemble-worst"))

plots[["ensembl_vs_bestworst_truek"]] <- 
  ggplot(df.long2 %>% filter(k == truenclust), aes(variable, value)) +
  geom_point(position = "jitter", alpha = 0.2, aes(color = ARI)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  viridis::scale_color_viridis(option = "C", direction = -1) +
  labs(color = "ARI best/worst", x = "",
       y = "ARI difference between ensembl and best/worst individual method, for true k", 
       title = "") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_x_discrete(labels = c("diff.best" = "ensemble-best", "diff.worst" = "ensemble-worst"))

pdf(gsub("rds$", "pdf", outrds), width = 12, height = 12)
plots[["ensembl_vs_first_allk"]]
plots[["ensembl_vs_first_truek"]]
plots[["ensembl_vs_bestworst_allk"]]
plots[["ensembl_vs_bestworst_truek"]]
dev.off()

saveRDS(plots, file = outrds)

date()
sessionInfo()
