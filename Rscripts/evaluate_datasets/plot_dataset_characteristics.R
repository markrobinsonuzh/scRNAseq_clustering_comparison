args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scefull)  
print(scefiltered)
print(outrds)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(scran)
  library(ggplot2)
  library(cowplot)
  library(tibble)
})

scefull <- readRDS(scefull)
scefiltered <- readRDS(scefiltered)

p3 <- plotQC(scefiltered, type = "highest-expression", n = 50)

df <- data.frame(ave.counts = log10(rowMeans(counts(scefull))), 
                 gene = rownames(scefull), 
                 row.names = rownames(scefull),
                 stringsAsFactors = FALSE) %>%
  dplyr::mutate(keep = ifelse(gene %in% rownames(scefiltered), "yes", "no"))
if (max(df$ave.counts[df$keep == "no"]) <= min(df$ave.counts[df$keep == "yes"])) {
  xint <- mean(min(df$ave.counts[df$keep == "yes"]), max(df$ave.counts[df$keep == "no"]))
  p4 <- ggplot(df, aes(x = ave.counts, group = keep)) + 
    geom_histogram(aes(alpha = keep), bins = 500) + 
    geom_vline(xintercept = xint, color = "red", linetype = "dashed") + 
    scale_alpha_manual(values = c(no = 0.3, yes = 0.9), name = "retain") + 
    xlab(expression(log[10]~"(average count)")) + ylab("Number of genes") + 
    theme(legend.position = "right")
} else {
  p4 <- ggplot(df, aes(x = ave.counts, group = keep)) + 
    geom_histogram(aes(alpha = keep), bins = 500) + 
    scale_alpha_manual(values = c(no = 0.3, yes = 0.9), name = "retain") + 
    xlab(expression(log[10]~"(average count)")) + ylab("Number of genes") + 
    theme(legend.position = "right")
}

expl_vars <- c("phenoid", "log10_total_counts", "log10_total_features", "pct_dropout",
               "pct_counts_top_200_features", "log10_counts_feature_controls",
               "pct_counts_feature_controls")
p5 <- plotQC(scefiltered, type = "explanatory-variables", variables = expl_vars)

p1 <- plotTSNE(scefiltered, colour_by = "phenoid")
p6 <- plotTSNE(scefiltered, colour_by = "total_features", size_by = "total_counts")

colData(scefiltered)$silhouette <- 
  cluster::silhouette(x = as.numeric(as.factor(scefiltered$phenoid)), 
                      dist = dist(prcomp(t(logcounts(scefiltered)))$x[, 1:10]))[, "sil_width"]
p2 <- ggplot(as.data.frame(colData(scefiltered)), 
             aes(x = phenoid, y = silhouette, color = phenoid)) + 
  geom_boxplot(outlier.size = -1) + geom_point(position = position_jitter(width = 0.2))

plot.qc <- plot_grid(p1, p6, p3, p4, p5, p2, ncol = 2, labels = "auto")
pdf(gsub("rds$", "pdf", outrds), width = 10, height = 12)
print(plot.qc)
dev.off()

saveRDS(NULL, outrds)

date()
sessionInfo()

  