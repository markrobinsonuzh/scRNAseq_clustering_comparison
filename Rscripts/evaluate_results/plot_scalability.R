args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

methods <- strsplit(methods, ",")[[1]]

print(scalabilitydir)  ## Name of .rds file containing a SingleCellExperiment object
print(methods)  ## Clustering method
print(outrds)  ## Name of .rds file where results will be written

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## Load colors
source("Rscripts/Colorscheme.R") 

scalability <- do.call(rbind, lapply(methods, function(m) {
  readRDS(paste0(scalabilitydir, "/scalability_", m, ".rds")) %>%
    dplyr::select(dataset, method, nbrcells, run, elapsed)
}))
scalabilitysum <- scalability  %>%
  dplyr::group_by(dataset, method, nbrcells) %>%
  dplyr::summarize(elapsed = median(elapsed, na.rm = TRUE))

png(gsub("rds$", "png", outrds), width = 10, height = 8, unit = "in", res = 400)
ggplot(scalability %>% dplyr::filter(!is.na(elapsed)), 
       aes(x = nbrcells, y = elapsed, color = method)) +
  geom_smooth(data = scalabilitysum, se = FALSE) + 
  geom_point(size = 2) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, NA)) + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  facet_wrap(~ method, scales = "free") + 
  manual.scale + 
  xlab("Number of cells") + ylab("Elapsed time (s)")
dev.off()

saveRDS(NULL, outrds)
date()
sessionInfo()
