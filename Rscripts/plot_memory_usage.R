args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(memusetxt)
print(outrds)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

memuse <- read.delim(memusetxt, header = FALSE, sep = "") %>% 
  dplyr::rename(file = V1, used = V2, usedMb = V3, trigger = V4, 
                triggerMb = V5, maxuse = V6, maxuseMb = V7)
y <- memuse %>% tidyr::separate(file, into = c("method", "dtype"), sep = ":") %>% 
  tidyr::separate(method, into = c("tmp", "tmp2", "tmp3", "filtering",
                                   "dataset", "method"), sep = "_") %>% 
  dplyr::mutate(method = gsub("\\.Rout", "", method))

pdf(gsub("rds$", "pdf", outrds), width = 14, height = 10)
ggplot(y %>% dplyr::filter(dtype == "Vcells"), 
       aes(x = method, y = maxuseMb, fill = method)) + 
  geom_bar(stat = "identity") + 
  xlab("") + ylab("Max used Vcells (Mb)") + 
  facet_wrap(dataset ~ filtering, scales = "free_y") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
