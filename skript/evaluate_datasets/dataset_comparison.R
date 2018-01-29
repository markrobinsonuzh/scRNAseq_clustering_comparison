#################################################
# comparison between data sets
#################################################
# this script compares the diffrent datasets
# It foloows the splatter worflow

# load libraries
suppressPackageStartupMessages( library(MultiAssayExperiment) )

suppressPackageStartupMessages( library(splatter) )

# compare sets
# import data as sceset
# file paths

source("FILESraw.R")
# create vector

# load data sets

data <- vector("list", length(files))

names(data) <- names(files)
for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
names(data) <-  c("Kumar", "Trapnell", "Zheng", "Koh", "simDataKumar", "simDataKumar2")
# compare datasets
comparison <-splatter::compareSCEs(data)
# available plots
names(comparison$Plots)
# create summary
panel <- makeCompPanel(comparison, labels=c("a","b","c","d","e","f", "g", "h"), title="")
cowplot::save_plot("results/QC_data/comp_panel.png", panel, nrow = 4, ncol = 3)



# Appendix
#splatter::diffSCEs()
#panel <- makeOverallPanel(comparison, difference)
#difference <-diffSCEs(data, ref = "Simple")
#cowplot::save_plot("results/QC_Data/comparison_panel.png", panel, ncol = 4, nrow = 7)
