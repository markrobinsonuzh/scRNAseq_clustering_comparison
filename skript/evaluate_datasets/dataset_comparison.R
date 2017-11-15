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
  data[[i]] <- sceset
  
}

# compare datasets
comparison <- compareSCESets(data)
# 
# available plots
names(comparison$Plots)
comparison$Plots$Means
# create summary
panel <- makeCompPanel(comparison)
cowplot::save_plot("results/QC_data/comp_panel.png", panel, nrow = 4, ncol = 3)




# Appendix
panel <- makeOverallPanel(comparison, difference)
difference <- diffSCESets(data, ref = "Simple")
cowplot::save_plot("results/QC_Data/comparison_panel.png", panel, ncol = 4, nrow = 7)
