#------------------------------------
# Color scheme for plots
#------------------------------------

# Create a custom color scale
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(ggplot2)
})

## Color set, from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
colors <- c(CIDR = "#332288", FlowSOM = "#6699CC", PCAHC = "#88CCEE", 
            PCAKmeans = "#44AA99", pcaReduce = "#117733",
            RtsneKmeans = "#999933", Seurat = "#DDCC77", SC3svm = "#661100", 
            SC3 = "#CC6677", TSCAN = "grey34", ascend = "orange", SAFE = "black",
            monocle = "red", RaceID2 = "blue")

manual.scale <- scale_colour_manual(name = "", values = colors)
