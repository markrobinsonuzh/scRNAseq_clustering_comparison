#------------------------------------
# Color scheme for plots
#------------------------------------

# Create a custom color scale
library(RColorBrewer)
library(ggplot2)
# Methods
methods <- c("CIDR", "FlowSOM", "PCAHC", "PCAKmeans", "pcaReduce", "RtsneKmeans",
             "Seurat", "SC3svm","SC3" , "TSCAN", "ascend", "SAFE") 

## Color set, from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol12qualitative <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733",
                      "#999933", "#DDCC77", "#661100", "#CC6677",
                      "grey34",  "orange", "black")
# set 3
# colors <- brewer.pal(length(methods), "Set3")
colors <- tol12qualitative

# manual scale for ggplot
names(colors) <- methods
manual.scale <- scale_colour_manual(name = "method", values = colors)
