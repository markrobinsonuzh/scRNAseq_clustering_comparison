#------------------------------------
# Color scheme for plots
#------------------------------------

#Create a custom color scale
library(RColorBrewer)
library(ggplot2)
# methods, add additional methods here
methods <- c("CIDR" , "FlowSOM" , "PCAHC", "PCAKmeans","pcaReduce","RaceID", "RtsneKmeans",
             "SC3" ,    "SC3svm" , "Seurat"  ,  "SIMLR"  ,  "TSCAN"   ) 

# color set , from https://www.r-bloggers.com/the-paul-tol-21-color-salute/
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
#set 3
#colors <- brewer.pal(length(methods) ,"Set3")
colors <- tol12qualitative

# manual scale for ggplot
names(colors) <- levels(methods )
manual.scale <- scale_colour_manual(name = "method",values = colors)