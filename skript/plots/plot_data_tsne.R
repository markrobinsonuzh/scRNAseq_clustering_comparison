####################################################
# Plot datasets in tSNE space
####################################################
# load libraries
source("skript/helper_files/Helper_functions.R")

library(gridExtra)
library(ggfortify)
library(cowplot)
library(scater)
library(dplyr)
require(Rtsne)

# load data sets
# data directory
DATA_DIR = "data"

source("FILES.R")
data <- load_data(files, DATA_DIR)
# load cell labels
labels <- load_labels(data) 
labels <- lapply(labels, as.factor)
# plot
p <- lapply( data, function(x) {scater::plotTSNE(x, colour_by = "phenoid") }
            )
plot_grid(plotlist = p, labels="auto", )
p.grid <- plot_grid(plotlist=p,ncol = 2, nrow=3,labels="auto" )
save_plot("results/QC_data/plot_data_tsne.pdf", p.grid,
          base_width=12, base_height  =12)
### 
