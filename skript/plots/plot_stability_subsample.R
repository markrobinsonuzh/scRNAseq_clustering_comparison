#######################################################
# PLots for Stability analysis of methods using bootstrap
# 
#######################################################

# load librareis
library(plyr)
library(dplyr)
library(mclust)
library(cowplot)

# load helper files
source("skript/helper_files/Helper_functions.R")
# results directory
DATA_DIR <-  "results/stability_analysis"

files <- list(
  
  CIDR = file.path(DATA_DIR, "res.cidr.txt"),
  Linnorm = file.path(DATA_DIR, "res.linnorm.txt"),
  SIMLR = file.path(DATA_DIR, "res.simlr.txt"),
  pcaReduce = file.path(DATA_DIR,"res.pcareduce.txt"),
  Seurat = file.path(DATA_DIR,"res.seurat.txt"),
  TSCAN = file.path(DATA_DIR,"res.tscan.txt"),
  RaceID = file.path(DATA_DIR,"res.raceid.txt"),
  tSNEkmeans = file.path(DATA_DIR,"res.rtsnekmeans.txt"),
  ZINBWaVE = file.path(DATA_DIR,"res.zinbwave.txt")
  
  
)


# tidy up
tbl<-lapply(files, read_file  ) %>% ldply ( data.frame) 
names(tbl) <-  c("method", "ari")
# create different plots

p2 <-   ggplot(tbl)+
  geom_boxplot(aes(x=method,y=ari))+
  labs(x="method", y="ARI")



p4 <- ggplot(tbl, aes(x=method, y=ari, fill=method))+
  geom_dotplot(binaxis = "y", stackdir = "center",  dotsize = 0.5, stackratio=1)+
  labs(x="Method", y="ARI")+
  theme_gray()+
  guides(fill = "none") 

# save plot p4

save_plot(plot=p4,filename= "results/plots/stability_subsample_boot.pdf", base_width = 12)
