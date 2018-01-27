#######################################################
# PLots for Stability analysis of methods using bootstrap
# 
#######################################################

# load librareis
library(dplyr)
library(mclust)
# load helper files
source("skript/helper_files/Helper_functions.R")
# results directory
DATA_DIR <-  "results/stability_analysis"

files <- list(
  
  cidr = file.path(DATA_DIR, "boot.cidr.txt"),
  linnorm = file.path(DATA_DIR, "boot.linnorm.txt"),
  simlr = file.path(DATA_DIR, "boot.simlr.txt"),
  pcareduce = file.path(DATA_DIR,"boot.pcareduce.txt"),
  seurat = file.path(DATA_DIR,"boot.seurat.txt"),
  tscan = file.path(DATA_DIR,"boot.tscan.txt"),
  raceid = file.path(DATA_DIR,"boot.raceid.txt")
  
)
# tidy up
tbl<-lapply(files, read_file  ) %>% ldply ( data.frame) 
names(tbl) <-  c("method", "ari", "NA")
# plot boxplot
p1 <- ggplot(tbl)+
  geom_point(aes(x=method,y=ari)+
  labs(x="method", y="ARI") )

p2 <-   ggplot(tbl)+
  geom_boxplot(aes(x=method,y=ari))+
  labs(x="method", y="ARI")



p4 <- ggplot(tbl, aes(x=method, y=ari, fill=method))+
  geom_dotplot(binaxis = "y", stackdir = "center",  dotsize = 0.6, stackratio=1)+
  labs(x="method", y="ARI")+
  theme_gray()+
  guides(fill = "none") 



save_plot(plot=p4,filename= "results/plots/stability_boot.pdf", base_width = 10)
