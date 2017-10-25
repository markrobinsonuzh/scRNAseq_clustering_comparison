#######################################################
### Plot the ARI for a  range of clusters k 
### For all Methods and datasets
######################################################

# load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library("gridExtra")

## define the data directories
DATA_DIR <-  "results/run_results"
## read in ARI results from Rdata files
# files directories per dataset

files.ari <- list(
  kumar2015 = file.path(DATA_DIR, "ari_krange_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "ari_krange_trapnell2014.rda"),
  xue2013 = file.path(DATA_DIR, "ari_krange_xue2013.rda"),
  koh2016 = file.path(DATA_DIR, "ari_krange_koh2016.rda")
)
# function to plot 
plot_ari <- function(x){
  # load files
  tmp <- lapply(x[[1]], function(x) get(load(x)))
  #remove the label column, sort to long format
  tmp <- ldply(tmp[[1]], as.data.frame)  
  tmp$par <- as.character(tmp$par)
  tmp$par <- as.numeric(tmp$par)
  tmp.k <- tmp%>%subset( .id %in% c("cidr","pcaReduce","RtSNEkmeans","SC3","SIMLR"))  
  
  tmp.eps <- subset( tmp, .id  %in% c("dbscan") )
  
  tmp.res <- subset( tmp, .id %in% c("Seurat") )
  # plot the ARIs per dataset
  
  p1 <- ggplot(data = tmp.k, aes(x = par, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    facet_grid(.id~.)+
    guides(colour = "none")+
    labs(x="k")
  
  
  p2 <- ggplot(data = tmp.eps, aes(x = par, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    labs(x="epsilon")
  
  p3 <- ggplot(data = tmp.res, aes(x = par, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    labs(x="resolution")
  
  #pgrid <- plot_grid(p1,p2, p3, ncol=2)
  p <- grid.arrange(p1,                             # First row with one plot spaning over 2 columns
                    arrangeGrob(p2,p3, nrow = 2), # Second row with 2 plots in 2 different columns
                    ncol = 2)                       # Number of rows
  return(p)
  
}
# plot all the data, store in list
p.all <- lapply(files.ari, plot_ari)
# save it to different files

lapply(names(p.all), 
       function(x)ggsave(filename=paste0("results/plots/plot_ari_krange_",x,".pdf"), plot=p.all[[x]]))

### Appendix
