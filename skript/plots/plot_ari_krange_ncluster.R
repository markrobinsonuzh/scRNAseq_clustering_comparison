#######################################################
### Plot the ARI for a  range of clusters k 
### For all Methods and datasets
######################################################
# plots ARi k range results by parameters
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

files.ari.krange <- list(
  kumar2015 = file.path(DATA_DIR, "ari_krange_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "ari_krange_trapnell2014.rda"),
  zhengmix2106 = file.path(DATA_DIR, "ari_krange_zhengmix2016.rda"),
  koh2016 = file.path(DATA_DIR, "ari_krange_koh2016.rda"),
  simDataKumar = file.path(DATA_DIR, "ari_krange_simDataKumar.rda")
)
# function to plot 
plot_ari_krange <- function(files.ari.krange){
  # load files
  tmp <- lapply(files.ari.krange[[1]], function(x) get(load(x)))
  #remove the label column, sort to long format
  tmp <- ldply(tmp[[1]], as.data.frame)  
  tmp$par <- as.character(tmp$par)
  tmp$par <- as.numeric(tmp$par)
  tmp.k <- tmp%>%subset( .id %in% c("pcaReduce", "RtSNEkmeans", "SC3", "SIMLR", "cidr" ,  "tscan", "linnorm"))  
  # plot the ARIs per dataset
  
  p1 <- ggplot(data = tmp.k, aes(x = ncluster, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    facet_grid(.id~.)+
    guides(colour = "none")+
    labs(x="k")

  

  #pgrid <- plot_grid(p1,p2, p3, ncol=2)
  p <- grid.arrange(p1)                       # Number of rows
  return(p)
  
}

# plot all the data, store in list
p.all <- lapply(files.ari.krange, plot_ari_krange)

# save plot per datafile
lapply(names(p.all), 
       function(x)ggsave(filename=paste0("results/plots/plot_ari_krange_ncluster_",x,".pdf"), plot=p.all[[x]]))
# in single plot

p.grid <- plot_grid(plotlist = p.all ,labels="auto" )
save_plot("results/plots/plot_ari_krange_ncluster_all.pdf", p.grid, base_height=10)

### Appendix
