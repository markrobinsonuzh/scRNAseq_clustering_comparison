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
library(ggpubr)

## define the data directories
DATA_DIR <-  "results/run_results"
## read in ARI results from Rdata files
# files directories per dataset

files.ari.krange <- list(
  kumar2015 = file.path(DATA_DIR, "ari_krange_kumar2015.rda"),
  trapnell2014 = file.path(DATA_DIR, "ari_krange_trapnell2014.rda"),
  zhengmix2106 = file.path(DATA_DIR, "ari_krange_zhengmix2016.rda"),
  koh2016 = file.path(DATA_DIR, "ari_krange_koh2016.rda"),
  simDataKumar = file.path(DATA_DIR, "ari_krange_simDataKumar.rda"),
  simDataKumar2 = file.path(DATA_DIR, "ari_krange_simDataKumar2.rda")
  
)

# set intercept at position of ground truth for datasets
xintercept <- list(
  kumar2015 =3,
  trapnell2014 = 3,
  zhengmix2106 = 4,
  koh2016 = 10,
  simDataKumar = 4,
  simDataKumar2 = 4
  
)
# set limits for x axis for datasets
xlim<- list(
  kumar2015 =c(2,10),
  trapnell2014 = c(2,10),
  zhengmix2106 = c(2,10),
  koh2016 = c(2,15),
  simDataKumar = c(2,10),
  simDataKumar2 = c(2,10)
  
)
#________________________________________________________________________________________________________________________
# function to plot each method seperately
plot_ari_krange <- function(files.ari.krange, xintercept, xlim){
  # load files
  tmp <- lapply(files.ari.krange[[1]], function(x) get(load(x)))
  #remove the label column, sort to long format
  tmp <- ldply(tmp[[1]], as.data.frame)  
  tmp$par <- as.character(tmp$par)
  tmp$par <- as.numeric(tmp$par)
  tmp.k <- tmp%>%subset( .id %in% c("pcaReduce", "RtSNEkmeans",  "SIMLR", "cidr" ,  "tscan", "linnorm", "Seurat", "raceid", "zinbwave", "SC3"))  
  # plot the ARIs per dataset
  
  p <- ggplot(data = tmp.k, aes(x = ncluster, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    theme_bw()+
    facet_grid(.id~.)+
    guides(colour = "none")+
    labs(x="k")+
    geom_vline(xintercept = xintercept)+ 
    xlim(xlim)
                       
  return(p)
  
}
#________________________________________________________________________________________________________________________

# plot all the data, store in list
p.all <- mapply(plot_ari_krange,files.ari.krange,  xintercept,xlim, SIMPLIFY = FALSE)
# save plot per datafile
lapply(names(p.all), 
       function(x)ggsave(filename=paste0("results/plots/plot_ari_krange_ncluster_",x,".pdf"), plot=p.all[[x]]))
# in single plot

p.grid <- plot_grid(plotlist = p.all ,labels="auto" )
save_plot("results/plots/plot_ari_krange_ncluster_all.pdf", p.grid, base_height=15)
#________________________________________________________________________________________________________________________

### plot methods in one plot , per data set
plot_ari_krangemerged <- function(files.ari.krange){
  # load files
  tmp <- lapply(files.ari.krange[[1]], function(x) get(load(x)))
  #remove the label column, sort to long format
  tmp <- ldply(tmp[[1]], as.data.frame)  
  tmp$par <- as.character(tmp$par)
  tmp$par <- as.numeric(tmp$par)
  tmp.k <- tmp%>%subset( .id %in% c("pcaReduce", "RtSNEkmeans", "SC3", "SIMLR", "cidr" ,  "tscan", "linnorm", "Seurat", "raceid", "zinbwave"))  
  # plot the ARIs per dataset
  p <- ggplot(data = tmp.k, aes(x = ncluster, y = ARI, colour = .id))+       
    geom_line(aes(group=.id))+
    geom_point()+
    theme_bw()+
    guides(colour = "legend")+
    labs(x="k")
          
  return(p)
  
}
# plot list
p.allmerged <- lapply(files.ari.krange, plot_ari_krangemerged)
# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)

p.allmerged <- ggarrange(plotlist=p.allmerged , common.legend = TRUE, labels="auto")
# in single plot

save_plot("results/plots/plot_ari_krange_ncluster_allmerged.pdf",p.allmerged, base_height = 8)

