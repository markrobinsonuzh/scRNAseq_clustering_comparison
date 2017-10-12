#######################################################
### Plot the ARI for a  range of clusters k 
### For all Methods and datasets
######################################################

# load libraries
library(dplyr)
library(plyr)
library(ggplot2)
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
tmp <- ldply(tmp[[1]], data.frame)  %>% filter(!(.id=="labels") )
# plot the ARIs per dataset
ggplot(data = tmp, aes(x = factor(par), y = ARI, colour = .id))+       
  geom_line(aes(group=.id))+
  geom_point()+
  facet_grid(.id~.)
return(tmp)
}
tmp <- lapply(files.ari, plot_ari)
tmp <- tmp%>%lapply(mutate(parameter=grepl("ˆ.")))

# plot the data
pdf("plot_ari_krange.pdf")


lapply(files.ari, plot_ari)
dev.off()
### Appendix
