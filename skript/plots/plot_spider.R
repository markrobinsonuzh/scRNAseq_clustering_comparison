##################################################
### plot heatmap for differences  in ari scores ##
##################################################
# Takes the ari_"dataset".rda files , changes data to matrix and plots the differences with pheatmap.
# saves plots in resuls/plots directory

## load libraries
library(tibble)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(psych)
library(fmsb)

## define the data directories
DATA_DIR <-  "results/run_results"
datatype1<- "filtered"
datatype2 <- "default"
datatype3 <- "unfiltered"
datatype4 <- "optimalk"
datatype5 <- "smooth"
## read in cluster results from Rdata files
# files directories per dataset
files_ari_filtered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar.rda"))
)
files_ari_default <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_simDataKumar.rda"))
)
files_ari_unfiltered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar.rda"))
)
files_ari_optimalk <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar.rda"))
)
files_ari_smooth <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_simDataKumar.rda"))
)

#_______________________________________________________________________
# function to convert stored Ari results to matrix
ari2matrix <- function(files.ari) {
  tmp<- list(
    kumar2015 = NULL,
    trapnell2014 =NULL,
    koh2016 = NULL,
    zhengmix2016 =NULL,
    simDataKumar = NULL
  )
  tmp <- lapply(files.ari, function(x) get(load(x)))
  # create name vector
  names.methods <- lapply(tmp,names) %>% ldply( data.frame)
  names(names.methods) <- c("data", "method")
  # create table with data , remove the column with the "ground truth" (label)
  tmp <- unlist(tmp)
  tmp <- ldply(tmp, data.frame, .id=) 
  # combine data and names, remove the labels
  tmp <- cbind(tmp , names.methods)%>% subset(!(method=="labels") ) 
  # which methods are missing?
  print(table(tmp$method, tmp$data)==0)
  # table
  #xtabs <-  xtabs(X..i..~ data+ method, data=tmp)
  tbl <- acast(tmp, data~method, value.var="X..i..")
  # order
  tbl[, order(colnames(tbl))]
  
  return(tbl)
  
  
}
mat <- ari2matrix(files_ari_filtered )%>%as.data.frame%>%t%>%as.data.frame
class(mat)
radarchart(mat)
radarchart(mat[1:3,1:3], axistype=1, plty=1 )
?radarchart
# Example from https://github.com/ricardo-bion/ggradar
library(ggplot2)
library(ggradar)
suppressPackageStartupMessages(library(dplyr))
library(scales)
install.packages("ggradar")
mtcars %>%
  add_rownames( var = "group" ) %>%
  mutate_each(funs(rescale), -group) %>%
  tail(4) %>% select(1:10) -> mtcars_radar

ggradar(mtcars_radar) 
#_______________________________________________________________________
plot_ari_stars <- function(data, main) {
  require(GGally)
  
  # convert results
  ari.mat <- ari2matrix(data )
  p1 <- stars(t(ari.mat), key.labels = rownames(ari.mat),  key.loc = c(8,2), main=main )
  #p2 <- stars(ari.mat, labels = rownames(ari.mat),main=main)
}
pdf("results/plots/plot_ari_stars_all.pdf")
p1 <- plot_ari_stars(files_ari_filtered ,main="ARI filtered")
p2 <- plot_ari_stars(files_ari_unfiltered , main="rank ARI unfiltered")
p3 <- plot_ari_stars(files_ari_default , main="rank ARI default")
p4 <- plot_ari_stars(files_ari_optimalk , main="rank ARI optimalk")
p5 <- plot_ari_stars(files_ari_smooth , main="rank ARI smooth")
dev.off()
plot_ari_stars <- function(data, main) {
  require(GGally)
  
  # convert results
  ari.mat <- ari2matrix(data )
  #p1 <- stars(t(ari.mat), key.labels = rownames(ari.mat),  key.loc = c(8,2), main=main )
  p2 <- stars(ari.mat, key.labels = colnames(ari.mat),  key.loc = c(5,2), main=main)
}
pdf("results/plots/plot_ari_stars_perdata.pdf")
p1 <- plot_ari_stars(files_ari_filtered ,main="ARI filtered")
p2 <- plot_ari_stars(files_ari_unfiltered , main="rank ARI unfiltered")
p3 <- plot_ari_stars(files_ari_default , main="rank ARI default")
p4 <- plot_ari_stars(files_ari_optimalk , main="rank ARI optimalk")
p5 <- plot_ari_stars(files_ari_smooth , main="rank ARI smooth")
dev.off()

### Appendix

?spider
example(spider)
example(radar)
?
# Data must be given as the data frame, where the first cases show maximum.
maxmin <- data.frame(
  total=c(5, 1),
  phys=c(15, 3),
  psycho=c(3, 0),
  social=c(5, 1),
  env=c(5, 1) )
# data for radarchart function version 1 series, minimum value must be omitted from above.
RNGkind("Mersenne-Twister")
set.seed(123)
dat <- data.frame(
  total=runif(3, 1, 5),
  phys=rnorm(3, 10, 2),
  psycho=c(0.5, NA, 3),
  social=runif(3, 1, 5),
  env=c(5, 2.5, 4))
dat <- rbind(maxmin,dat)
op <- par(mar=c(1, 2, 2, 1),mfrow=c(2, 2))
radarchart(dat, axistype=1, seg=5, plty=1, vlabels=c("Total\nQOL", "Physical\naspects", 
                                                     "Phychological\naspects", "Social\naspects", "Environmental\naspects"), 
           title="(axis=1, 5 segments, with specified vlabels)", vlcex=0.5)

