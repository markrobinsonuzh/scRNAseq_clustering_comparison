########################################################
### plot radarcharts of ARI scores for the run modes ##
########################################################
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
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype1,"_simDataKumar2.rda"))
  
)
files_ari_default <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype2,"_simDataKumar2.rda"))
  
)
files_ari_unfiltered <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype3,"_simDataKumar2.rda"))
  
)
files_ari_optimalk <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype4,"_simDataKumar2.rda"))
  
)
files_ari_smooth <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_koh2016.rda") ),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_zhengmix2016.rda") ),
  simDataKumar = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR, paste0("ari_single_",datatype5,"_simDataKumar2.rda"))
  
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

#==================
# Plot 2: Same plot with custom features
plot_radarchart <- function(files, nobs, rownames, colnames, title){
  require(fmsb)
  # matrix
  mat <- ari2matrix(files )%>%as.data.frame%>%t%>%as.data.frame
  # add rows with max and min
  mat <- rbind(rep(1,5) , rep(0,5) , mat)
  # rename
  rownames(mat) <- rownames
  colnames(mat) <- colnames
  color = brewer.pal(nobs, "Set3")

  radarchart( mat , axistype=1 , maxmin=TRUE,
            #custom polygon
            pcol=color,  plwd=2 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,length.out=5), cglwd=0.8,
            #custom labels
            vlcex=1 ,
            centerzero=FALSE,
            title=title
  )
  #legend(x=0.9, y=0.5, legend = rownames(mat[-c(1,2),]), bty = "n", pch=20 , col=color , text.col = "black", cex=0.8, pt.cex=2)
}
#__________________
# matrix
mat <- ari2matrix(files_ari_filtered )%>%as.data.frame%>%t%>%as.data.frame
# add rows with max and min
mat <- rbind(rep(1,5) , rep(0,5) , mat)
# rename
row_names <- c("max","min","CIDR", "Linnorm"  , "pcaReduce", "RaceID", "tSNEkmeans", "SC3", "Seurat", "SIMLR", "SIMLRlarge", "TSCAN", "ZINBWaVE" )
col_names <- c("Koh", "Kumar", "simDataKumar", "simDataKumar2", "Trapnell", "Zheng")
# plot
# radarchart for several observation

#________________
pdf("results/plots/plot_radar_filtered.pdf")
plot_radarchart(files=files_ari_filtered , nobs=11, rownames=row_names, colnames = col_names, title="" )
dev.off()
pdf("results/plots/plot_radar_unfiltered.pdf")

plot_radarchart(files=files_ari_unfiltered , nobs=11, rownames=row_names, colnames = col_names, title="" )
dev.off()
pdf("results/plots/plot_radar_default.pdf")
plot_radarchart(files=files_ari_default , nobs=11, rownames=row_names, colnames = col_names, title="" )
dev.off()
pdf("results/plots/plot_radar_optimalk.pdf")
plot_radarchart(files=files_ari_optimalk , nobs=11, rownames=row_names, colnames = col_names, title="" )
dev.off()

pdf("results/plots/plot_radar_all.pdf")
par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(1,1,1,1),xpd=NA)
plot_radarchart(files=files_ari_filtered , nobs=11, rownames=row_names, colnames = col_names, title="filtered" )
plot_radarchart(files=files_ari_default , nobs=11, rownames=row_names, colnames = col_names,   title="default" )
plot_radarchart(files=files_ari_unfiltered , nobs=11, rownames=row_names, colnames = col_names,  title="unfiltered" )
color = brewer.pal(11, "Set3")
legend(x=2.5, y=1.0, legend = rownames(mat[-c(1,2),]), bty = "n", pch=20 , col=color , text.col = "black", cex=1, pt.cex=2,y.intersp=0.9)
dev.off()

