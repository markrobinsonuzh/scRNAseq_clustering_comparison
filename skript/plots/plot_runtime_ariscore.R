###############################################
### import and evaluate runtime results    ####
###############################################

# This file reads in the text files from the cluster runs for each method as well as the labels with the ground truth for the respective dataset. 
# The cluster and labels are then saved as a single list per dataset and stored as a "datasetname".rdata file in the run_results directory
# To do: if method not available, load labels from general file...
# load libs
require(ggplot2)
require(reshape2)
require(cowplot)
library(tibble)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
require(cowplot)


# load the helper function
source("skript/helper_files/Helper_functions.R")
# define which method to load
# methods: "pcaReduce","dbscan", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "SNNCliq", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm"
METHOD <- c("pcaReduce", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm")   

# file paths to the clustering results, change the path according to the processed datasets
DATA_DIR <-  "results"
DATASET <- c( "kumar2015" ,"trapnell2014" ,"zhengmix2016" , "koh2016" , "simDataKumar","simDataKumar2")
datatype <- "filtered"
#-__________________________________________________________________________________________________________________
# function for loading the data
return_runtime_single <- function( DATASET, METHOD,DATA_DIR, datatype  ){
  require(plyr)
  require(dplyr)
  ## files systime
  files_systime <- file.path(DATA_DIR,datatype, METHOD,paste0(METHOD,"_systime_",DATASET,".txt"))%>%as.list()
  # assign names
  names(files_systime) <- METHOD
  # load the .csv files
  time <- vector("list", length(files_systime))
  names(time) <- names(files_systime) 
  # exception handling
  read_file <- function(file){
    if ( file.exists(file)) {
      tryCatch( read.csv(file,sep="") , error=function(e) NA)
    } else {
      return(NA)
    }
  }
  # load stuff
  time <- lapply(files_systime, read_file)%>% unlist%>% as.vector%>%as.numeric
  names(time) <- names(files_systime)
  return(time )
  
}
#load the ARI files , as long table
## define the data directories
DATA_DIR <-  "results/run_results"
files_ari <- list(
  kumar2015 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_kumar2015.rda") ),
  trapnell2014 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_trapnell2014.rda") ),
  koh2016 = file.path(DATA_DIR, paste0("ari_single_",datatype,"_koh2016.rda")),
  zhengmix = file.path(DATA_DIR, paste0("ari_single_",datatype,"_zhengmix2016.rda")),
  simDataKumar = file.path(DATA_DIR,paste0("ari_single_",datatype,"_simDataKumar.rda")),
  simDataKumar2 = file.path(DATA_DIR,paste0("ari_single_",datatype,"_simDataKumar2.rda"))
  
)


# TABLE ARI scores
tbl_ari <- function(files_ari, datatype){

  
  # load dataset
  tmp <- lapply(files_ari, function(x) get(load(x)))
  # create name vector
  names.methods <- lapply(tmp,names) %>% ldply( data.frame)
  names(names.methods) <- c("data", "method")
  # create table with data , remove the column with the "ground truth" (label)
  tmp <- unlist(tmp)
  tmp <- ldply(tmp, data.frame, .id=) 
  # combine data and names, remove the labels
  tmp <- cbind(tmp , names.methods)%>% subset(!(method=="labels") ) 
  tmp$data <-  as.factor(tmp$data)
  levels(tmp$data) <- c("Koh", "Kumar", "simDataKumar", "simDataKumar2", "Trapnell", "Zheng")
  tmp$method <- as.factor(tmp$method)
  
  
  
  return(tmp)
}

#___________________________________________________________________________________________________________________

# load time data
t.tbl<- sapply(DATASET,  return_runtime_single, METHOD=METHOD,DATA_DIR="results",datatype="default")
# tidy up
t.tbl <- cbind( method= rownames(t.tbl),t.tbl) %>% as.data.frame
t.tbl<- melt(t.tbl, id.vars = c("method"), variable.name = "dataset", value.name = ".time")
# rename levels of datasets
levels(t.tbl$dataset) <- c("Kumar", "Trapnell", "Zheng", "Koh", "simDataKumar", "simDatakumar2")

# load ARI scores

ari.tbl <- tbl_ari(files_ari, "filtered")



# merge tables
m.tbl <- merge(ari.tbl,t.tbl, by.x=c("data", "method"), by.y=c("dataset", "method"), all.x = FALSE)
m.tbl$method <- factor(m.tbl$method)
m.tbl$.time <- as.numeric(m.tbl$.time)
levels(m.tbl$method) <- c("CIDR" , "Linnorm" ,"pcaReduce" ,"RaceID" ,"tSNEkmeans" ,"SC3"  ,   "Seurat"         
, "SIMLR" ,"SIMLRlargescale", "TSCAN","ZINBWaVE")
# averages of ARI and time
avg.tbl <- m.tbl%>%group_by(method)%>%summarise(avg.ari=mean(X..i..), avg.t = mean(.time))
p1 <- ggplot(avg.tbl, aes(x=avg.ari,y=avg.t))+
  geom_point(aes(color=method) , size=4)+
  labs(x="ARI", y="runtime (s)")+
  scale_color_brewer(palette = "Set3")+
  scale_y_log10()

p1 <- ggplot(avg.tbl, aes(x=avg.ari,y=avg.t))+
  geom_point(aes(color=method) , size=4)+
  labs(x="ARI", y="runtime (s)")+
  scale_color_brewer(palette = "Set3")+
  scale_y_continuous()+
  scale_x_continuous(limits=c(0,1))


#_____________________________________________________________
save_plot(plot=p1,filename= "results/plots/plot_avgaritime_filtered.pdf", base_width = 10)



