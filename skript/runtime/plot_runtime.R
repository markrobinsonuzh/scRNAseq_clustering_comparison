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
# load the helper function
source("skript/helper_files/Helper_functions.R")
# define which method to load
# methods: "pcaReduce","dbscan", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "SNNCliq", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm"
METHOD <- c("pcaReduce", "RtSNEkmeans", "SC3", "SIMLR","SIMLRlargescale", "cidr" , "Seurat", "zinbwave", "tscan","raceid", "linnorm")   

# file paths to the clustering results, change the path according to the processed datasets
DATA_DIR <-  "results"
DATASET <-"simDataKumar2"  # "kumar2015" ,"trapnell2014" ,"zhengmix2016" , "koh2016" , "simDataKumar","simDataKumar2"
datatype <- "filtered"

#-------------------------------------------------------------------------------------------------------------------
# store as single .rda object , per dataset
save_runtime_single(  METHOD,DATA_DIR, DATASET, datatype )

#-------------------------------------------------------------------------------------------------------------------
# function for loading the data
return_runtime_single <- function( METHOD,DATA_DIR, DATASET, datatype  ){
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
  time <- lapply(files_systime, read_file)%>% unlist%>% as.vector
  names(time) <- names(files_systime)
  return(time )
  # save the clusters as a rda file 
}

#-------------------------------------------------------------------------------------------------------------------

# save as vectors, for each data file: "kumar2015" ,"trapnell2014" ,"zhengmix2016" , "koh2016" , "simDataKumar"
kumar.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="kumar2015", datatype )
trapnell.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="trapnell2014", datatype )
zhengmix.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="zhengmix2016", datatype )
koh.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="koh2016", datatype )
simDataKumar.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="simDataKumar", datatype )
simDataKumar2.time <-  return_runtime_single(  METHOD,DATA_DIR, DATASET="simDataKumar2", datatype )

#tidy up
tbl <- cbind( kumar.time, trapnell.time, zhengmix.time, koh.time , simDataKumar.time,simDataKumar2.time )
tbl <- cbind( method= rownames(tbl), tbl) %>% as.data.frame

tbl<- melt(tbl, id.vars = c("method"), variable.name = "dataset", value.name = ".time")

tbl$.timelog <- log10(as.integer(tbl$.time))
tbl$.time <- as.integer(tbl$.time)
# rename levels of datasets
levels(tbl$dataset) <- c("Kumar", "Trapnell", "Zheng", "Koh", "simDataKumar", "simDatakumar2")
#-------------------------------------------------------------------------------------------------------------------
p1 <- ggplot(tbl)+
  geom_bar(aes(x=dataset,y=.timelog,fill=method),
           stat='identity',position='dodge')+
  labs(x="dataset", y="log10( runtime (s))")+
  scale_fill_brewer(palette = "Set3")
p1 <- ggplot(tbl)+
  geom_bar(aes(x=dataset,y=.time,fill=method),
           stat='identity',position='dodge')+
  labs(x="dataset", y="runtime (s)")+
  scale_fill_brewer(palette = "Set3")+
  scale_y_log10(breaks=c(0,10,100,1000),labels=c(0,10,100,1000))




p2 <- ggplot(tbl)+
  geom_bar(aes(x=dataset,y=.time,fill=method),
           stat='identity',position='dodge')+
  labs(x="data set", y="runtime (s)")+
  scale_fill_brewer(palette = "Set3")  


save_plot(plot=p1,filename= "results/plots/runtimeslog.pdf", base_width = 10)
save_plot(plot=p2,filename= "results/plots/runtimes.pdf", base_width = 10)



