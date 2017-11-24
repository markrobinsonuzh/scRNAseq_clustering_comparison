#####################
# RACEID
#####################
# RaceID is an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. 
# The method is based on transcript counts obtained with unique molecular identifies.

source("skript/helper_files/Helper_functions.R")
# source file paths: fileterd , raw etc.
source("FILES.R")
# source method raceid
source("skript/run_methods/run_functions/run_function_raceid.R")
# load data sets

data <- load_data(files, DATA_DIR)

# load cell labels
labels <- load_labels(data) 

# RUN RaceID
# set the maxexpr parameter for UMI count data. discarding genes with at least maxexpr transcripts in at least a single cell after normalization or downsampling.
#the cln parameter is the number of cluster k in kmeans. Default is do.gap true and cln to 0, then gap stat is used to find k.
par.mintotal <- list(
  kumar2015 = 200,
  trapnell2014 = 200,
  zhengmix2016 = 200,
  koh2016 = 200,
  simDataKumar=200
)

par.maxexpr <-  list(
  kumar2015 = Inf,
  trapnell2014 = Inf,
  zhengmix2016 = Inf,
  koh2016 = Inf,
  simDataKumar=Inf
)
do.gap<-  list(
  kumar2015 = TRUE,
  trapnell2014 = TRUE,
  zhengmix2016 = TRUE,
  koh2016 = TRUE,
  simDataKumar=TRUE
)

par.cln <-  list(
  kumar2015 = 0,
  trapnell2014 = 0,
  zhengmix2016 = 0,
  koh2016 = 0,
  simDataKumar=0
)



# define datatype
datatype <- "default"

# run RACEID
run_function_raceid( data=data, labels=labels, datatype=datatype ,par.mintotal=par.mintotal, par.maxexpr=par.maxexpr,do.gap=do.gap,cln=par.cln) 


### Appendix



#plot clustering stats
plotgap(sc)
plotjaccard(sc)
plotsilhouette(sc)
# plot distances
clustheatmap(sc,final=FALSE,hmethod="single")
#findoutliers
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40), outdistquant=.75)
plotsensitivity(sc) # set threshold
plotbackground(sc) # check model fit
cdiff <- clustdiffgenes(sc,pvalue=.01) # not important
