################################
# Run all, single run scripts ##
################################
rm(list=ls())
set.seed(1234)
# source the functions
FILE_DIR = "skript/run_methods/run_method/"


source( paste0(FILE_DIR , "meth_cidrv2.R") ) 
source( paste0(FILE_DIR , "meth_linnormv2.R") ) 
source( paste0(FILE_DIR , "meth_pcaReducev2.R") ) 
source( paste0(FILE_DIR , "meth_RtSNEkmeansv2.R") ) 
source( paste0(FILE_DIR , "meth_SIMLR_largescalev2.R") ) 
source( paste0(FILE_DIR , "meth_SIMLRv2.R") ) 
source( paste0(FILE_DIR , "meth_TSCANv2.R") ) 
source( paste0(FILE_DIR , "meth_RaceIDv2.R") ) 
source( paste0(FILE_DIR , "meth_zinbwavev2.R") ) 
source( paste0(FILE_DIR , "meth_seuratv2.R") ) 
source( paste0(FILE_DIR , "meth_SC3v2.R") ) 
# run the methods on the respetive datasets
# which dataset "kumar2015"  ,  "trapnell2014" ,"zhengmix2016" ,"koh2016"  ,"simDataKumar", "simDataKumar2"
dataset <- c( "kumar2015"  ,  "trapnell2014" ,"zhengmix2016" ,"koh2016"  ,"simDataKumar", "simDataKumar2"    )
# on defualt
datatype <- "default"

#run
analyze_cidr(datatype, dataset)
analyze_linnorm(datatype, dataset)
analyze_pcareduce(datatype,dataset)
analyze_rtsnekmeans(datatype, dataset) 
analyze_simlrlargescale(datatype, dataset)
analyze_simlr(datatype, dataset)
analyze_tscan(datatype, dataset)
analyze_raceid(datatype, dataset)
analyze_seurat(datatype, dataset)
analyze_zinbwave(datatype, dataset)
analyze_sc3(datatype, dataset)



# on filtered
datatype <- "filtered"
analyze_cidr(datatype, dataset)
analyze_linnorm(datatype, dataset)
analyze_pcareduce(datatype,dataset)  
analyze_rtsnekmeans(datatype, dataset) 
analyze_simlrlargescale(datatype, dataset)
analyze_simlr(datatype, dataset)
analyze_tscan(datatype, dataset)
analyze_raceid(datatype, dataset)
analyze_seurat(datatype, dataset)
analyze_zinbwave(datatype, dataset)
analyze_sc3(datatype, dataset)

# on unfiltered
datatype <- "unfiltered"
analyze_cidr(datatype, dataset)
analyze_linnorm(datatype, dataset)
analyze_pcareduce(datatype,dataset)
analyze_rtsnekmeans(datatype, dataset) 
analyze_simlrlargescale(datatype, dataset)
analyze_simlr(datatype, dataset)
analyze_tscan(datatype, dataset)
analyze_raceid(datatype, dataset)
analyze_seurat(datatype, dataset)
analyze_zinbwave(datatype, dataset)
analyze_sc3(datatype, dataset)

# on optimal k
datatype <- "optimalk"
analyze_cidr(datatype, dataset)
analyze_linnorm(datatype, dataset)
analyze_pcareduce(datatype,dataset)
analyze_rtsnekmeans(datatype, dataset) 
analyze_simlrlargescale(datatype, dataset)
analyze_simlr(datatype, dataset)
analyze_tscan(datatype, dataset)
analyze_raceid(datatype, dataset)
analyze_seurat(datatype, dataset)
analyze_zinbwave(datatype, dataset)
analyze_sc3(datatype, dataset)

# on smooth data
datatype <- "smooth"
analyze_cidr(datatype, dataset)
analyze_linnorm(datatype, dataset)
analyze_pcareduce(datatype,dataset)
analyze_rtsnekmeans(datatype, dataset) 
analyze_simlrlargescale(datatype, dataset)
analyze_simlr(datatype, dataset)
analyze_tscan(datatype, dataset)
analyze_raceid(datatype, dataset)
analyze_seurat(datatype, dataset)
analyze_zinbwave(datatype, dataset)
analyze_sc3(datatype, dataset)


