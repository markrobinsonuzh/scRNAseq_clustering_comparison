################################
# Run all, single run scripts ##
################################
rm(list=ls())
# source the functions
FILE_DIR = "skript/run_methods/run_method/"

source( paste0(FILE_DIR , "meth_cidrv2.R") ) #ok
source( paste0(FILE_DIR , "meth_linnormv2.R") ) #ok
source( paste0(FILE_DIR , "meth_pcaReducev2.R") ) #ok
source( paste0(FILE_DIR , "meth_RtSNEkmeansv2.R") ) #ok
source( paste0(FILE_DIR , "meth_SIMLR_largescalev2.R") ) #ok
source( paste0(FILE_DIR , "meth_TSCANv2.R") ) #ok
source( paste0(FILE_DIR , "meth_RaceIDv2.R") ) #ok
source( paste0(FILE_DIR , "meth_seuratv2.R") ) #DLL fulll
source( paste0(FILE_DIR , "meth_SC3v2.R") ) # DLL fulll
# run the methods on the respetive datasets
# on defualt
datatype <- "default"
analyze_cidr(datatype)
analyze_linnorm(datatype)
analyze_pcareduce(datatype)
analyze_raceid(datatype)
analyze_rtsnekmeans(datatype) #error
analyze_seurat(datatype)
analyze_simlrlargescale(datatype)
analyze_tscan(datatype)
analyze_sc3(datatype)


# on filtered
datatype <- "filtered"
analyze_cidr(datatype)
analyze_linnorm(datatype)
analyze_pcareduce(datatype)
analyze_raceid(datatype)
analyze_rtsnekmeans(datatype) #error
analyze_seurat(datatype)
analyze_simlrlargescale(datatype)
analyze_tscan(datatype)
analyze_sc3(datatype)

# on unfiltered
datatype <- "unfiltered"
analyze_cidr(datatype)
analyze_linnorm(datatype)
analyze_pcareduce(datatype)
analyze_raceid(datatype)
analyze_rtsnekmeans(datatype) #error
analyze_seurat(datatype)
analyze_simlrlargescale(datatype)
analyze_tscan(datatype)
analyze_sc3(datatype)




