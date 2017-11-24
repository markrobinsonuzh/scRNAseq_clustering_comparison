################################
# Run all, single run scripts ##
################################
rm(list=ls())

FILE_DIR = "skript/run_methods/run_method/"

source( paste0(FILE_DIR , "meth_cidr.R") ) #ok
source( paste0(FILE_DIR , "meth_dbscan.R") ) #ok
source( paste0(FILE_DIR , "meth_linnorm.R") ) #ok
source( paste0(FILE_DIR , "meth_pcaReduce.R") ) #ok
source( paste0(FILE_DIR , "meth_RtSNEkmeans.R") ) #ok
source( paste0(FILE_DIR , "meth_SIMLR_largescale.R") ) #ok
source( paste0(FILE_DIR , "meth_SIMLR.R") ) #
source( paste0(FILE_DIR , "meth_TSCAN.R") ) #ok
source( paste0(FILE_DIR , "meth_zinbwave.R") ) #ok
source( paste0(FILE_DIR , "meth_RaceID.R") ) #ok

source( paste0(FILE_DIR , "meth_SC3.R") ) # DLL fulll

source( paste0(FILE_DIR , "meth_seurat.R") ) #DLL fulll


