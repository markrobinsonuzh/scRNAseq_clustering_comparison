# Foramt data
setwd("~/Desktop/masterthesis/data")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(scater))

# GSE41265
# import data
res <- readRDS("~/Desktop/masterthesis/data/GSE60749-GPL13112.rds")
# extract tpm , length-scaled count scaled tpm, pheno data

cts <- assays(experiments(res)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(res)[["gene"]])[["TPM"]]
phn <- pData(res)

# create scaterobject
res<- newSCESet(
  countData = cts, 
  tpmData = tpms,
  phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn))
)



# Quality control


