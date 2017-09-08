###########################################################
# Load and process data set: Koh 2016 SRP073808
###########################################################
pdf("results/QC_data/QC_SRP073808.pdf")


# load packages
source("skript/helper_files/Helper_functions.R")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
# load multiassay dataset from IMLS conquer
maex <- readRDS("data/SRP073808.rds")
# summary of data set

print(maex)
dim(assay(maex)) # 65218x651


# extract the gene-level TPM and the counts obtaines by scaling the TPMs

if (all(c("count_lstpm", "TPM") %in% names(assays(experiments(maex)[["gene"]])))) {
  cts <- assays(experiments(maex)[["gene"]])[["count_lstpm"]]
  tpms <- assays(experiments(maex)[["gene"]])[["TPM"]]
} else {
  cts <- assays(experiments(maex)[["gene"]])[["count"]]
  tpms <- NULL
}

# extract the phenotype data

phn <- colData(maex)
# extract the cell annotation to phenoid col

phn[, "phenoid"] <- as.character(phn[,"LibraryName"])


table(phn[, "phenoid"])


# Create an SCEobject 
sceset <- newSCESet(countData = cts, 
                    phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
set_exprs(sceset, "tpm") <- tpms


#################
# QC 
##################
#Plot overview for each cell

plot(sceset, block1 = paste("phenoid"))


# Calculate quality metrics and filter
# use ERCC spike inn as control features
if (sum(rowSums(counts(sceset)[grep("^ERCC-", rownames(sceset)), ] > 0) > 
        min(table(pData(sceset)[, paste("phenoid", collapse = ".")]))) > 1) {
  withcontrols <- TRUE
  sceset <- calculateQCMetrics(sceset, feature_controls = grep("^ERCC-", rownames(sceset)))
} else {
  sceset <- calculateQCMetrics(sceset)
  withcontrols <- FALSE
}


#keep only those features that are observed in at least 1 cell.
keep_features <- rowSums(counts(sceset) > 0) >= 1



sceset <- sceset[keep_features, ]
dim(sceset) # 48986   x   651 
# Plot QC metrics






###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "data/sceset_SRP073808.rda")

###### save session info
sink(file = "results/QC_data/session_info_SRP073808.txt")
sessionInfo()
sink()

dev.off()


## Appendix
### Automatic Cell filtering
sceset <-plotPCA(sceset,
                 size_by = "total_features", 
                 pca_data_input = "pdata",
                 detect_outliers = TRUE,
                 return_SCESet = TRUE)

# Dealing with confounders
