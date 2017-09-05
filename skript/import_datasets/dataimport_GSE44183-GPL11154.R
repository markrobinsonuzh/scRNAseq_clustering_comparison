###########################################################
# Load and process data set: Xue 2013 GSE44183-GPL11154
###########################################################

pdf("~/Desktop/masterthesis/results/QC_data/QC_GSE44183-GPL11154.pdf")

# load packages
source("~/Desktop/masterthesis/skript/helper_functions/Helper_functions.R")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
# load multiassay dataset from IMLS conquer
maex <- readRDS("~/Desktop/masterthesis/data/GSE44183-GPL11154.rds")
# summary of data set
print(maex)
dim(assay(maex)) # 65218 x 29


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
phn[, "phenoid" ] <- as.character(interaction(as.data.frame(phn[, c("source_name_ch1")])))
table(phn[, "phenoid"])


# Create an SCEobject 
sceset <- newSCESet(countData = cts, 
                   phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))


set_exprs(sceset, "tpm") <- tpms
### check matrices
get_exprs(sceset,"counts")[1:3,1:6] # count scaled length scaled tpms

get_exprs(sceset,"tpm")[1:3,1:6] # tpms

#################
# QC 
##################
#Plot overview for each cell
tryCatch({
  plot(sceset, block1 = paste("source_name_ch1", collapse = "."))
}, error = function(e) print(e))

# Calculate quality metrics and filter
# use ERCC spike inn as control features
if (sum(rowSums(counts(sceset)[grep("^ERCC-", rownames(sceset)), ] > 0) > 
        min(table(pData(sceset)[, paste("source_name_ch1", collapse = ".")]))) > 1) {
  withcontrols <- TRUE
  sceset <- calculateQCMetrics(sceset, feature_controls = grep("^ERCC-", rownames(sceset)))
} else {
  sceset <- calculateQCMetrics(sceset)
  withcontrols <- FALSE
}


#keep only those features that are observed in at least 1 cell.
keep_features <- rowSums(counts(sceset) > 0) >= 1

sceset <- sceset[keep_features, ]

# Plot QC metrics
plot_QC(sceset)


#### rename cell labels
# leave it for now

#### save cell labels

labels <- phenoData(sceset)$phenoid
dir_labels <-  paste0("~/Desktop/masterthesis/results/dataset_labels/labels_xue2013.txt")
write.table(labels , file=dir_labels, sep="\t")


###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "~/Desktop/masterthesis/data/sceset_GSE44183-GPL11154.rda")


###### save session info
sink(file = "~/Desktop/masterthesis/results/QC_data/session_info_GSE44183-GPL11154.txt")
sessionInfo()
sink()

dev.off()


## Appendix
### Automatic Cell filtering
sceset <-plotPCA(sceset,
                 size_by = "total_features", 
                 pca_data_input = "pdata",
                 detect_outliers = TRUE,
                 return_SCESet = TRUE)

