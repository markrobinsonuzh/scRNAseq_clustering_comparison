###########################################################
# Load and process data set: Kumar 2014 GSE60749-GPL13112
###########################################################
# load packages
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
# load multiassay dataset from IMLS conquer
maex <- readRDS("~/Desktop/masterthesis/data/GSE60749-GPL13112.rds")
# summary of data set
print(maex)
dim(assay(maex)) #45686 x 268


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
phn[, paste("source_name_ch1", collapse = ".")] <- as.character(interaction(as.data.frame(phn[, c("source_name_ch1","characteristics_ch1.1")])))
phn$phenoid <- phn$source_name_ch1
table(phn[, paste("phenoid", collapse = ".")])

# Create an SCEobject 
sceset <- newSCESet(countData = cts, 
                   phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
set_exprs(sceset, "tpm") <- tpms


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
plotQC(sceset, type = "highest-expression")

plotQC(sceset, type = "find-pcs", variable = "total_features")
plotQC(sceset, type = "find-pcs", variable = "total_features", plot_type = "pairs-pcs")

expl_vars <- c("source_name_ch1", "log10_total_counts", "log10_total_features", "pct_dropout",
               "pct_counts_top_200_features", ifelse(withcontrols, "log10_counts_feature_controls", NA),
               ifelse(withcontrols, "pct_counts_feature_controls", NA))

plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)])
plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)], method = "pairs")

plotPhenoData(sceset, aes_string(x = "log10_total_counts", y = "total_features",
                                 colour = paste("source_name_ch1", collapse = "."))) + 
  guides(colour = guide_legend( byrow = TRUE)) + theme(legend.position = "bottom")


if (withcontrols) {
  plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_feature_controls", 
                                   color = paste("source_name_ch1", collapse = "."))) + 
    guides(colour = guide_legend( byrow = TRUE)) + theme(legend.position = "bottom")
}


tryCatch({
  plotQC(sceset, type = "exprs-freq-vs-mean")
}, error = function(e) NULL)

### Automatic Cell filtering
sceset <-plotPCA(sceset,
                  size_by = "total_features", 
                  pca_data_input = "pdata",
                  detect_outliers = TRUE,
                  return_SCESet = TRUE)

# Dealing with confounders

###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "~/Desktop/masterthesis/data/sceset_GSE60749-GPL13112.rds")

###### save session info
sink(file = "~/Desktop/masterthesis/results/session_info_GSE60749-GPL13112.txt")
sessionInfo()
sink()
