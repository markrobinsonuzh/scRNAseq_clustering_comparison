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
print(plotQC(sceset, type = "highest-expression")) # Top 50 genes account for 18.4 % of the counts, profile relatively flat

print(plotQC(sceset, type = "find-pcs", variable = "total_features")) # PC1 shows high correlation with total number observerd features per cell
print(plotQC(sceset, type = "find-pcs", variable = "total_features", plot_type = "pairs-pcs")) # PC1 and PC75, no correlation between other PCs


# plot the variance explained for some explanatory vars.
expl_vars <- c("phenoid", "log10_total_counts", "log10_total_features", "pct_dropout",
               "pct_counts_top_200_features", ifelse(withcontrols, "log10_counts_feature_controls", NA),
               ifelse(withcontrols, "pct_counts_feature_controls", NA))

print(plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)])) # the cell type and culture condition is the main source of variation

# Pairs plot, showing correlation between some explanatory variables

print(plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)], method = "pairs")) # negativ association between dropout and total features, clear association between total counts and total features 


# PLot 
print(plotPhenoData(sceset, aes_string(x = "log10_total_counts", y = "total_features",
                                       colour = paste("source_name_ch1", collapse = "."))) + 
        guides(colour = guide_legend( byrow = TRUE)) + theme(legend.position = "bottom")





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
#sceset <-plotPCA(sceset,
#                 size_by = "total_features", 
#                 pca_data_input = "pdata",
#                 detect_outliers = TRUE,
#                 return_SCESet = TRUE)

# Dealing with confounders
