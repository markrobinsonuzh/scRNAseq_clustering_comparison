###########################################################
# Load and process data set: Trapnell2014 GSE52529-GPL16791
###########################################################
# load packages
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(mvoutlier))
# load multiassay dataset from IMLS conquer
maex <- readRDS("~/Desktop/masterthesis/data/GSE52529-GPL16791.rds")
# summary of data set
print(maex)
dim(assay(maex)) # 65218 x 288


# extract the gene-level TPM and the counts obtaines by scaling the TPMs, else counts

if (all(c("count_lstpm", "TPM") %in% names(assays(experiments(maex)[["gene"]])))) {
  cts <- assays(experiments(maex)[["gene"]])[["count_lstpm"]]
  tpms <- assays(experiments(maex)[["gene"]])[["TPM"]]
} else {
  cts <- assays(experiments(maex)[["gene"]])[["count"]]
  tpms <- NULL
}


# extract the phenotype data

phn <- colData(maex)
phn@listData$phenoid <- array(data=NA,dim = phn@nrows)

# extract the cell annotation to phenoid col
phn[, "phenoid"] <- as.character(interaction(as.data.frame(phn[, "source_name_ch1"])))

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
dim(sceset) # 41112 x 288

# Plot QC metrics
# plot frequency of total counts from top 50 expressed genes

plotQC(sceset, type = "highest-expression") # Top 50 genes account for 18.4 % of the counts, profile relatively flat

plotQC(sceset, type = "find-pcs", variable = "total_features") # PC1 shows high correlation with total number observerd features per cell
plotQC(sceset, type = "find-pcs", variable = "total_features", plot_type = "pairs-pcs") # PC1 and PC75, no correlation between other PCs


# plot the variance explained for some explanatory vars.
expl_vars <- c("source_name_ch1", "log10_total_counts", "log10_total_features", "pct_dropout",
               "pct_counts_top_200_features", ifelse(withcontrols, "log10_counts_feature_controls", NA),
               ifelse(withcontrols, "pct_counts_feature_controls", NA))

plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)]) # the cell type and culture condition is the main source of variation

# Pairs plot, showing correlation between some explanatory variables

plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)], method = "pairs") # negativ association between dropout and total features, clear association between total counts and total features 


# PLot 
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
# 44 cells removed
knitr::kable(
  as.data.frame(table(sceset$outlier)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by automatic filter (FALSE)'
)


dim(sceset) # 
# Dealing with confounders


# dim of reduced data set
dim(sceset[fData(sceset)$outlier, pData(sceset)$outlier==FALSE])

###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "~/Desktop/masterthesis/data/sceset_GSE52529-GPL16791.rds")
###### save session info
sink(file = "~/Desktop/masterthesis/results/session_info_GSE52529-GPL16791.txt")
sessionInfo()
sink()
