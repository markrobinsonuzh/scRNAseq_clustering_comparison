###########################################################
# Load and process data set: Kumar 2014 GSE60749-GPL13112
###########################################################
pdf("results/QC_data/QC_GSE60749-GPL13112.pdf")


# load packages
source("skript/helper_functions/Helper_functions.R")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
# load multiassay dataset from IMLS conquer
maex <- readRDS("data/GSE60749-GPL13112.rds")
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
phn[, "phenoid"] <- as.character(interaction(as.data.frame(phn[, c("source_name_ch1","characteristics_ch1.1")])))
phn[, "phenoid"] <-  revalue(phn[, "phenoid"] , c("Dgcr8 knockout mouse embryonic stem cells.culture conditions: serum+LIF"="Dgcr8 knockout mouse serum+LIF", "v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF"="v6.5 mouse 2i+LIF","v6.5 mouse embryonic stem cells.culture conditions: serum+LIF"="v6.5 mouse serum+LIF"))


table(phn[, "phenoid"])


# Create an SCEobject 
sceset <- newSCESet(countData = cts, 
                   phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
set_exprs(sceset, "tpm") <- tpms


#################
# QC 
##################
#Plot overview for each cell

  plot(sceset, block1 = paste("source_name_ch1", collapse = "."))


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
head(colnames(pData(sceset)))
pData(sceset)

#keep only those features that are observed in at least 1 cell.
keep_features <- rowSums(counts(sceset) > 0) >= 1



sceset <- sceset[keep_features, ]

# Plot QC metrics

plot_QC(sceset)





###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "data/sceset_GSE60749-GPL13112.rda")

###### save session info
sink(file = "results/QC_data/session_info_GSE60749-GPL13112.txt")
sessionInfo()
sink()

dev.off()






## Appendix
tsne <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], exprs_values = "tpm", return_SCESet = TRUE)


pdf("results/QC_data/Kumar2014tsne.pdf")
print(scater::plotReducedDim(tsne, colour_by = "phenoid") + 
        guides(fill = guide_legend(byrow = TRUE)) )
dev.off()





tsne <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], exprs_values = "tpm", return_SCESet = TRUE)
phenoid <- as.data.frame(pData(tsne)$phenoid)

tsne.data <- as.data.frame(tsne@reducedDimension)
tsne.data <- cbind(tsne.data, phenoid)
colnames(tsne.data) <- c("Dimension 1","Dimension 2", "phenoid")

pdf("results/QC_data/Kumar2014tsne.pdf")
ggplot(data = tsne.data , mapping = aes(x=`Dimension 1`,y=`Dimension 2`))+
  geom_point(aes(colour=phenoid), size=2)+labs(colour=phenoid)+scale_colour_manual(values=c("#000000", "#E69F00", "#0072B2"))

dev.off()


$### Automatic Cell filtering
sceset <-plotPCA(sceset,
                 size_by = "total_features", 
                 pca_data_input = "pdata",
                 detect_outliers = TRUE,
                 return_SCESet = TRUE)

# Dealing with confounders
