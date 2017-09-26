###########################################################
# Load and process data set: Xue 2013 GSE44183-GPL11154
###########################################################

pdf("results/QC_data/QC_GSE44183-GPL11154.pdf")

# load packages
source("skript/helper_files/Helper_functions.R")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
# load multiassay dataset from IMLS conquer
maex <- readRDS("data/GSE44183-GPL11154.rds")
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


#keep only those features that are observed in at least 1 cell.
keep_features <- rowSums(counts(sceset) > 0) >= 1

sceset <- sceset[keep_features, ]

# Plot QC metrics
plot_QC(sceset)


#### rename cell labels
# leave it for now

#### save cell labels

labels <- phenoData(sceset)$phenoid
dir_labels <-  paste0("results/dataset_labels/labels_xue2013.txt")
write.table(labels , file=dir_labels, sep="\t")


###### save as SCEobject and the tpm Matrix as Text file
res <- sceset
save(res,file = "data/sceset_GSE44183-GPL11154.rda")


###### save session info
sink(file = "results/QC_data/session_info_GSE44183-GPL11154.txt")
sessionInfo()
sink()

dev.off()

tsne <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], exprs_values = "tpm", return_SCESet = TRUE)

pdf("results/QC_data/Xue2013pca.pdf")
print(scater::plotReducedDim(tsne, colour_by = "source_name_ch1") + 
        guides(fill = guide_legend(byrow = TRUE)) )
dev.off()



tsne <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], exprs_values = "tpm", return_SCESet = TRUE)
phenoid <- as.data.frame(pData(tsne)$phenoid)

tsne.data <- as.data.frame(tsne@reducedDimension)
tsne.data <- cbind(tsne.data, phenoid)
colnames(tsne.data) <- c("Dimension 1","Dimension 2", "phenoid")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000001","#000002","#000003","#000005")

pdf("results/QC_data/Xue2013tsne.pdf")
ggplot(data = tsne.data , mapping = aes(x=`Dimension 1`,y=`Dimension 2`))+
  geom_point(aes(colour=phenoid), size=2)+labs(colour=phenoid)+scale_colour_manual(values=cbbPalette)

dev.off()



## Appendix
### Automatic Cell filtering
#sceset <-plotPCA(sceset,
#                 size_by = "total_features", 
#                 pca_data_input = "pdata",
#                 detect_outliers = TRUE,
#                 return_SCESet = TRUE)

