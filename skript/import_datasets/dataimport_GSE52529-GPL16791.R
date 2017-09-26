###########################################################
# Load and process data set: Trapnell2014 GSE52529-GPL16791
###########################################################

pdf("results/QC_data/QC_GSE52529-GPL16791.pdf")

# load packages
source("skript/helper_functions/Helper_functions.R")

suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(mvoutlier))
suppressPackageStartupMessages(library(plyr))
library(cowplot)

# load multiassay dataset from IMLS conquer
maex <- readRDS("data/GSE52529-GPL16791.rds")
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

# rename cell levels
phn[, "phenoid"] <- mapvalues(phn[, "phenoid"], from = c(grep("T0",phn[, "phenoid"], value=TRUE),grep("T24",phn[, "phenoid"], value=TRUE),grep("T48",phn[, "phenoid"], value=TRUE)) , to=c(array("T0", dim=96),array("T24", dim=96),array("T48", dim=96)) )

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


#keep only those features that are observed in at least 1 cell.
keep_features <- rowSums(counts(sceset) > 0) >= 1



sceset <- sceset[keep_features, ]
dim(sceset) # 41112 x 288

# Plot QC metrics
# plot frequency of total counts from top 50 expressed genes

plot_QC(sceset)



#### save cell labels

labels <- phenoData(sceset)$phenoid
dir_labels <-  paste0("results/dataset_labels/labels_trapnell2014.txt")
write.table(labels , file=dir_labels, sep="\t")


###### save as SCEobject 
res <- sceset
save(res,file = "data/sceset_GSE52529-GPL16791.rda")
###### save session info
sink(file = "results/QC_data/session_info_GSE52529-GPL16791.txt")
sessionInfo()
sink()
dev.off()



# Appendix


tsne <- scater::plotTSNE(sceset[!duplicated(tpm(sceset)), ], exprs_values = "tpm", return_SCESet = TRUE)
phenoid <- as.data.frame(pData(tsne)$phenoid)

tsne.data <- as.data.frame(tsne@reducedDimension)
tsne.data <- cbind(tsne.data, phenoid)
colnames(tsne.data) <- c("Dimension 1","Dimension 2", "phenoid")

pdf("results/QC_data/Trapnell2014tsne.pdf")
ggplot(data = tsne.data , mapping = aes(x=`Dimension 1`,y=`Dimension 2`))+
  geom_point(aes(colour=phenoid), size=2)+labs(colour=phenoid)+scale_colour_manual(values=c("#000000", "#E69F00", "#0072B2"))

dev.off()

dim(sceset) # 41112    x  288 
# Dealing with confounders


# dim of reduced data set
dim(sceset[fData(sceset)$outlier, pData(sceset)$outlier==FALSE])

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
