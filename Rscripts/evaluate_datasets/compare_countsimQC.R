#######################################################################
### Compare simulations with full Kumar dataset using countsimQC 
##########################################################################
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(datasets)
print(filtering)
print(outrds)

suppressPackageStartupMessages({
  require(SingleCellExperiment)
  require(splatter)
  require(countsimQC)
  require(DESeq2)
})

# load the datasets
data <- list(
  Kumarfull     = readRDS(paste0("/Volumes/Shared/data/seq/scRNAseq_clustering_comparison/data/sce_full/sce_full_Kumar.rds")),
  SimKumar4easy = readRDS(paste0("/Volumes/Shared/data/seq/scRNAseq_clustering_comparison/data/sce_full/sce_full_SimKumar4easy.rds")),
  SimKumar4hard = readRDS(paste0("/Volumes/Shared/data/seq/scRNAseq_clustering_comparison/data/sce_full/sce_full_SimKumar4hard.rds")),
  SimKumar8hard = readRDS(paste0("/Volumes/Shared/data/seq/scRNAseq_clustering_comparison/data/sce_full/sce_full_SimKumar8hard.rds"))
)

# create List

ddsList<- lapply( data, function(x){DESeqDataSetFromMatrix(countData = round( counts(x), 0), 
                                                 colData = colData(x),
                                                 design = ~ as.factor(phenoid))}
)


# Report
countsimQCReport(ddsList = ddsList, outputFile = "countsim_report.html", 
                 outputDir = "./comparison_countsim", outputFormat = "html_document", 
                 showCode = FALSE, forceOverwrite = TRUE,
                 savePlots = TRUE, 
                 description = "This is a comparison of the Kumar full and the 
                 simulated counts thereof. Simulation was done using the Splatter package, 
                 three groups were simulated.", 
                 maxNForCorr = 25, maxNForDisp = 50, 
                 calculateStatistics = TRUE, subsampleSize = 25,
                 kfrac = 0.01, kmin = 5, 
                 permutationPvalues = FALSE, nPermutations = NULL)


