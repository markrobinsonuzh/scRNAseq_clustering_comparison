## Compare simulations with full Kumar dataset using countsimQC 

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(datadir)
print(outdir)
print(outfile)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
  library(countsimQC)
  library(DESeq2)
})

## load the datasets
data <- list(
  Kumarfull = readRDS(paste0(datadir, "/sce_full/sce_full_Kumar.rds")),
  SimKumar4easy = readRDS(paste0(datadir, "/sce_full/sce_full_SimKumar4easy.rds")),
  SimKumar4hard = readRDS(paste0(datadir, "/sce_full/sce_full_SimKumar4hard.rds")),
  SimKumar8hard = readRDS(paste0(datadir, "/sce_full/sce_full_SimKumar8hard.rds"))
)

## create List of DESeqDataSets
ddsList <- lapply(data, function(x){ 
  DESeqDataSetFromMatrix(countData = round(counts(x), 0), 
                         colData = colData(x),
                         design = ~ as.factor(phenoid))
})

## Generate report
countsimQCReport(ddsList = ddsList, outputFile = outfile, 
                 outputDir = outdir, outputFormat = "html_document", 
                 showCode = FALSE, forceOverwrite = TRUE,
                 savePlots = FALSE, 
                 description = "This is a comparison of the unfiltered Kumar data set and the 
                 simulated data sets derived from one of the cell subpopulations of the Kumar data set. 
                 Simulation was done using the Splatter package, generating three 
                 cell populations for each of the simulated data sets.", 
                 maxNForCorr = 25, maxNForDisp = 50, 
                 calculateStatistics = TRUE, subsampleSize = 25,
                 kfrac = 0.01, kmin = 5, 
                 permutationPvalues = FALSE, nPermutations = NULL)

date()
sessionInfo()

