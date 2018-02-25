## --------------------------------- NOTE!! --------------------------------- ##
## This script is only a convenience script to generate json files with
## parameter settings for each data set/method combination. Note that each time
## a parameter file is updated, this will trigger the regeneration of the
## corresponding clustering results. Thus, run only the lines that are modified.
## -------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(rjson)
})

## General dataset-specific (but method-agnostic) parameters (range of cluster numbers etc)
## -------------------------------------------------------------------------- ##
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_Kumar.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_Trapnell.json")
write(toJSON(list(range_clusters = 2:15)), file = "parameter_settings/sce_full_Koh.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_Zhengmix4eq.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_Zhengmix4uneq.json")
write(toJSON(list(range_clusters = 2:15)), file = "parameter_settings/sce_full_Zhengmix8eq.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_SimKumarEasy.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_SimKumarHard.json")

write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_Kumar.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_Trapnell.json")
write(toJSON(list(range_clusters = 2:15)), file = "parameter_settings/sce_filteredExpr_Koh.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_Zhengmix4eq.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_Zhengmix4uneq.json")
write(toJSON(list(range_clusters = 2:15)), file = "parameter_settings/sce_filteredExpr_Zhengmix8eq.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_SimKumarEasy.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filteredExpr_SimKumarHard.json")

## CIDR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/CIDR.json")

## Dataset-specific
write(toJSON(list()), file = "parameter_settings/sce_full_Kumar_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Trapnell_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Koh_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Zhengmix_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_SimKumar_CIDR.json")

write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Kumar_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Trapnell_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Koh_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Zhengmix_CIDR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_SimKumar_CIDR.json")

## TSCAN parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/TSCAN.json")

## Dataset-specific
write(toJSON(list()), file = "parameter_settings/sce_full_Kumar_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Trapnell_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Koh_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Zhengmix_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_full_SimKumar_TSCAN.json")

write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Kumar_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Trapnell_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Koh_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Zhengmix_TSCAN.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_SimKumar_TSCAN.json")

## RtsneKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(initial_dims = 50)), file = "parameter_settings/RtsneKmeans.json")

## Dataset-specific
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_full_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_full_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_full_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_full_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_full_SimKumar_RtsneKmeans.json")

write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_filteredExpr_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_filteredExpr_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_filteredExpr_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_filteredExpr_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, dims = 3)),
      file = "parameter_settings/sce_filteredExpr_SimKumar_RtsneKmeans.json")

## pcaReduce parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nbt = 100)), file = "parameter_settings/pcaReduce.json")

## Dataset-specific
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_full_Kumar_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_full_Trapnell_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_full_Koh_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_full_Zhengmix_pcaReduce.json")
write(toJSON(list(q = 30)),
      file = "parameter_settings/sce_full_SimKumar_pcaReduce.json")

write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_filteredExpr_Koh_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_pcaReduce.json")
write(toJSON(list(q = 30)), 
      file = "parameter_settings/sce_filteredExpr_SimKumar_pcaReduce.json")

## Seurat parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(range_resolutions = seq(0.3, 1.5, by = 0.1))), 
      file = "parameter_settings/Seurat.json")

## Dataset-specific
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_full_Kumar_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_full_Trapnell_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_full_Koh_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_full_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_full_SimKumar_Seurat.json")

write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)),
      file = "parameter_settings/sce_filteredExpr_Kumar_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_filteredExpr_Koh_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)),
      file = "parameter_settings/sce_filteredExpr_SimKumar_Seurat.json")

## SIMLR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLR.json")

## Dataset-specific
write(toJSON(list()), file = "parameter_settings/sce_full_Kumar_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Trapnell_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Koh_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Zhengmix_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_full_SimKumar_SIMLR.json")

write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Kumar_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Trapnell_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Koh_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Zhengmix_SIMLR.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_SimKumar_SIMLR.json")

## SIMLRlargescale parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLRlargescale.json")

## Dataset-specific
write(toJSON(list()), file = "parameter_settings/sce_full_Kumar_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Trapnell_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Koh_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Zhengmix_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_full_SimKumar_SIMLRlargescale.json")

write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Kumar_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Trapnell_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Koh_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Zhengmix_SIMLRlargescale.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_SimKumar_SIMLRlargescale.json")

## Linnorm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/Linnorm.json")

## Dataset-specific
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_full_Kumar_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_full_Trapnell_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_full_Koh_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.1, num_PC = 3)), 
      file = "parameter_settings/sce_full_Zhengmix_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_full_SimKumar_Linnorm.json")

write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_filteredExpr_Koh_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.1, num_PC = 3)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
      file = "parameter_settings/sce_filteredExpr_SimKumar_Linnorm.json")

## RaceID parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/RaceID.json")

## Dataset-specific
write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_full_Kumar_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr1 = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_full_Trapnell_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr1 = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_full_Koh_RaceID.json")
write(toJSON(list(mintotal = 200, minexpr = 1, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_full_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_full_SimKumar_RaceID.json")

write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_filteredExpr_Koh_RaceID.json")
write(toJSON(list(mintotal = 200, minexpr = 1, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
      file = "parameter_settings/sce_filteredExpr_SimKumar_RaceID.json")

## SC3 parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SC3.json")

## Dataset-specific
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Kumar_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Trapnell_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Koh_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Zhengmix_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_SimKumar_SC3.json")

write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Koh_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_SC3.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_SimKumar_SC3.json")

## SC3svm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SC3svm.json")

## Dataset-specific
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Kumar_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Trapnell_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Koh_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Zhengmix_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_SimKumar_SC3svm.json")

write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Trapnell_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Koh_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_Zhengmix_SC3svm.json")
write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filteredExpr_SimKumar_SC3svm.json")

## FlowSOM parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 50)), file = "parameter_settings/FlowSOM.json")

## Dataset-specific
write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_full_Kumar_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_full_Trapnell_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_full_Koh_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_full_Zhengmix_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_full_SimKumar_FlowSOM.json")

write(toJSON(list(xdim = 10, ydim = 10)), 
      file = "parameter_settings/sce_filteredExpr_Kumar_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)),
      file = "parameter_settings/sce_filteredExpr_Trapnell_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)),
      file = "parameter_settings/sce_filteredExpr_Koh_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)),
      file = "parameter_settings/sce_filteredExpr_Zhengmix_FlowSOM.json")
write(toJSON(list(xdim = 10, ydim = 10)),
      file = "parameter_settings/sce_filteredExpr_SimKumar_FlowSOM.json")

## PCAKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 10)), file = "parameter_settings/PCAKmeans.json")

## Dataset-specific
write(toJSON(list()), file = "parameter_settings/sce_full_Kumar_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Trapnell_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Koh_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_full_Zhengmix_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_full_SimKumar_PCAKmeans.json")

write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Kumar_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Trapnell_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Koh_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_Zhengmix_PCAKmeans.json")
write(toJSON(list()), file = "parameter_settings/sce_filteredExpr_SimKumar_PCAKmeans.json")
