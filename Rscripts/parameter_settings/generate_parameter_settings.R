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
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_Zheng.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_full_SimKumar.json")

write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filtered_Kumar.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filtered_Trapnell.json")
write(toJSON(list(range_clusters = 2:15)), file = "parameter_settings/sce_filtered_Koh.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filtered_Zheng.json")
write(toJSON(list(range_clusters = 2:10)), file = "parameter_settings/sce_filtered_SimKumar.json")

## ZINB-WaVE parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/zinbwave.json")

## Dataset-specific
write(toJSON(list(n_genes = 1000)),
      file = "parameter_settings/sce_full_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000)), 
      file = "parameter_settings/sce_full_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000)), 
      file = "parameter_settings/sce_full_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200)), 
      file = "parameter_settings/sce_full_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000)), 
      file = "parameter_settings/sce_full_SimKumar_zinbwave.json")

write(toJSON(list(n_genes = 1000)),
      file = "parameter_settings/sce_filtered_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000)),
      file = "parameter_settings/sce_filtered_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000)),
      file = "parameter_settings/sce_filtered_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200)),
      file = "parameter_settings/sce_filtered_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000)),
      file = "parameter_settings/sce_filtered_SimKumar_zinbwave.json")

## TSCAN parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/TSCAN.json")

## Dataset-specific
write(toJSON(list(minexpr_percent = 0.5)),
      file = "parameter_settings/sce_full_Kumar_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)),
      file = "parameter_settings/sce_full_Trapnell_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)),
      file = "parameter_settings/sce_full_Koh_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)),
      file = "parameter_settings/sce_full_Zhengmix_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)),
      file = "parameter_settings/sce_full_SimKumar_TSCAN.json")

write(toJSON(list(minexpr_percent = 0)),
      file = "parameter_settings/sce_filtered_Kumar_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)),
      file = "parameter_settings/sce_filtered_Trapnell_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)),
      file = "parameter_settings/sce_filtered_Koh_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)),
      file = "parameter_settings/sce_filtered_Zhengmix_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)),
      file = "parameter_settings/sce_filtered_SimKumar_TSCAN.json")

## CIDR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/CIDR.json")

## Dataset-specific
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_full_Kumar_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_full_Trapnell_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_full_Koh_CIDR.json")
write(toJSON(list(nPC = 4)),
      file = "parameter_settings/sce_full_Zhengmix_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_full_SimKumar_CIDR.json")

write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_filtered_Kumar_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_filtered_Trapnell_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_filtered_Koh_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_filtered_Zhengmix_CIDR.json")
write(toJSON(list(nPC = 4)), 
      file = "parameter_settings/sce_filtered_SimKumar_CIDR.json")

## RtsneKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(initial_dims = 50)), file = "parameter_settings/RtsneKmeans.json")

## Dataset-specific
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_full_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_full_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_full_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_full_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_full_SimKumar_RtsneKmeans.json")

write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_filtered_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_filtered_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_filtered_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_filtered_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30)),
      file = "parameter_settings/sce_filtered_SimKumar_RtsneKmeans.json")

## pcaReduce parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/pcaReduce.json")

## Dataset-specific
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_full_Kumar_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_full_Trapnell_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_full_Koh_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_full_Zhengmix_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)),
      file = "parameter_settings/sce_full_SimKumar_pcaReduce.json")

write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_filtered_Kumar_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_filtered_Trapnell_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_filtered_Koh_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_filtered_Zhengmix_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), 
      file = "parameter_settings/sce_filtered_SimKumar_pcaReduce.json")

## Seurat parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(range_resolutions = seq(0.3, 1.5, by = 0.1))), 
      file = "parameter_settings/Seurat.json")

## Dataset-specific
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_full_Kumar_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_full_Trapnell_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.7)), 
      file = "parameter_settings/sce_full_Koh_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_full_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_full_SimKumar_Seurat.json")

write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)),
      file = "parameter_settings/sce_filtered_Kumar_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_filtered_Trapnell_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.7)), 
      file = "parameter_settings/sce_filtered_Koh_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)), 
      file = "parameter_settings/sce_filtered_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:10, resolution = 0.6)),
      file = "parameter_settings/sce_filtered_SimKumar_Seurat.json")

## SIMLR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLR.json")

## Dataset-specific
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Kumar_SIMLR.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Trapnell_SIMLR.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Koh_SIMLR.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Zhengmix_SIMLR.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_SimKumar_SIMLR.json")

write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_filtered_Kumar_SIMLR.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Trapnell_SIMLR.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Koh_SIMLR.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Zhengmix_SIMLR.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_SimKumar_SIMLR.json")

## SIMLRlargescale parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLRlargescale.json")

## Dataset-specific
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Kumar_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Trapnell_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Koh_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_Zhengmix_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_full_SimKumar_SIMLRlargescale.json")

write(toJSON(list(normalize = FALSE)), 
      file = "parameter_settings/sce_filtered_Kumar_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Trapnell_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Koh_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Zhengmix_SIMLRlargescale.json")
write(toJSON(list(normalize = FALSE)),
      file = "parameter_settings/sce_filtered_SimKumar_SIMLRlargescale.json")

## Linnorm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/Linnorm.json")

## Dataset-specific
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_full_Kumar_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_full_Trapnell_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_full_Koh_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.1, BE_strength = 0.5)), 
      file = "parameter_settings/sce_full_Zhengmix_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_full_SimKumar_Linnorm.json")

write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_filtered_Kumar_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_filtered_Trapnell_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_filtered_Koh_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.1, BE_strength = 0.5)), 
      file = "parameter_settings/sce_filtered_Zhengmix_Linnorm.json")
write(toJSON(list(minNonZeroPortion = 0.75, BE_strength = 0.5)), 
      file = "parameter_settings/sce_filtered_SimKumar_Linnorm.json")

## RaceID parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/RaceID.json")

## Dataset-specific
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_full_Kumar_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_full_Trapnell_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_full_Koh_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = 500, cln = 0)), 
      file = "parameter_settings/sce_full_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_full_SimKumar_RaceID.json")

write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_filtered_Kumar_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_filtered_Trapnell_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_filtered_Koh_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = 500, cln = 0)), 
      file = "parameter_settings/sce_filtered_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = 3000, minexprs = 5, minnumber = 1, maxexpr = Inf, cln = 0)), 
      file = "parameter_settings/sce_filtered_SimKumar_RaceID.json")

## SC3 parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SC3.json")

## Dataset-specific
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Kumar_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Trapnell_SC3.json")
write(toJSON(list(ks = 2:15, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Koh_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_Zhengmix_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_full_SimKumar_SC3.json")

write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filtered_Kumar_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filtered_Trapnell_SC3.json")
write(toJSON(list(ks = 2:15, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filtered_Koh_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filtered_Zhengmix_SC3.json")
write(toJSON(list(ks = 2:10, pct_dropout_max = 90)), 
      file = "parameter_settings/sce_filtered_SimKumar_SC3.json")

## FlowSOM parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 50)), file = "parameter_settings/FlowSOM.json")

## Dataset-specific
write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_full_Kumar_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_full_Trapnell_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_full_Koh_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_full_Zhengmix_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_full_SimKumar_FlowSOM.json")

write(toJSON(list(xdim = 15, ydim = 15)), 
      file = "parameter_settings/sce_filtered_Kumar_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)),
      file = "parameter_settings/sce_filtered_Trapnell_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)),
      file = "parameter_settings/sce_filtered_Koh_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)),
      file = "parameter_settings/sce_filtered_Zhengmix_FlowSOM.json")
write(toJSON(list(xdim = 15, ydim = 15)),
      file = "parameter_settings/sce_filtered_SimKumar_FlowSOM.json")

## PCAKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 50)), file = "parameter_settings/PCAKmeans.json")

## Dataset-specific
write(toJSON(list()), 
      file = "parameter_settings/sce_full_Kumar_PCAKmeans.json")
write(toJSON(list()), 
      file = "parameter_settings/sce_full_Trapnell_PCAKmeans.json")
write(toJSON(list()), 
      file = "parameter_settings/sce_full_Koh_PCAKmeans.json")
write(toJSON(list()), 
      file = "parameter_settings/sce_full_Zhengmix_PCAKmeans.json")
write(toJSON(list()), 
      file = "parameter_settings/sce_full_SimKumar_PCAKmeans.json")

write(toJSON(list()), 
      file = "parameter_settings/sce_filtered_Kumar_PCAKmeans.json")
write(toJSON(list()),
      file = "parameter_settings/sce_filtered_Trapnell_PCAKmeans.json")
write(toJSON(list()),
      file = "parameter_settings/sce_filtered_Koh_PCAKmeans.json")
write(toJSON(list()),
      file = "parameter_settings/sce_filtered_Zhengmix_PCAKmeans.json")
write(toJSON(list()),
      file = "parameter_settings/sce_filtered_SimKumar_PCAKmeans.json")
