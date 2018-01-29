## --------------------------------- NOTE!! --------------------------------- ##
## This script is only a convenience script to generate json files with
## parameter settings for each data set/method combination. Note that each time
## a parameter file is updated, this will trigger the regeneration of the
## corresponding clustering results. Thus, run only the lines that are modified.
## -------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(rjson)
})

## General parameters (range of cluster numbers etc)
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
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_full_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_full_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_full_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200)), file = "parameter_settings/sce_full_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_full_SimKumar_zinbwave.json")

write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_filtered_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_filtered_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_filtered_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200)), file = "parameter_settings/sce_filtered_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000)), file = "parameter_settings/sce_filtered_SimKumar_zinbwave.json")

## TSCAN parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(minexpr_percent = 0.5)), file = "parameter_settings/sce_full_Kumar_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)), file = "parameter_settings/sce_full_Trapnell_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)), file = "parameter_settings/sce_full_Koh_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)), file = "parameter_settings/sce_full_Zhengmix_TSCAN.json")
write(toJSON(list(minexpr_percent = 0.5)), file = "parameter_settings/sce_full_SimKumar_TSCAN.json")

write(toJSON(list(minexpr_percent = 0)), file = "parameter_settings/sce_filtered_Kumar_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)), file = "parameter_settings/sce_filtered_Trapnell_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)), file = "parameter_settings/sce_filtered_Koh_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)), file = "parameter_settings/sce_filtered_Zhengmix_TSCAN.json")
write(toJSON(list(minexpr_percent = 0)), file = "parameter_settings/sce_filtered_SimKumar_TSCAN.json")

## CIDR parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_full_Kumar_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_full_Trapnell_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_full_Koh_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_full_Zhengmix_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_full_SimKumar_cidr.json")

write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_filtered_Kumar_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_filtered_Trapnell_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_filtered_Koh_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_filtered_Zhengmix_cidr.json")
write(toJSON(list(nPC = 4)), file = "parameter_settings/sce_filtered_SimKumar_cidr.json")

## SIMLR parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Kumar_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Trapnell_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Koh_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Zhengmix_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_SimKumar_SIMLR.json")

write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_filtered_Kumar_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Trapnell_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Koh_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Zhengmix_SIMLR.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_SimKumar_SIMLR.json")

## RtsneKmeans parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_SimKumar_RtsneKmeans.json")

write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Kumar_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Trapnell_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Koh_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Zhengmix_RtsneKmeans.json")
write(toJSON(list(perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_SimKumar_RtsneKmeans.json")

## pcaReduce parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_full_Kumar_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_full_Trapnell_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_full_Koh_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_full_Zhengmix_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_full_SimKumar_pcaReduce.json")

write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_filtered_Kumar_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_filtered_Trapnell_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_filtered_Koh_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_filtered_Zhengmix_pcaReduce.json")
write(toJSON(list(nbt = 1, q = 30)), file = "parameter_settings/sce_filtered_SimKumar_pcaReduce.json")

## DBSCAN parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_full_Kumar_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_full_Trapnell_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_full_Koh_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_full_Zhengmix_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_full_SimKumar_DBSCAN.json")

write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_filtered_Kumar_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_filtered_Trapnell_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_filtered_Koh_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_filtered_Zhengmix_DBSCAN.json")
write(toJSON(list(eps = , Pts = )), file = "parameter_settings/sce_filtered_SimKumar_DBSCAN.json")

## Seurat parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_full_Kumar_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_full_Trapnell_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_full_Koh_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_full_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_full_SimKumar_Seurat.json")

write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )),
      file = "parameter_settings/sce_filtered_Kumar_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_filtered_Trapnell_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_filtered_Koh_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )), 
      file = "parameter_settings/sce_filtered_Zhengmix_Seurat.json")
write(toJSON(list(min.cells = , min.genes = , dims.use = , resolution = )),
      file = "parameter_settings/sce_filtered_SimKumar_Seurat.json")

## SIMLRlargescale parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Kumar_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Trapnell_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Koh_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_Zhengmix_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_full_SimKumar_SIMLRlargescale.json")

write(toJSON(list(k = 10, normalize = FALSE)), 
      file = "parameter_settings/sce_filtered_Kumar_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Trapnell_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Koh_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_Zhengmix_SIMLRlargescale.json")
write(toJSON(list(k = 10, normalize = FALSE)),
      file = "parameter_settings/sce_filtered_SimKumar_SIMLRlargescale.json")

## Linnorm parameters
## -------------------------------------------------------------------------- ##
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
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_full_Kumar_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_full_Trapnell_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_full_Koh_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_full_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_full_SimKumar_RaceID.json")

write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_filtered_Kumar_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_filtered_Trapnell_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_filtered_Koh_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_filtered_Zhengmix_RaceID.json")
write(toJSON(list(mintotal = , minexprs = , minnumber = , maxexpr = , cln = , do.gap = )), 
      file = "parameter_settings/sce_filtered_SimKumar_RaceID.json")

## SC3 parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_full_Kumar_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_full_Trapnell_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_full_Koh_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_full_Zhengmix_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_full_SimKumar_SC3.json")

write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_filtered_Kumar_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_filtered_Trapnell_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_filtered_Koh_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_filtered_Zhengmix_SC3.json")
write(toJSON(list(ks = , estimate_k = , pct_dropout_max = )), 
      file = "parameter_settings/sce_filtered_SimKumar_SC3.json")
