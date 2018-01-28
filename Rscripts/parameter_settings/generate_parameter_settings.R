## --------------------------------- NOTE!! --------------------------------- ##
## This script is only a convenience script to generate json files with
## parameter settings for each data set/method combination. Note that each time
## a parameter file is updated, this will trigger the regeneration of the
## corresponding clustering results. Thus, run only the lines that are modified.
## -------------------------------------------------------------------------- ##

## ZINB-WaVE parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(n_genes = 1000, k = 3)), 
      file = "parameter_settings/sce_full_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 4)), 
      file = "parameter_settings/sce_full_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 9)), 
      file = "parameter_settings/sce_full_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200, k = 4)), 
      file = "parameter_settings/sce_full_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 4)), 
      file = "parameter_settings/sce_full_SimKumar_zinbwave.json")

write(toJSON(list(n_genes = 1000, k = 3)), 
      file = "parameter_settings/sce_filtered_Kumar_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 3)), 
      file = "parameter_settings/sce_filtered_Trapnell_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 9)), 
      file = "parameter_settings/sce_filtered_Koh_zinbwave.json")
write(toJSON(list(n_genes = 200, k = 4)), 
      file = "parameter_settings/sce_filtered_Zhengmix_zinbwave.json")
write(toJSON(list(n_genes = 1000, k = 4)), 
      file = "parameter_settings/sce_filtered_SimKumar_zinbwave.json")

## TSCAN parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(minexpr_percent = 0.5, range_clusters = 2:10, n_clusters = 3)),
      file = "parameter_settings/sce_full_Kumar_tscan.json")
write(toJSON(list(minexpr_percent = 0.5, range_clusters = 2:10, n_clusters = 3)),
      file = "parameter_settings/sce_full_Trapnell_tscan.json")
write(toJSON(list(minexpr_percent = 0.5, range_clusters = 2:15, n_clusters = 9)),
      file = "parameter_settings/sce_full_Koh_tscan.json")
write(toJSON(list(minexpr_percent = 0.5, range_clusters = 2:10, n_clusters = 4)),
      file = "parameter_settings/sce_full_Zhengmix_tscan.json")
write(toJSON(list(minexpr_percent = 0.5, range_clusters = 2:10, n_clusters = 4)),
      file = "parameter_settings/sce_full_SimKumar_tscan.json")

write(toJSON(list(minexpr_percent = 0, range_clusters = 2:10, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Kumar_tscan.json")
write(toJSON(list(minexpr_percent = 0, range_clusters = 2:10, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Trapnell_tscan.json")
write(toJSON(list(minexpr_percent = 0, range_clusters = 2:15, n_clusters = 9)),
      file = "parameter_settings/sce_filtered_Koh_tscan.json")
write(toJSON(list(minexpr_percent = 0, range_clusters = 2:10, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_Zhengmix_tscan.json")
write(toJSON(list(minexpr_percent = 0, range_clusters = 2:10, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_SimKumar_tscan.json")

## CIDR parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 3)),
      file = "parameter_settings/sce_full_Kumar_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 3)),
      file = "parameter_settings/sce_full_Trapnell_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 0)),
      file = "parameter_settings/sce_full_Koh_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 4)),
      file = "parameter_settings/sce_full_Zhengmix_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 4)),
      file = "parameter_settings/sce_full_SimKumar_cidr.json")

write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Kumar_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Trapnell_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 0)),
      file = "parameter_settings/sce_filtered_Koh_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_Zhengmix_cidr.json")
write(toJSON(list(n_pc = 4, range_clusters = NULL, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_SimKumar_cidr.json")

## SIMLR parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(normalize = FALSE, n_clusters = 3)),
      file = "parameter_settings/sce_full_Kumar_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 3)),
      file = "parameter_settings/sce_full_Trapnell_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 0)),
      file = "parameter_settings/sce_full_Koh_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 4)),
      file = "parameter_settings/sce_full_Zhengmix_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 4)),
      file = "parameter_settings/sce_full_SimKumar_simlr.json")

write(toJSON(list(normalize = FALSE, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Kumar_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 3)),
      file = "parameter_settings/sce_filtered_Trapnell_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 0)),
      file = "parameter_settings/sce_filtered_Koh_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_Zhengmix_simlr.json")
write(toJSON(list(normalize = FALSE, n_clusters = 4)),
      file = "parameter_settings/sce_filtered_SimKumar_simlr.json")

## RtsneKmeans parameters
## -------------------------------------------------------------------------- ##
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Kumar_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Trapnell_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:15, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Koh_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_Zhengmix_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_full_SimKumar_RtsneKmeans.json")

write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Kumar_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Trapnell_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:15, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Koh_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_Zhengmix_RtsneKmeans.json")
write(toJSON(list(range_clusters = 2:10, perplexity = 30, initial_dims = 50)),
      file = "parameter_settings/sce_filtered_SimKumar_RtsneKmeans.json")

## pcaReduce parameters
## -------------------------------------------------------------------------- ##
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
