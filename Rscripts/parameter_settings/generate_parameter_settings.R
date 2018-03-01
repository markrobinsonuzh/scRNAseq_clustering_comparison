## --------------------------------- NOTE!! --------------------------------- ##
## This script is only a convenience script to generate json files with
## parameter settings for each data set/method combination. Note that each time
## a parameter file is updated, this will trigger the regeneration of the
## corresponding clustering results. Thus, run only the lines that are modified.
## -------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(rjson)
})

filterings <- c("full", "filteredExpr50", "filteredExpr10", "filteredHVG50", 
                "filteredHVG10")
datasets <- c("Kumar", "Trapnell", "Koh", "Zhengmix4eq", "Zhengmix4uneq", 
              "Zhengmix8eq", "SimKumar4easy", "SimKumar4hard", "SimKumar8hard")

## General dataset-specific (but method-agnostic) parameters (range of cluster numbers etc)
## -------------------------------------------------------------------------- ##
for (f in filterings) {
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_Kumar.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_Trapnell.json"))
  write(toJSON(list(range_clusters = 2:15)), 
        file = paste0("parameter_settings/sce_", f, "_Koh.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_Zhengmix4eq.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_Zhengmix4uneq.json"))
  write(toJSON(list(range_clusters = 2:15)), 
        file = paste0("parameter_settings/sce_", f, "_Zhengmix8eq.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_SimKumar4easy.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_SimKumar4hard.json"))
  write(toJSON(list(range_clusters = 2:15)), 
        file = paste0("parameter_settings/sce_", f, "_SimKumar8hard.json"))
}

## CIDR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/CIDR.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_CIDR.json"))
  }
}

## TSCAN parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/TSCAN.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_TSCAN.json"))
  }
}

## RtsneKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(initial_dims = 50)), file = "parameter_settings/RtsneKmeans.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(perplexity = 30, dims = 3)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_RtsneKmeans.json"))
  }
}

## pcaReduce parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nbt = 100)), file = "parameter_settings/pcaReduce.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(q = 30)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_pcaReduce.json"))
  }
}

## Seurat parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(range_resolutions = seq(0.3, 1.5, by = 0.1))), 
      file = "parameter_settings/Seurat.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(min.cells = 0, min.genes = 0, dims.use = 1:30)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
  }
}

## SIMLR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLR.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SIMLR.json"))
  }
}

## SIMLRlargescale parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(k = 10)), file = "parameter_settings/SIMLRlargescale.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SIMLRlargescale.json"))
  }
}

## Linnorm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/Linnorm.json")

## Dataset-specific
for (f in filterings) {
  for (d in setdiff(datasets, c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq"))) {
    write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_Linnorm.json"))
  }
}

for (f in filterings) {
  for (d in c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq")) {
    write(toJSON(list(minNonZeroPortion = 0.1, num_PC = 3)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_Linnorm.json"))
  }
}

## RaceID parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/RaceID.json")

## Dataset-specific
for (f in filterings) {
  for (d in setdiff(datasets, c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq"))) {
    write(toJSON(list(mintotal = 3000, minexpr = 5, minnumber = 1, maxexpr = Inf)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_RaceID.json"))
  }
}

for (f in filterings) {
  for (d in c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq")) {
    write(toJSON(list(mintotal = 200, minexpr = 1, minnumber = 1, maxexpr = Inf)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_RaceID.json"))
  }
}

## SC3 parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SC3.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SC3.json"))
  }
}

## SC3svm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SC3svm.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(pct_dropout_min = 10, pct_dropout_max = 90)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SC3svm.json"))
  }
}

## FlowSOM parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 50)), file = "parameter_settings/FlowSOM.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(xdim = 10, ydim = 10)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_FlowSOM.json"))
  }
}

## PCAKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 10)), file = "parameter_settings/PCAKmeans.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_PCAKmeans.json"))
  }
}

## PCAHC parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 10)), file = "parameter_settings/PCAHC.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_PCAHC.json"))
  }
}
