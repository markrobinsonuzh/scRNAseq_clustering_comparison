## --------------------------------- NOTE!! --------------------------------- ##
## This script is only a convenience script to generate json files with
## parameter settings for each data set/method combination. Note that each time
## a parameter file is updated, this will trigger the regeneration of the
## corresponding clustering results. Thus, run only the lines that are modified.
## -------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(rjson)
})

filterings <- c("filteredExpr10", "filteredHVG10", "filteredM3Drop10")
datasets <- c("Kumar", "Trapnell", "Koh", "Zhengmix4eq", "Zhengmix4uneq", 
              "Zhengmix8eq", "SimKumar4easy", "SimKumar4hard", "SimKumar8hard",
              "KohTCC", "KumarTCC", "TrapnellTCC")

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
  write(toJSON(list(range_clusters = 2:10)),
        file = paste0("parameter_settings/sce_", f, "_KumarTCC.json"))
  write(toJSON(list(range_clusters = 2:10)), 
        file = paste0("parameter_settings/sce_", f, "_TrapnellTCC.json"))
  write(toJSON(list(range_clusters = 2:15)), 
        file = paste0("parameter_settings/sce_", f, "_KohTCC.json"))
}

## SAFE parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/SAFE.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SAFE.json"))
    
  }
}

## CIDR parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/CIDR.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    if (f == "filteredExpr10" && d == "Zhengmix4uneq") {
      write(toJSON(list(range_clusters = 2:9)),
            file = paste0("parameter_settings/sce_", f, "_", d, "_CIDR.json"))
    } else {
      write(toJSON(list()), 
            file = paste0("parameter_settings/sce_", f, "_", d, "_CIDR.json"))
    }
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
write(toJSON(list(range_resolutions = seq(0.3, 1.5, by = 0.1),
                  min.cells = 0, min.genes = 0, dims.use = 1:30)), 
      file = "parameter_settings/Seurat.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    if (f == "filteredExpr10" && d == "Zhengmix8eq") {
      write(toJSON(list(range_resolutions = c(0.01, 0.3, 0.35, 0.4, 0.5, 1.3, 1.4, 1.6, 1.7, 1.8))),
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    } else if (f == "filteredHVG10" && d == "Zhengmix8eq") {
      write(toJSON(list(range_resolutions = c(0.1, 0.3, 1, 1.1, 1.3, 1.4, 1.6, 2))),
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    } else if (f == "filteredM3Drop10" && d == "Zhengmix8eq") {
      write(toJSON(list(range_resolutions = c(0.1, 0.7, 1, 1.3, 1.6, 2.5, 3))),
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    } else if (d %in% c("Zhengmix4eq", "Zhengmix4uneq")) {
      write(toJSON(list(range_resolutions = c(0.05, 0.1, 0.2, seq(0.3, 1.5, by = 0.1)))),
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    } else if (f == "filteredM3Drop10" && d == "SimKumar8hard") {
      write(toJSON(list(range_resolutions = sort(c(1.13, 1.15, 1.18, seq(0.3, 1.5, by = 0.1))))),
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    } else if (d %in% c("Koh", "KohTCC")) {
      if (d == "Koh" && f == "filteredM3Drop10") {
        write(toJSON(list(range_resolutions = seq(0.3, 3.6, by = 0.1))),
              file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
      } else if (d == "KohTCC" && f == "filteredExpr10") {
        write(toJSON(list(range_resolutions = sort(c(1.73, 1.76, seq(0.3, 2.1, by = 0.1))))),
              file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
      } else {
        write(toJSON(list(range_resolutions = seq(0.3, 2.1, by = 0.1))),
              file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
      }
    } else {
      write(toJSON(list()), 
            file = paste0("parameter_settings/sce_", f, "_", d, "_Seurat.json"))
    }
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

# ## Linnorm parameters
# ## -------------------------------------------------------------------------- ##
# ## General
# write(toJSON(list()), file = "parameter_settings/Linnorm.json")
# 
# ## Dataset-specific
# for (f in filterings) {
#   for (d in setdiff(datasets, c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq"))) {
#     write(toJSON(list(minNonZeroPortion = 0.75, num_PC = 3)), 
#           file = paste0("parameter_settings/sce_", f, "_", d, "_Linnorm.json"))
#   }
# }
# 
# for (f in filterings) {
#   for (d in intersect(datasets, c("Zhengmix4eq", "Zhengmix4uneq", "Zhengmix8eq"))) {
#     write(toJSON(list(minNonZeroPortion = 0.1, num_PC = 3)), 
#           file = paste0("parameter_settings/sce_", f, "_", d, "_Linnorm.json"))
#   }
# }

## RaceID parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/RaceID.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(mintotal = 1, minexpr = 0, minnumber = 1, maxexpr = Inf)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_RaceID.json"))
  }
}

## RaceID2 parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list()), file = "parameter_settings/RaceID2.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list(mintotal = 1, minexpr = 0, minnumber = 1, maxexpr = Inf)), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_RaceID2.json"))
  }
}

## SC3 parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(pct_dropout_min = 0, pct_dropout_max = 100, 
                  gene_filter = FALSE)), file = "parameter_settings/SC3.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_SC3.json"))
  }
}

## SC3svm parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(pct_dropout_min = 0, pct_dropout_max = 100,
                  gene_filter = FALSE)), file = "parameter_settings/SC3svm.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
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
    if (d %in% c("Kumar", "KumarTCC", "Trapnell", "TrapnellTCC")) {
      write(toJSON(list(xdim = 5, ydim = 5)), 
            file = paste0("parameter_settings/sce_", f, "_", d, "_FlowSOM.json"))
    } else if (d %in% c("Koh", "KohTCC", "SimKumar4easy", "SimKumar4hard", "SimKumar8hard")) {
      write(toJSON(list(list(xdim = 8, ydim = 8))), 
            file = paste0("parameter_settings/sce_", f, "_", d, "_FlowSOM.json"))
    } else {
      write(toJSON(list(xdim = 15, ydim = 15)), 
            file = paste0("parameter_settings/sce_", f, "_", d, "_FlowSOM.json"))
    }
  }
}

## PCAKmeans parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 30)), file = "parameter_settings/PCAKmeans.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_PCAKmeans.json"))
  }
}

## ascend parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 30)), file = "parameter_settings/ascend.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_ascend.json"))
  }
}

## PCAHC parameters
## -------------------------------------------------------------------------- ##
## General
write(toJSON(list(nPC = 30)), file = "parameter_settings/PCAHC.json")

## Dataset-specific
for (f in filterings) {
  for (d in datasets) {
    write(toJSON(list()), 
          file = paste0("parameter_settings/sce_", f, "_", d, "_PCAHC.json"))
  }
}
