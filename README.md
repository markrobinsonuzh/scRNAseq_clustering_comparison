# scRNAseq clustering comparison
A study to compare methods for clustering of scRNA-seq data

## Instructions
- Clone this repository
- Modify the `R` and `Rscript` variables in the top of `Makefile` to point to the version of `R` that you want to use.
- Install the necessary `R` packages and dependencies:

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("RColorBrewer", "ggplot2", "dplyr", "tidyr", "VCCRI/CIDR", 
           "flowCore", "FlowSOM", "Rtsne", "scater", "SC3", 
           "Seurat", "TSCAN", "IMB-Computational-Genomics-Lab/ascend",
           "BiocParallel", "class", "SingleCellExperiment",
           "JustinaZ/pcaReduce", "clue", "rjson", 
           "plyr", "purrr", "reshape2", "parallel", "cowplot", 
           "pheatmap", "ggthemes", "viridis", "data.tree", "ggtree",
           "splatter", "cluster", "scran", "tibble", "grid", 
           "gridExtra", "VennDiagram", "ggalluvial", "M3Drop",
           "MultiAssayExperiment", "rjson", "hadley/multidplyr",
           "mclust", "ape", "clusterExperiment"))
```

- Run `make setup` to set up the directory structure and download the raw data (total size approximately 2.5GB).
- Run `make` to perform the analysis. Note that this will take a lot of time to run through completely, and you may want to run it initially with a subset of the methods (modify `include_methods.mk`) and/or data sets (modify `include_datasets.mk`). 

## Important!
- Only run `make setup` once, at the start of the analysis. Running it again will cause all the data sets to be updated, and consequently all results will be out of date and will need to be rerun.

## Adding a data set
- To add a data set to the comparison, construct a script that will generate a `SingleCellExperiment` object for each of the filterings, and put it in the corresponding folder in `data`.
- Then add the name of the data set to `include_datasets.mk`.

## Adding a method
- To add a clustering method to the comparison, construct a script containing a function that defines how to apply the method to a data set (see scripts in `Rscripts/clustering` for examples).
- Then add the name of the method to `include_methods.mk`.
