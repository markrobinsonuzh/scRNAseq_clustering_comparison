# scRNAseq clustering comparison
This repository contains the code for our study comparing methods for unsupervised clustering of scRNA-seq data

## Instructions
If you want to reproduce the results of our study, you need to go through the steps below. Please note that running all the analyses will take a considerable amount of time.

- Clone this repository
- Modify the `R` and `Rscript` variables in the top of [Makefile](Makefile) to point to the version of `R` that you want to use.
- Install the necessary `R` packages and dependencies:

```
install.packages("BiocManager")
BiocManager::install(c("RColorBrewer", "ggplot2", "dplyr", 
                       "tidyr", "VCCRI/CIDR", 
                       "flowCore", "FlowSOM", "Rtsne", "scater", "SC3", 
                       "Seurat", "TSCAN", 
                       "IMB-Computational-Genomics-Lab/ascend",
                       "BiocParallel", "class", "SingleCellExperiment",
                       "JustinaZ/pcaReduce", "clue", "rjson", 
                       "plyr", "purrr", "reshape2", "parallel", "cowplot", 
                       "pheatmap", "ggthemes", "viridis", "data.tree", "ggtree",
                       "splatter", "cluster", "scran", "tibble", "grid", 
                       "gridExtra", "VennDiagram", "ggalluvial", "M3Drop",
                       "MultiAssayExperiment", "rjson", "hadley/multidplyr",
                       "mclust", "ape", "clusterExperiment", "csoneson/countsimQC"))
```

- Run `make setup` to set up the directory structure and download the raw data (total size approximately 2.5GB).
- Run `make` to perform the analysis. Note that this will take a lot of time to run through completely, and you may want to run it initially with a subset of the methods (modify [include_methods.mk](include_methods.mk)) and/or data sets (modify [include_datasets.mk](include_datasets.mk)). 

## Important!
- Only run `make setup` once, at the start of the analysis. Running it again will cause all the data sets to be updated, and consequently all results will be out of date and will need to be rerun.

## Adding a data set
- To add a data set to the comparison, construct a script that will generate a `SingleCellExperiment` object for each of the filterings, and put it in the corresponding folder in `data`.
- Generate `rds` files containing the parameter values to use for each method (or empty lists if there are no hyperparameters). See [Rscripts/parameter_settings/generate_parameter_settings.R](Rscripts/parameter_settings/generate_parameter_settings.R) for examples.
- Then add the name of the data set to [include_datasets.mk](include_datasets.mk).

## Adding a method
- To add a clustering method to the comparison, construct a script containing a function that defines how to apply the method to a data set (see scripts in [Rscripts/clustering](Rscripts/clustering) for examples).
- Generate `rds` files containing the parameter values to use for each data set (or empty lists if there are no hyperparameters). See [Rscripts/parameter_settings/generate_parameter_settings.R](Rscripts/parameter_settings/generate_parameter_settings.R) for examples.
- Then add the name of the method to [include_methods.mk](include_methods.mk).

## Downloading the data and results directly
The unfiltered and filtered data sets, as well as all the clustering results, can be downloaded as a compressed archive from [here](http://imlspenticton.uzh.ch/robinson_lab/DuoClustering2018/DuoClustering2018.tar.gz) (4.93GB). 

## countsimQC report 
The [countsimQC](https://github.com/csoneson/countsimQC) report, comparing the characteristics of the simulated data sets to the underlying real data set, can be found [here](http://htmlpreview.github.io/?https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/output/countsimQC/Kumar_countsimQC.html).