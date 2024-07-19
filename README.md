[![DOI](https://zenodo.org/badge/98410072.svg)](https://zenodo.org/badge/latestdoi/98410072)

# scRNAseq clustering comparison
This repository contains the code for our study comparing methods for unsupervised clustering of scRNA-seq data:

- Duò A, Robinson MD and Soneson C: [A systematic performance evaluation of clustering methods for single-cell RNA-seq data](https://f1000research.com/articles/7-1141/v2). F1000Research 7:1141 (2018).

Please note that the purpose of this repository is to provide a record of the exact code we used for the analyses in this paper. Thus, to run it you need the same versions of software and packages as we used (see below), and we can not guarantee that it will execute properly or generate the same results with new versions of any of the software. If you want to add a method to the comparison, we recommend using the [`DuoClustering2018`](https://bioconductor.org/packages/DuoClustering2018/) package to retrieve the data and our results (see below). 

Updates to the code aimed at addressing differences between the package versions used in our analyses and more recent versions are included in the `updates` branch. Please note that this is intended only as a help to users and that these changes have not been tested, and we can not guarantee that the code in this branch will run without errors, nor that it contains all necessary updates. 

## Downloading the data and results directly
The data and clustering results from the paper are available from the [`DuoClustering2018`](https://bioconductor.org/packages/DuoClustering2018/) Bioconductor package. This is the recommended way to access and extend our benchmark (rather than rerunning all the code in this repository), and the package contains vignettes to show, e.g., how to add new methods. Install the package using the `BiocManager` CRAN package (note that you need the development version of Bioconductor):

```
install.packages("BiocManager")
BiocManager::install("DuoClustering2018")
```

The unfiltered and filtered data sets, as well as all the clustering results from v1 of the paper, can also be downloaded as a compressed archive from [here](https://zenodo.org/records/12772318) (4.93GB). 


## Instructions for running the code
If you want to reproduce the results of our study, you need to go through the steps below. Please note that running all the analyses will take a considerable amount of time (and that the preprocessed data and clustering output is available as described above).

- Clone this repository
- Modify the `R` and `Rscript` variables in the top of the [Makefile](Makefile) to point to the version of `R` that you want to use.
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
                       "plyr", "purrr", "reshape2", "parallel", 
                       "pheatmap", "ggthemes", "viridis", "data.tree", 
                       "ggtree", "cowplot", "grid", "M3Drop",
                       "splatter", "cluster", "scran", "tibble", 
                       "gridExtra", "VennDiagram", "ggalluvial", 
                       "MultiAssayExperiment", "rjson", 
                       "hadley/multidplyr",
                       "mclust", "ape", "clusterExperiment", 
                       "csoneson/countsimQC", "monocle"))
```

- If you want to include RaceID2, download the `RaceID2_StemID_class.R` script from [https://github.com/dgrun/StemID](https://github.com/dgrun/StemID) (we used the version from March 3, 2017) and place it in the [Rscripts/clustering](Rscripts/clustering) directory. 
- Run `make setup` to set up the directory structure and download the raw data (total size approximately 2.7GB).
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

## countsimQC report 
The [countsimQC](https://github.com/csoneson/countsimQC) report, comparing the characteristics of the simulated data sets to the underlying real data set, can be found [here](http://imlspenticton.uzh.ch/robinson_lab/DuoClustering2018/Kumar_countsimQC.html).

## R package versions
The following R packages were used for the clustering and evaluation in our paper:

```
ADPclust_0.7
amap_0.8-16
ape_5.1
ascend_0.5.0
bindrcpp_0.2.2
Biobase_2.40.0
BiocGenerics_0.26.0
BiocParallel_1.14.1
bit_1.1-14
bit64_0.9-7
cidr_0.1.5
class_7.3-14
clue_0.3-55
cluster_2.0.7-1
clusterExperiment_2.0.2
cowplot_0.9.2
data.table_1.11.4
data.tree_0.7.5
DelayedArray_0.6.0
dplyr_0.7.5
flexmix_2.3-14
flowCore_1.46.1
FlowSOM_1.12.0
fpc_2.1-11
futile.logger_1.4.3
GenomeInfoDb_1.16.0
GenomicRanges_1.32.3
ggalluvial_0.6.0
ggplot2_2.2.1
ggthemes_3.5.0
ggtree_1.12.0
gridExtra_2.3
igraph_1.2.1
IRanges_2.14.10
lattice_0.20-35
locfit_1.5-9.1
M3Drop_1.4.0
MASS_7.3-50
Matrix_1.2-14
matrixStats_0.53.1
mclust_5.4
mnormt_1.5-5
multidplyr_0.0.0.9000
numDeriv_2016.8-1
pcaMethods_1.72.0
pcaReduce_1.0
permute_0.9-4
pheatmap_1.0.10
plyr_1.8.4
purrr_0.2.5
RColorBrewer_1.1-2
reshape2_1.4.3
rjson_0.2.19
Rtsne_0.13
S4Vectors_0.18.2
SC3_1.8.0
scater_1.8.0
scran_1.6.9
Seurat_2.3.1
SingleCellExperiment_1.2.0
stringr_1.3.1
SummarizedExperiment_1.10.1
tidyr_0.8.1
TSCAN_1.18.0
tsne_0.1-3
vegan_2.5-2
VennDiagram_1.6.20
viridis_0.5.1
viridisLite_0.3.0
```

## Disclaimer
All methods included in the evaluation except [SAFE](https://www.biorxiv.org/content/early/2018/03/28/215723) and [RaceID2]() are implemented as R packages (for version, see above). To ensure full reproducibility of our analyses, we have included the version of the SAFE-clustering code that was used for our evaluations in the repository. This includes the [gpmetis](gpmetis) and [shmetis](shmetis) executables, as well as the folder [Rscripts/clustering/SAFE_2.1_Linux](Rscripts/clustering/SAFE_2.1_Linux). These files were obtained from [https://yunliweb.its.unc.edu//safe/](https://yunliweb.its.unc.edu//safe/) on July 2, 2018. Minor modifications were made in the [Rscripts/clustering/SAFE_2.1_Linux/SAFE_modified.R](Rscripts/clustering/SAFE_2.1_Linux/SAFE_modified.R) and [Rscripts/clustering/SAFE_2.1_Linux/individual_clustering_modified.R](Rscripts/clustering/SAFE_2.1_Linux/individual_clustering_modified.R) scripts (clearly marked in the respective scripts), in order to allow setting the number of cores for execution, and to disable the internal cell/gene filtering for SC3 and Seurat, since we provide pre-filtered data in our evaluation. 

In order to run RaceID2, you need to download the `RaceID2_StemID_class.R` script from [https://github.com/dgrun/StemID](https://github.com/dgrun/StemID) (we used the version from March 3, 2017) and place it in the [Rscripts/clustering](Rscripts/clustering) directory. 
