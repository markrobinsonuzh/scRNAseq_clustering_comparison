# scRNAseq clustering comparison
A study to compare methods for clustering of scRNA-seq data

## Instructions
- Clone this repository
- Modify the `R` and `Rscript` variables in the top of `Makefile` to point to the version of `R` that you want to use.
- Run `make setup` to set up the directory structure and download the raw data (total size approximately 2.5GB).
- Run `make` to perform the analysis.

## Notes
- Only run `make setup` once, at the start of the analysis. Running it again will cause all the data sets to be updated, and consequently all results will be out of date and will need to be rerun.

## Adding a data set
- To add a data set to the comparison, construct a script that will generate a `SummarizedExperiment` object for each of the filterings, and put it in the corresponding folder in `data`.
- Then add the name of the data set to `include_datasets.mk`.

## Adding a method
- To add a clustering method to the comparison, construct a script containing a function that defines how to apply the method to a data set (see scripts in `Rscripts/clustering` for examples).
- Then add the name of the method to `include_methods.mk`.
