# Overview Table for datasets
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

filterings <- strsplit(filterings, ",")[[1]]
datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets

print(filterings)
print(pctkeep)
print(datasets)
print(datadir)
print(silhouetterds)
print(outcsv)

suppressPackageStartupMessages({
  library(scater)
  library(parallel)
  library(dplyr)
})

datasets_filtered <- c(outer(paste0(datadir, "/sce_filtered", filterings, pctkeep, "/sce_filtered", 
                                    filterings, pctkeep), 
                             paste0(datasets, ".rds"),
                             FUN = paste, sep = "_"))
names(datasets_filtered) <- gsub("\\.rds", "", basename(datasets_filtered))

datasets_full <- paste0(datadir, "/sce_full/sce_full_", datasets, ".rds")
names(datasets_full) <- gsub("\\.rds", "", basename(datasets_full))

all_data <- lapply(c(datasets_filtered, datasets_full), function(x) {
  readRDS(x)
})

## number of cells, number of features per dataset, median counts per cell,
## median genes per cell, number of suppopulations
tbl <- do.call(rbind, lapply(names(all_data), function(nm) {
  x <- all_data[[nm]]
  data.frame(dataset = nm,
             ncells = ncol(x),
             ngenes = nrow(x),
             mediancount = median(colSums(counts(x)), na.rm = TRUE),
             mediangenes = median(colSums(counts(x) != 0), na.rm = TRUE),
             npop = length(unique(colData(x)$phenoid)),
             stringsAsFactors = FALSE)
}))

## avg.silhoutte width
silh <- readRDS(file = silhouetterds)
avg.s <- sapply(silh, function(x) x$avg.width)

## table
tbl <- tbl %>% dplyr::left_join(data.frame(dataset = names(avg.s),
                                           avg_silh = avg.s,
                                           stringsAsFactors = FALSE),
                                by = "dataset")

write.csv(tbl, file = outcsv)

date()
sessionInfo()