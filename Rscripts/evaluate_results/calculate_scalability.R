args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(scefile)  ## Name of .rds file containing a SingleCellExperiment object
print(method)  ## Clustering method
print(outrds)  ## Name of .rds file where results will be written

suppressPackageStartupMessages({
  library(rjson)
  library(SingleCellExperiment)
})
source(paste0("Rscripts/clustering/apply_", method, ".R"))

## Read data
sce0 <- readRDS(scefile)

## Define subsampling sizes
nbr_cells <- round(c(0.01, 0.05, seq(0.1, 1, length.out = 7)) * ncol(sce0))

## Run clustering
set.seed(1234)
scalability <- do.call(rbind, lapply(nbr_cells, function(nc) {  ## for each number of cells
  cat(paste0("nbr cells = ", nc, "\n"))
  
  ## Subsample
  sce <- sce0[, sample(x = seq_len(ncol(sce0)), size = nc, replace = FALSE)]
  
  ## Load parameter files. General dataset and method parameters as well as
  ## dataset/method-specific parameters
  params <- c(fromJSON(file = paste0("parameter_settings/", 
                                     gsub("\\.rds$", "_", basename(scefile)), method, ".json")),
              fromJSON(file = paste0("parameter_settings/", method, ".json")), 
              fromJSON(file = paste0("parameter_settings/", 
                                     gsub("\\.rds$", ".json", basename(scefile))))
  )
  ## If any parameter is repeated, take the most specific
  if (any(duplicated(names(params)))) {
    warning("Possibly conflicting settings")
    params <- params[!duplicated(names(params))]
  }
  
  ## Run clustering
  if (method == "Seurat") {
    res <- get(paste0("apply_", method))(sce = sce, params = params, resolution = 1)
  } else {
    res <- get(paste0("apply_", method))(sce = sce, params = params, 
                                         k = length(unique(sce$phenoid)))
  }
  
  ## Tabulate timing
  data.frame(dataset = gsub("\\.rds$", "", basename(scefile)), 
             method = method,
             nbrcells = nc,
             k = length(unique(res$cluster)),
             user.self = res$st[["user.self"]],
             sys.self = res$st[["sys.self"]],
             user.child = res$st[["user.child"]],
             sys.child = res$st[["sys.child"]],
             elapsed = res$st[["elapsed"]],
             stringsAsFactors = FALSE, row.names = NULL)
}))

## Save results
saveRDS(scalability, file = outrds)

gc()
sessionInfo()