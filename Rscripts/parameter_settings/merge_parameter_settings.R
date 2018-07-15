args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

filterings <- strsplit(filterings, ",")[[1]]
names(filterings) <- filterings
datasets <- strsplit(datasets, ",")[[1]]
names(datasets) <- datasets
methods <- strsplit(methods, ",")[[1]]
names(methods) <- methods

print(paramdir)
print(datasets)
print(filterings)
print(methods)
print(outrds)

suppressPackageStartupMessages({
  library(rjson)
})

## Load parameter files. General dataset and method parameters as well as
## dataset/method-specific parameters
allparams <- lapply(filterings, function(f) {
  lapply(datasets, function(d) {
    lapply(methods, function(m) {
      params <- c(fromJSON(file = paste0(paramdir, "/sce_", f, "_", d, "_", m, ".json")),
                  fromJSON(file = paste0(paramdir, "/", m, ".json")), 
                  fromJSON(file = paste0(paramdir, "/sce_", f, "_", d, ".json")))
      ## If any parameter is repeated, take the most specific
      if (any(duplicated(names(params)))) {
        warning("Possibly conflicting settings")
        params <- params[!duplicated(names(params))]
      }
      params
    })
  })
})
allparams <- unlist(unlist(allparams, recursive = FALSE), recursive = FALSE)
names(allparams) <- paste0("sce_filtered", gsub("\\.", "_", names(allparams)))

saveRDS(allparams, outrds)
date()
sessionInfo()
