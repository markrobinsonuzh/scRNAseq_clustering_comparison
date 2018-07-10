## -------------------------------------------------------------------------- ##
## Ensemble between methods, per run, all k
## -------------------------------------------------------------------------- ##
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(clusteringsummary)
print(ncores)
print(outrds)

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(clue)
  library(reshape2)
  library(parallel)
})

## Load clustering summary
df <- readRDS(file = clusteringsummary)

## Sample Seurat results, keep only one resolution per k 
set.seed(42)
df_sub <- df %>% dplyr::filter(!(method %in% c("Seurat"))) 
df_seurat <- df %>% dplyr::filter(method == "Seurat") %>%
  dplyr::group_by(dataset, method, k) %>% 
  dplyr::filter(resolution == sample(unique(resolution), 1)) %>%
  dplyr::ungroup()
df_sub <- plyr::rbind.fill(df_sub, df_seurat)

## Help function for computing ensemble clusters
## df: dataframe with variables method, k, dataset, cell, trueclass, run
## methods: character vector with method names for ensemble
## output: tibble longformat with method cluster, ensemble cluster, trueclass, dataset
## ------------------------------------------------------------------
helper_ensemble <- function(methods, df) {
  ## Process each dataset separately
  l <- vector("list", length(unique(df$dataset)))
  names(l) <- unique(df$dataset)
  for (i in unique(df$dataset)) {
    print(paste(i, paste(methods, collapse = "+")))
    
    ## Extract only results for the current dataset
    res <- df %>% dplyr::filter(dataset == i) 
    
    ## Process each k separately
    combined_l <- vector("list", length(unique(res$k)))
    names(combined_l) <- unique(res$k)
    
    for (u in unique(res$k)) {
      ## Extract only results for the current k
      res.k <- res %>% dplyr::filter(k == u)
      
      ## Skip if results for one method not exist
      if (sum(unique(res.k$method) %in% methods) != length(methods)) {
        next
      } else {
        ## Wide format (each method/run combination in a separate column)
        res.w <- dcast(res.k %>% filter(method %in% methods), 
                       trueclass + cell ~ method + run, 
                       value.var = "cluster")
        res2 <- res.w %>% dplyr::select(-trueclass) %>% 
          tibble::column_to_rownames("cell") %>% as.matrix
        ## If all values are NA, replace them with 0
        if (all(is.na(res2))) {
          res2[] <- 0
        } 
        
        ## Process each run separately
        runs <- unique(res.k$run)
        m <- matrix(NA, nrow = nrow(res2), ncol = length(runs))
        rownames(m) <- rownames(res2)
        colnames(m) <- runs
        for (j in seq_along(runs)) {
          ## Extract only results for the current run
          run <- grep(paste0("_", runs[j], "$"), colnames(res2))
          res3 <- res2[, run] 
          re <- res3 %>% plyr::alply(2, clue::as.cl_partition) # partition
          
          ## Replace NA cluster assignment by zeros
          re <- lapply(re, function(x) { 
            x <- as.matrix(x$.Data)
            x[is.na(x)] <- 0
            x
          })
          
          ## Make ensemble and generate consensus clustering
          re <- clue::as.cl_ensemble(re)
          if (all(sapply(re, length) == 0)) {
            m[, runs[j]] <- NA
          } else {
            re <- clue::cl_consensus(re, method = "SE", control = list(nruns = 50)) # no NAS!
            clusters <- clue::cl_class_ids(re)
            m[, runs[j]] <- clusters[rownames(m)]
          }
        }

        out <- data.frame(res.w, stringsAsFactors = FALSE) %>%
          dplyr::mutate(dataset = i, k = u, method = paste(methods, collapse = ".")) %>%
          dplyr::left_join(data.frame(m, stringsAsFactors = FALSE, check.names = FALSE) %>% 
                           tibble::rownames_to_column("cell"), by = "cell")
        combined_l[[as.character(u)]] <- out
      }
    }
    
    l[[i]] <- plyr::rbind.fill(combined_l)
    
  }
  res.df <-  plyr::rbind.fill(l)
  res.df <- reshape2::melt(data = res.df, 
                           id.vars = c("dataset", "trueclass", "cell", "method", "k"), 
                           measure.vars = as.character(runs),
                           value.name = "cons_cluster", variable.name = "run")
  return(res.df)
}

## Construct all combinations of two methods
comb.ensmbl <- list(l1 = unique(df_sub$method),
                    l2 = unique(df_sub$method)) %>% 
  purrr::cross() %>% purrr::map((paste))

names(comb.ensmbl) <- sapply(comb.ensmbl, function(x) paste0(x, collapse = "."))

## Remove pairings of a method with itself
comb.ensmbl <- comb.ensmbl %>% purrr::discard(function(x) x[1] == x[2])

## Run 
out <- mclapply(comb.ensmbl, helper_ensemble, df = df_sub, 
                mc.preschedule = FALSE, mc.cores = ncores)
out <- plyr::rbind.fill(out)

## Save
saveRDS(out, file = outrds)

date()
sessionInfo()
