## Apply SAFE clustering

source("Rscripts/clustering/SAFE_2.1_Linux/individual_clustering_modified.R")
source("Rscripts/clustering/SAFE_2.1_Linux/SAFE_modified.R")

apply_SAFE <- function(sce, params, k) {
  tryCatch({
    dat <- counts(sce)
    st <- system.time({
      indc <- individual_clustering(inputTags = dat, datatype = "count", 
                                    SC3 = TRUE, gene_filter = FALSE, 
                                    svm_num_cells = 1000, CIDR = TRUE, 
                                    Seurat = TRUE, nPC = NULL, 
                                    resolution = 0.9, seurat_min_cell = 0, 
                                    resolution_min = 1.2, tSNE = TRUE, 
                                    var_genes = NULL, SEED = NULL,
                                    n_cores = 1)
      safe <- SAFE(cluster_results = indc, k_min = k, k_max = k, 
                   MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
                   cspc_cell_max = NULL)
      
      cluster <- safe$optimal_clustering
      names(cluster) <- colnames(dat)
    })
    ## Determine number of clusters automatically
    safe <- SAFE(cluster_results = indc, k_min = 2, k_max = NULL, 
                 MCLA = TRUE, HGPA = FALSE, CSPA = FALSE, 
                 cspc_cell_max = NULL)
    est_k <- safe$optimal_k
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = est_k)
  }, error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA), 
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
