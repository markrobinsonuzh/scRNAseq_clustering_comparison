# NOTE! This scripts has been modified from the one downloaded from
# https://yunliweb.its.unc.edu//safe/download.php, in order to set the number of
# cores to use for the clustering. All modifications are indicated in the code.
# 
# !/usr/bin/env Rscript
# 11/10/2017

# @title Single-cell Aggregated clustering From Ensemble (SAFE)
#
# @description This function performs single-cell clustering using four state-of-the-art methods, 
# SC3, CIDR, Seurat and tSNE+kmeans.
# 
# @param inputTags a G*N matrix with G genes and N cells.
# @param datatype defines the type of data, which could be "count", "CPM", "RPKM" and "FPKM". 
# Default is "count".
# @param SC3 a boolean variable that defines whether to cluster cells using SC3 method. 
# Default is "TRUE".
# @param gene_filter a boolean variable defines whether to perform gene filtering 
# before SC3 clustering, when \code{SC3 = TRUE}.
# @param svm_max, if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.
# @param CIDR a boolean parameter that defines whether to cluster cells using CIDR method. 
# Default is "TRUE".
# @param Seurat is a boolean variable that defines whether to cluster cells using Seurat method. 
# Default is "TRUE".
# @param nPC defines the number of princple compoents used in Seurat clustering, when \code{Seurat = TRUE}. 
# Default value is esimated by \code{nPC} of \code{CIDR}.
# @param resolution defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}.
# @param seurat_min_cell defines the mimimum number of cells in input dataset below which 
# \code{resolution} is set to 1.2, when \code{Seurat = TRUE}.
# @param resolution_min defines the resolution used in Seurat clustering for small dataset, 
# when \code{Seurat == TRUE} and cell number of input file < \code{seurat_min_cell}.
# @param tSNE is a boolean variable that defines whether to cluster cells using tSNE method.
# Default is "TRUE".
# @param var_genes defines the number of variable genes used by tSNE analysis, when \code{tSNE = TRUE}.
# @param SEED sets the seed of the random number generator. Setting the seed to a fixed value can 
# produce reproducible clustering results. 
#
# @return a matrix of indiviudal clustering results is returned.
#
# @author Yuchen Yang <yyuchen@email.unc.edu>, Ruth Huh <rhuh@live.unc.edu>, 
# Houston Culpepper <hculpepp@live.unc.edu>, Yun Li <yunli@med.unc.edu>
# @citation Yuchen Yang, Ruth Huh, Houston Culpepper, Yun Li. SAFE (Single-cell Aggregated clustering From Ensemble): Cluster ensemble for single-cell RNA-seq data. 2017

## Modification: add "n_cores" argument
individual_clustering <- function(inputTags, datatype = "count", 
                                  SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, 
                                  Seurat = TRUE, nPC = NULL, resolution = 0.9, seurat_min_cell = 200, resolution_min = 1.2, 
                                  tSNE = TRUE, var_genes = NULL, SEED = 1, n_cores){
  
  ## Modification: Only set seed to SEED if not NULL, otherwise generate a random seed
  if (is.null(SEED)) {
    SEED <- round(1e6*runif(1))
  }
  set.seed(SEED)
  
  cluster_results <- NULL
  inputTags = as.matrix(inputTags)
  
  
  ##### SC3
  if(SC3 == TRUE){
    library(scater)  #https://github.com/davismcc/scater
    library(SC3)  #https://github.com/hemberg-lab/SC3
    exp_cell_exprs <- NULL
    sc3OUTPUT <- NULL
    
    # cell expression
    if (datatype == "count") {
      ### If the input data is original count data, it would be normalized by the total cound number and then log2 transformed
      exp_cell_exprs <- SingleCellExperiment(assays = list(counts = as.matrix(inputTags)))
      normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
      logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
      #exp_cell_exprs <- inputTags
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
      ### If the input data is CPM, FPKM, RPKM or TPM, it would be log2 transformed
      exp_cell_exprs <- SingleCellExperiment(assays = list(normcounts = as.matrix(inputTags)))
      logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
    }
    
    rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
    exp_cell_exprs <- exp_cell_exprs[!duplicated(rownames(exp_cell_exprs)), ]
    
    ### Estimating optimal number of clustering
    exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
    optimal_K = metadata(exp_cell_exprs)$sc3$k_estimation
    
    ### Clustering by SC3 at the optimal K
    ## Modification: set n_cores and pct_dropout_min/max
    if (ncol(inputTags) < svm_num_cells){
      exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, rand_seed = SEED, n_cores = n_cores, pct_dropout_min = 0, pct_dropout_max = 100)
    }
    
    ### Runing SVM
    if (ncol(inputTags) >= svm_num_cells){
      ## Modification: set n_cores and pct_dropout_min/max
      exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, 
                            svm_max = svm_num_cells, svm_num_cells = svm_num_cells, rand_seed = SEED,
                            n_cores = n_cores, pct_dropout_min = 0, pct_dropout_max = 100)
      exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
    }
    
    ### Exporting SC3 results
    p_Data <- colData(exp_cell_exprs)
    col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
    sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
    cluster_results = rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = T))
  }
  
  
  ##### CIDR
  if(CIDR == TRUE){
    library(cidr)  #https://github.com/VCCRI/CIDR
    cidrOUTPUT <- NULL
    
    cidrOUTPUT <- scDataConstructor(inputTags)
    cidrOUTPUT <- determineDropoutCandidates(cidrOUTPUT)
    cidrOUTPUT <- wThreshold(cidrOUTPUT)
    ## Modification: set n_cores
    cidrOUTPUT <- scDissim(cidrOUTPUT, threads = n_cores)
    cidrOUTPUT <- scPCA(cidrOUTPUT)
    cidrOUTPUT <- nPC(cidrOUTPUT)
    
    ### nPC is not necessarily the same between CIDR and Seurat. Here we use the same nPC value for both two methods for convenience
    if(!is.null(nPC)) {
      cidrOUTPUT@nPC <- nPC
    } else {
      nPC <- cidrOUTPUT@nPC
    }
    
    ### Clustering by CIDR.
    # The optimal clustering number is determined automatically
    cidrOUTPUT <- scCluster(cidrOUTPUT, nPC = nPC)
    cluster_results = rbind(cluster_results, matrix(c(cidrOUTPUT@clusters), nrow = 1, byrow = T))
  }
  
  
  ##### Seurat
  if (Seurat == TRUE){
    library(dplyr)
    library(Seurat)  #http://satijalab.org/seurat/get_started.html
    library(Matrix)
    seuratOUTPUT <- NULL
    
    # Initialize the Seurat object with the raw data (non-normalized data)
    # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
    ## Modification: remove cell and gene filtering
    seuratOUTPUT <- CreateSeuratObject(raw.data = inputTags, min.cells = 0, min.genes = 0, project = "single-cell clustering")
    
    # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
    if (datatype == "count"){
      seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM"){
      raw.data <- GetAssayData(object = seuratOUTPUT, slot = "raw.data")
      normalized.data <- log(raw.data+1)
      colnames(x = normalized.data) <- colnames(x = raw.data)
      rownames(x = normalized.data) <- rownames(x = raw.data)
      seuratOUTPUT <- SetAssayData(object = seuratOUTPUT, assay.type = "RNA",slot = "data", new.data = normalized.data)
    }
    
    # Detection of variable genes across the single cells
    ## Modification: remove search for variable genes since we provide pre-filtered data
    # seuratOUTPUT = FindVariableGenes(object = seuratOUTPUT, mean.function = ExpMean, dispersion.function = LogVMR, 
    #                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    
    # Regress out unwanted sources of variation
    seuratOUTPUT <- ScaleData(object = seuratOUTPUT, vars.to.regress = c("nUMI"))
    
    ### Perform linear dimensional reduction
    ## Modification: use all genes
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, pc.genes = rownames(seuratOUTPUT@data), do.print = FALSE)
    
    if (length(inputTags[1,]) >= seurat_min_cell){
      ### Determine statistically significant principal components
      # NOTE: This process can take a long time for big datasets, comment out for expediency.
      # More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
      # Here we chooes the same number of PCs used in CIDR
      
      ### Clustering the cells by Seurat
      seuratOUTPUT <- Seurat::FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC, algorithm = 3,
                                           resolution = resolution, print.output = F, random.seed = SEED)
    } else {
      resolution <- resolution_min
      seuratOUTPUT <- Seurat::FindClusters(object = seuratOUTPUT, reduction.type = "pca", dims.use = 1:nPC, algorithm = 3,
                                           resolution = resolution_min, print.output = F, random.seed = SEED)
    }
    
    ### Complementing the missing data
    cells_dropout <- NULL
    num_genes <- colSums(seuratOUTPUT@raw.data > 0)
    cells_dropout <- names(num_genes[which(num_genes <= 0)])
    if (length(cells_dropout != 0)){
      seurat_output <- matrix(NA, ncol = ncol(seuratOUTPUT@raw.data), byrow = T)
      colnames(seurat_output) <- colnames(seuratOUTPUT@raw.data)
      seurat_retained <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
      colnames(seurat_retained) <- colnames(seuratOUTPUT@data)
      for (i in 1:ncol(seurat_retained)){
        seurat_output[1,colnames(seurat_retained)[i]] <- seurat_retained[1,colnames(seurat_retained)[i]]
      }
    } else {
      seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@ident)))
    }
    
    cluster_results = rbind(cluster_results, matrix(c(as.numeric(seuratOUTPUT@ident)), nrow = 1, byrow = T))
  }
  
  
  ##### tSNE+kmeans
  if(tSNE == TRUE){
    library(Rtsne)
    library(ADPclust)
    input_lcpm <- NULL
    tsne_input <- NULL
    tsne_output <- NULL
    tsne_kmeansOUTPUT <- NULL
    adpOUTPUT <- NULL
    ### Data tranformation
    if (datatype == "count") {
      ### If the input data is original count data or CPM, it would be tranformed to CPM
      input_lcpm <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
      tsne_input <- input_lcpm
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
      ### If the input data is FPKM or RPKM, we use the transformed TPM data generated before as the input
      tsne_input <- log2(inputTags + 1)
    }
    
    ### Dimensionality reduction by Rtsne
    if(length(inputTags[1,]) >= 200) {
      perplexity = 30
    } else {
      perplexity = 10
    }
    
    if (is.null(var_genes)){
      tsne_output <- Rtsne(t(tsne_input), perplexity = perplexity)
    } else{
      se_genes = rep(NA, nrow(tsne_input))
      for (i in 1:nrow(tsne_input)){
        se_genes[i] = sqrt(var(tsne_input[i,])/length(tsne_input[i,]))
      }
      decreasing_rank = order(se_genes, decreasing = T)
      
      tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes],]), perplexity = perplexity)
    }
    
    ### Determining the optimal cluster number (k) by ADPclust
    
    mim_K = 2
    max_K = max(optimal_K, cidrOUTPUT@nCluster, max(as.numeric(seuratOUTPUT@ident)))
    
    adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto", nclust = mim_K:max_K)
    
    ### Clustering the cells by kmeans
    tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)
    cluster_results = rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = T))
  }
  
  return(cluster_results)
}
