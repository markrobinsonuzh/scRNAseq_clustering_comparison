###########################
# Helper Functions
###########################

##########################
# load cell labels from a SingleCellexpriment class
#####################
### Input: list of SingleCellexpriment class data
### Output:  list of labels
load_labels <- function(data) {
  
  labels <- vector("list", length(data))
  names(labels) <- names(data)
  
for (i in names(data)) {
  labels[[i]] <- as.character(colData(data[[i]])$phenoid)
  }
  return(labels)
}



############################################
# load SingleCellexpriment class data sets 
#################################################
### Input: list of filepaths, vector with data directory path
### Output:  list with SingleCellexpriment data sets

load_data <- function( files, DATA_DIR ) { 
  
  data <- labels <- vector("list", length(files))
  
  names(data) <-names(labels) <-  names(files)
  
  for (i in 1:length(data)){
    f <- files[[i]]
    load(f)
    data[[i]] <- res
    
  }
  return(data)
}
#####################
# for Data import
##################

plot_QC <- function(data=sceset){
  ### Frequency of top 50 genes
  print(plotQC(sceset, type = "highest-expression")) # Top 50 genes account for 18.4 % of the counts, profile relatively flat
  
  print(plotQC(sceset, type = "find-pcs", variable = "total_features")) # PC1 shows high correlation with total number observerd features per cell
  print(plotQC(sceset, type = "find-pcs", variable = "total_features", plot_type = "pairs-pcs")) # PC1 and PC75, no correlation between other PCs
  
  
  # plot the variance explained for some explanatory vars.
  expl_vars <- c("source_name_ch1", "log10_total_counts", "log10_total_features", "pct_dropout",
                 "pct_counts_top_200_features", ifelse(withcontrols, "log10_counts_feature_controls", NA),
                 ifelse(withcontrols, "pct_counts_feature_controls", NA))
  
  print(plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)])) # the cell type and culture condition is the main source of variation
  
  # Pairs plot, showing correlation between some explanatory variables
  
  print(plotQC(sceset, type = "explanatory-variables", variables = expl_vars[!is.na(expl_vars)], method = "pairs")) # negativ association between dropout and total features, clear association between total counts and total features 
  
  
  # PLot 
  print(plotPhenoData(sceset, aes_string(x = "log10_total_counts", y = "total_features",
                                   colour = paste("source_name_ch1", collapse = "."))) + 
    guides(colour = guide_legend( byrow = TRUE)) + theme(legend.position = "bottom")
  )
  
  
  if (withcontrols) {
   print( plotPhenoData(sceset, aes_string(x = "total_features", y = "pct_counts_feature_controls", 
                                     color = paste("source_name_ch1", collapse = "."))) + 
      guides(colour = guide_legend( byrow = TRUE)) + theme(legend.position = "bottom")
   )
  }
  
  
  print(plotQC(sceset, type = "exprs-freq-vs-mean"))
}


###########################
# For ARI script
###########################

# Store results of ARI in textfiles, per method, per dataset
store.ari <- function(ari.cluster, DATA_DIR, METHOD) {
  file_results <- paste0(DATA_DIR,METHOD ,"_ARI_", names(ari.cluster), ".txt")
  names(file_results) <- names(ari.cluster)
  
  for (i in names(ari.cluster)) {
    res_i <- ari.cluster[[i]]
    write.table(res_i, file = file_results[i], row.names = FALSE, quote = FALSE, sep = "\t")
  }
}

#### read in labels for clustering evaluation

read.labels <- function(files_labels) {
  
  labels <- vector("list", length(files_labels))
  names(labels) <- names(files_labels) 
  
  for (i in names(files_labels)) {
    lab_i <-  read.csv(files_labels[[i]], sep="\t")
    labels[[i]] <- as.vector(unlist(lab_i))
    
  }
  return(labels)
}
#### read in cluster for clustering evaluation

read.cluster <- function(files_clusters){
  cluster <- vector("list", length(files_clusters))
  names(cluster) <- names(files_clusters) 
  
  for (i in 1:length(files_clusters)) {
    clus_i <-  read.csv(files_clusters[[i]], sep="\t")
    cluster[[i]] <- as.vector(unlist(clus_i))
  }
  return(cluster)
}


###############################
#  PLOTS
###############################

#### MEAN SD smoothScatter PLOT function
meansdplot <- function (data,rank,ylim, xlim, title){
  mean = apply(data , 1,mean)
  sd = apply(data, 1,sd)
  runmed=runmed(sd, k=3)
  if (rank==FALSE){
    smoothScatter(mean, sd, main=title, xlab="mean", ylab="sd", pch=".",ylim = ylim , xlim=xlim)
  }
  else {
    rankmean=rank(mean)
    smoothScatter(rankmean, sd, main=title, xlab="rank(mean)", ylab="sd", pch=".",ylim = ylim , xlim=xlim)
  }
}


###############################
#  methods
###############################

######## save the clusters results to txt file

save_clusters <- function(res.cluster, dir_cluster){
  for (i in seq_len(length(res.cluster))) {
    res_i <- res.cluster[[i]]
    write.table(res_i, file = dir_cluster[i], row.names = FALSE, quote = FALSE, sep = "\t")
  }
}

######## save the system time results to txt file

save_systemtime <- function(sys.time, dir_systime){
  for (i in seq_len(length(sys.time))){
  sys_i <- sys.time[[i]]["elapsed"]
  write.table(sys_i, file=dir_systime[i], sep="\t")
  }
}

######## Determine number of cluster


