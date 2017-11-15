#####################
# SNN-Cliq
#####################
# SNN-Cliq uses shared nearest neighbour graphs to model the high dimensional data given by the expression matrix. The identification of clusters is done by merging the quasi cliques obtained by the graph.
# Paramteres supplied by the user are a threshold r which defines the connectivity of the quasi-cliques, a threshold m which defines the merging rate of the quasi-cliques. r and m is typically set to 0.7 and 0,5, respectively.
# The number nearest neighbors k has to be defined by the user, higher number of neighbors gives lower number of clusters. Working distances are euclidean but diferent can be used, defined by the distan argument.
require(scater)
source("skript/helper_files/Helper_functions.R")


# file paths
source("FILES.R")

# load data sets

data <- res.cluster <- vector("list", length(files))

names(data) <- names(res.cluster) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# create storage files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- input_matrix<- labels<- list

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}


# extract expression data
for (i in names(input_matrix)){
  input_matrix[[i]] <-exprs(data[[i]]) # use exprs slot of SCeset; log2, normalized count_lstpm
}


# RUN SNN-Cliq
# set workin directory to save shared nearest neighbourhood d matrix 
setwd("skript/run_methods/snn-cliq")
# Set parameters, k is number of nearest neighbours, r is cutoff for quasi-clique, m is threshold of overlapping quasicliques to be merged (standart threshold is half of the points)

distan <- "euclidean"
# set number of neirest neighbors, below 3 neibghors doesnt work
par.k <-  list(
  kumar2015 =ncol(data[[i]])* c(0.01,0.05,0.1,0.15),
  trapnell2014 = ncol(data[[i]])*c(0.01,0.05,0.1,0.15),
  zhengmix2016 = ncol(data[[i]])*c(0.01,0.05,0.1,0.15),
  koh2016 = ncol(data[[i]])*c(0.01,0.05,0.1,0.15)
)

par.r <- list(
  kumar2015 = 0.7,
  trapnell2014 = 0.7,
  zhengmix2016 = 0.7,
  koh2016 = 0.7
)
par.m <-  list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  zhengmix2016 = 0.5,
  koh2016 = 0.5
)
# SNN-clique function
run_snnclique <- function( input_matrix, par.k, par.m, par.r ) {
  
  for (i in names(input_matrix)){
    df.clus <-  matrix( nrow = ncol(input_matrix[[i]]), ncol = length(par.k[[i]]) ) 
    for ( j in seq_len( length(par.k[[i]]) ) ){
      # construct a graph 
      scRNA.seq.funcs::SNN(
        data = t(input_matrix[[i]]),
        outfile = "snn-cliq.txt",
        k = par.k[[i]][j],
        distance = distan
      )
      
      
      # find clusters in the graph
      
        snn.res <- 
          system(
            paste0(
              "python Cliq.py ", 
              "-i snn-cliq.txt ",
              "-o res-snn-cliq.txt ",
              "-r ", par.r[[i]],
              " -m ", par.m[[i]]
            ),
            intern = TRUE
          )
     
      #
      cat(paste(snn.res, collapse = "\n"))
      snn.res <- read.table("res-snn-cliq.txt")
      
      
      df.clus[,j]<- as.integer(snn.res[,1])
      # remove files that were created during the analysis
      #system("rm snn-cliq.txt res-snn-cliq.txt")
    }
    colnames(df.clus) <- c( paste0(par.k[[i]]) )
    res.cluster[[i]] <- df.clus
  }
  return( res.cluster )
}
# run method

res.cluster <- run_snnclique(input_matrix,par.k,par.m,par.r)
1# save clusters

setwd("~/Desktop/masterarbeit/scRNAseq_clustering_comparison")
dir_cluster <- paste0("results/SNNCliq/SNNCliq_krange_clus_", names(res.cluster), ".txt")

save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SNNCliq/SNNCliq_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("results/SNNCliq/SNNCliq_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "results/SNNCliq/session_info_SNNCliq_krange.txt")
sessionInfo()
sink()

# Appendix

