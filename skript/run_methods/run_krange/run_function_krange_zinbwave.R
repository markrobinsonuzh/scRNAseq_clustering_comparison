#####################
# ZINB-WaVE
#####################
# uses reduced raw count  matrix. 

require("zinbwave")
require(scater)
require(dplyr)
require(Rtsne)

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
# load cell labels
labels <- load_labels(data) 


# RUN ZINB-WaVE
par.k <- list(
  kumar2015 = c(2:10),
  trapnell2014 =  c(2:10),
  zhengmix2016 =  c(2:10),
  koh2016 =  c(2:15),
  simDataKumar = c(2:10),
  simDataKumar2 = c(2:10)
)
n.genes <- list(
  kumar2015 = 1000,
  trapnell2014 = 1000,
  zhengmix2016 = 200,
  koh2016 = 1000,
  simDataKumar = 1000,
  simDataKumar2 = 1000
  
  
)


# zinbwave function

run_zinbwave <- function(data, par.k, n.genes) {

list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->transformedExp-> res.zinb ->d -> tsne_data -> res.cluster
for (i in names(data)) {
  
  filter <- rowSums((counts(data[[i]]))>5)>5 
  data[[i]] <- data[[i]][filter,]# filter out genes with at least five counts
  counts( data[[i]] ) %>% log1p %>% rowVars -> vars
  names(vars) <- rownames( data[[i]] )
  vars <- sort(vars, decreasing = TRUE)
  data[[i]] <- data[[i]][names(vars)[1:n.genes[[i]]],]
  
  res.zinb[[i]] <- zinbFit( round(counts(data[[i]]),0) ,
                            K=2, epsilon=n.genes[[i]], verbose=TRUE,
                            nb.repeat.initialize = 2, 
                            maxiter.optimize = 25,
                            stop.epsilon.optimize = 1e-04)  # round data as it assumes whole counts

  
  d[[i]]<- dist(getW( res.zinb[[i]] ))
  
  df.clus <- matrix( nrow = ncol( data[[i]]), ncol = length(par.k[[i]]) )
  
  for ( j in seq_len(length(par.k[[i]])) ) {
    df.clus[,j]<- as.integer( kmeans( d[[i]], centers = par.k[[i]][j])$cluster )
    }
  colnames(df.clus) <-  c( paste0(par.k[[i]]) )
  res.cluster[[i]] <- df.clus
}
return(res.cluster)

}

#
res.cluster <- run_zinbwave(data, par.k, n.genes)

# save clusters

dir_cluster <- paste0("results/filtered/zinbwave/zinbwave_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)



# save experiment labels

file_names <-  paste0("results/filtered/zinbwave/zinbwave_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/filtered/zinbwave/session_info_zinbwave_krange.txt")
sessionInfo()
sink()

### Appendix


