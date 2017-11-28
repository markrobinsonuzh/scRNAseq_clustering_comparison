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
  kumar2015 = c(2:12),
  trapnell2014 =  c(2:12),
  zhengmix2016 =  c(2:12),
  koh2016 =  c(2:12)
)


# zinbwave function

run_zinbwave <- function(data, par.k) {

list<- vector("list", length(files))
names(list) <- names(files)
list->sys.time->transformedExp-> res.zinb ->d -> tsne_data -> res.cluster
for (i in names(data)) {
  
  filter <- rowSums((counts(data[[i]]))>5)>5 
  data[[i]] <- data[[i]][filter,]# filter out genes with at least five counts
  counts( data[[i]] ) %>% log1p %>% rowVars -> vars
  names(vars) <- rownames( data[[i]] )
  vars <- sort(vars, decreasing = TRUE)
  data[[i]] <- data[[i]][names(vars)[1:1000],]
  
  res.zinb[[i]] <- zinbFit( round(counts(data[[i]]),0) , K=2, epsilon=1000, verbose=TRUE)  # round data as it assumes whole counts
  d[[i]]<- dist(getW( res.zinb[[i]] ))
  tsne_data[[i]] <- Rtsne(d[[i]], is_distance = TRUE, pca = FALSE, 
                          perplexity=10, max_iter=5000)
  res.cluster[[i]] <- kmeans(d[[i]], centers=par.k[[i]] )$cluster
  
}

}

#
res.cluster <- run_zinbwave(data, par.k)

# save clusters

dir_cluster <- paste0("results/zinbwave/zinbwave_krange_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/zinbwave/zinbwave_krange_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/zinbwave/zinbwave_krange_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}


###### Save Session Info
sink(file = "results/zinbwave/session_info_zinbwave_krange.txt")
sessionInfo()
sink()

### Appendix


