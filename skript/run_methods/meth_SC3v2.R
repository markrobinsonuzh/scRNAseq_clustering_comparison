####################
# SC3
######################
# change sc 3 prepare
source("skript/helper_files/Helper_functions.R")

# load libraries

library("scater")
library("DESeq2")
library("SC3")

# file paths


DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda"),
  koh2016 = file.path(DATA_DIR,"sceset_SRP073808.rda")
)

# load data sets

data <- vector("list", length(files))
labels <- data
names(data) <- names(labels) <- names(files)

for (i in names(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in names(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$phenoid)
}


# RUN SC3
# list to store results
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- list




# QC data
# done in scater

# run the analysis
par.ks <- list(
  kumar2015=3,
  trapnell2015=3,
  xue2013=8,
  koh2016 = 10
  
)


for (i in names(data)){
  sys.time[[i]] <- system.time({
    data[[i]]<- sc3_prepare(data[[i]], ks=par.ks[i])        # uses the exprs slot of SCEset ; log2transformed, normalized data
    data[[i]]<- sc3_estimate_k(data[[i]])# optional
    data[[i]]<- sc3(data[[i]], ks = data[[i]]@sc3$k_estimation, biology = FALSE)
  })
  # store clusters
  p_data <- pData(data[[i]])
  res.cluster[[i]] <- p_data[ , grep("sc3_", colnames(p_data))]
  
}




# save clusters

dir_cluster <- paste0("results/SC3/sc3_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SC3/sc3_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

file_names <-  paste0("results/SC3/sc3_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}
###### Save Session Info
sink(file = "results/SC3/session_info_sc3.txt")
sessionInfo()
sink()

### Appendix
#

