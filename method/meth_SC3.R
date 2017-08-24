####################
# SC3
######################

# load libraries

library("scater")
library("DESeq2")
library("SC3")

# file paths


DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rds")
)

# load data sets

data <- vector("list", length(files))
labels <- data
names(data) <- names(labels) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in 1:length(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$source_name_ch1)
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
  xue2013=8
  
)


for (i in 1:length(data)){
  sys.time[[i]] <- system.time({
  data[[i]]<- sc3_prepare(data[[i]], ks=par.ks[i])        # uses the exprs slot of SCEset ; log2transformed, normalized data
  data[[i]]<- sc3_estimate_k(data[[i]])# optional
  data[[i]]<- sc3(data[[i]], ks = data[[i]]@sc3$k_estimation, biology = FALSE)
  })
  # store clusters
  p_data <- pData(data[[i]])
  res.cluster[[i]] <- p_data[ , grep("sc3_", colnames(p_data))]
  
}

# save clusters

file_names_clus <- paste0("~/Desktop/masterthesis/results/SC3/sc3_clus_", 
                     names(res.cluster), ".txt")

for (i in 1:length(file_names_clus)) {
  res_i <- res.cluster[[i]]
  write.table(res_i, file = file_names_clus[i], row.names = FALSE, quote = FALSE, sep = "\t")
}
# save systemtime
file_names_systime <-  paste0("~/Desktop/masterthesis/results/SC3/sc3_systime_",names(sys.time))
for (i in 1:length(sys.time)){
  sys_i <- sys.time[[i]]["elapsed"]
  write.table(sys_i, file=file_names_systime[i], sep="\t")
  
}
# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/SC3/sc3_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}
###### Save Session Info
sink(file = "~/Desktop/masterthesis/results/SC3/session_info_sc3.txt")
sessionInfo()
sink()

### Appendix
#


