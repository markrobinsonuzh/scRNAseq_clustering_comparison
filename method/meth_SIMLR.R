#####################
# SIMLR
#####################


#load libraries

library(SIMLR)
library(igraph)
library(scater)

set.seed(1234)

# file paths

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rds")
)

# load data sets

data <- labels<- vector("list", length(files))

names(labels) <- names(data) <- names(files)



for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}

# load cell labels
for(i in 1:length(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$source_name_ch1)
}


# RUN SIMLR
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- list
# Set paramaeters
par.c <-  list(
  kumar2015 <- 3,
  trapnell2014 <- 3,
  xue2013 <- 8
)

# 
for (i in 1:length(data)){
  
  res.cluster[[i]] = SIMLR(X = exprs(data[[i]]), c = par.c[[i]], cores.ratio = 0) # use exprs slot of SCeset; log2, normalized count_lstpm

}

# save clusters

file_names <- paste0("~/Desktop/masterthesis/results/SIMLR/SIMLR_clus_", 
                     names(res.cluster), ".txt")

for (i in 1:length(file_names)) {
  res_i <- res.cluster[[i]]$y$cluster    
  write.table(res_i, file = file_names[i], row.names = FALSE, quote = FALSE, sep = "\t")
}
# save systemtime
file_names <-  paste0("~/Desktop/masterthesis/results/SIMLR/SIMLR_systime_",names(sys.time))
for (i in 1:length(sys.time)){
  sys_i <- sys.time[[i]]["elapsed"]
  write.table(sys_i, file=file_names[i], sep="\t")
  
}

# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/SIMLR/SIMLR_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "~/Desktop/masterthesis/results/SIMLR/session_info_SIMLR.txt")
sessionInfo()
sink()

# Appendix


