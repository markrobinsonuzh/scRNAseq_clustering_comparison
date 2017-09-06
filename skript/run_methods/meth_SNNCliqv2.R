#####################
# SNN-Cliq
#####################
source("skript/helper_files/Helper_functions.R")


# file paths

DATA_DIR <- "data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rda"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rda"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rda")
)

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


# Set parameters, k is number of nearest neighbours, r is cutoff for quasi-clique, m is threshold of overlapping quasicliques to be merged (standart threshold is half of the points)

distan <- "euclidean"

par.k <-  list(
  kumar2015 = 15,
  trapnell2014 = 20,
  xue2013 = 3
)

par.r <- list(
  kumar2015 = 0.7,
  trapnell2014 = 0.7,
  xue2013 = 0.7
)
par.m <-  list(
  kumar2015 = 0.5,
  trapnell2014 = 0.5,
  xue2013 = 0.5
)

write.csv(x=NULL,file="skript/run_methods/snn-cliq.txt")
for (i in names(input_matrix)){
  # construct a graph 
scRNA.seq.funcs::SNN(
  data = t(input_matrix[[i]]),
  outfile = "snn-cliq.txt",
  k = par.k[[i]],
  distance = distan
)
  

# find clusters in the graph
setwd("skript/run_methods")
sys.time[[i]] <- system.time({
snn.res <- 
  system(
    paste0(
      "python snn-cliq/Cliq.py ", 
      "-i snn-cliq.txt ",
      "-o res-snn-cliq.txt ",
      "-r ", par.r[[i]],
      " -m ", par.m[[i]]
    ),
    intern = TRUE
  )
})
#
cat(paste(snn.res, collapse = "\n"))
snn.res <- read.table("res-snn-cliq.txt")

pData(data[[i]])$SNNCliq <- as.character(snn.res[,1])
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")
}


# save clusters
for (i in seq_len(length(res.cluster))){
  res.cluster[[i]] <- pData(data[[i]])$SNNCliq
}
dir_cluster <- paste0("results/SNNCliq/SNNCliq_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/SNNCliq/SNNCliq_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)


# save experiment labels

file_names <-  paste0("results/SNNCliq/SNNCliq_labels_",names(labels), ".txt")
for (i in 1:length(labels)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "results/SNNCliq/session_info_SNNCliqreduce.txt")
sessionInfo()
sink()

# Appendix


