#####################
# SNN-Cliq
#####################


# file paths

DATA_DIR <- "~/Desktop/masterthesis/data"
files <- list(
  kumar2015 = file.path(DATA_DIR, "sceset_GSE60749-GPL13112.rds"),
  trapnell2014 = file.path(DATA_DIR, "sceset_GSE52529-GPL16791.rds"),
  xue2013 = file.path(DATA_DIR, "sceset_GSE44183-GPL11154.rds")
)

# load data sets

data <- vector("list", length(files))

names(data) <- names(files)

for (i in 1:length(data)){
  f <- files[[i]]
  load(f)
  data[[i]] <- res
  
}
# create storage files
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- input_matrix<- labels<- list

# load cell labels
for(i in 1:length(data)) {
  labels[[i]] <- as.character(phenoData(data[[i]])@data$source_name_ch1)
}


# extract expression data
for (i in 1:(length(input_matrix))){
  input_matrix[[i]] <-exprs(data[[i]]) # use exprs slot of SCeset; log2, normalized count_lstpm
}


# RUN SNN-Cliq


# Set parameters, 

distan <- "euclidean"

par.k <-  list(
  kumar2015 <- 3,
  trapnell2014 <- 3,
  xue2013 <- 8
)

par.r <- list(
  kumar2015 <- 0.7,
  trapnell2014 <- 0.7,
  xue2013 <- 0.7
)
par.m <-  list(
  kumar2015 <- 0.5,
  trapnell2014 <- 0.5,
  xue2013 <- 0.5
)


for (i in 1:length(input_matrix)){
  # construct a graph 
scRNA.seq.funcs::SNN(
  data = t(input_matrix[[i]]),
  outfile = "snn-cliq.txt",
  k = par.k[[i]],
  distance = distan
)
  

# find clusters in the graph
setwd("~/Desktop/masterthesis/skript")
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
#
cat(paste(snn.res, collapse = "\n"))
snn.res <- read.table("res-snn-cliq.txt")

pData(data[[i]])$SNNCliq <- as.character(snn.res[,1])
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")
}





# save clusters

file_names <- paste0("~/Desktop/masterthesis/results/SNNCliq/SNNCliq_clus_", 
                     names(res.cluster), ".txt")

for (i in 1:length(file_names)) {
  res_i <- pData(data[[i]])$SNNCliq
  write.table(res_i, file = file_names[i], row.names = FALSE, quote = FALSE, sep = "\t")
}
# save systemtime
file_names <-  paste0("~/Desktop/masterthesis/results/SNNCliq/SNNCLiq_systime_",names(sys.time))
for (i in 1:length(sys.time)){
  sys_i <- sys.time[[i]]["elapsed"]
  write.table(sys_i, file=file_names[i], sep="\t")
  
}

# save experiment labels

file_names <-  paste0("~/Desktop/masterthesis/results/SNNCliq/SNNCliq_labels_",names(sys.time), ".txt")
for (i in 1:length(sys.time)){
  sys_i <- as.data.frame(labels[[i]])
  write.table(sys_i, file=file_names[i], sep="\t")
}

###### Save Session Info
sink(file = "~/Desktop/masterthesis/results/SNNCliq/session_info_SNNCliqreduce.txt")
sessionInfo()
sink()

# Appendix


