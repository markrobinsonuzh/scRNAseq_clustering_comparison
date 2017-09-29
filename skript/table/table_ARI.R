###################
# Table ARI
###################
library(stargazer)
library(dplyr)
library(plyr)

# define data directory, METHODS and DATASETS

RES_DIR <- "results/"

DATASET <- c(
  "kumar2015",
  "trapnell2014", 
  "xue2013",
  "koh2016"
)
METHOD <- c(
  "tSNEkmeans",
  "SIMLR",
  "SC3",
  "SNNCliq",
  "pcaReduce",
  "seurat",
  "dbscan",
  "cidr"
  
)

# create file directories
file_names <- vector("list", length = length(METHOD))
names(file_names) <- METHOD

for (i in seq_len(length(file_names))){
  
  file_names[[i]] <- file.path(paste0(RES_DIR, METHOD[i],"/",METHOD[i],"_ARI_", DATASET[seq_len(length(DATASET))],".txt" ))
  names(file_names[[i]]) <- DATASET
}


# load data
ARI <- sapply(unlist(file_names),read.csv)

#### Table
ari.m <- matrix(ARI,nrow = length(DATASET), dimnames = list(DATASET,METHOD)) 
ari.m[] <- lapply(ari.m, round,2)



tbl.ari <- stargazer(ari.m, title = "Adjusted rand index ")

### new table

ari.m2 <- ari.m[1:4,c(-4) ]
rownames(ari.m2) <- c("Kumar 2015", "Trapnell 2014", "Xue 2013", "Koh 2016")
colnames(ari.m2) <- c("tSNE kmeans", "SIMLR", "SC3" , "pcaReduce" ,"SEURAT" ,"dbscan", "CIDR")

 stargazer(ari.m2, title = "Results of comparison of clustering methods for data set Kumar 2015, Trapnell 2014, Xue 2013 and Koh 2016. Shown is the Adjusted Rand Index for the methods tSNE and kmeans, SIMLR, SC3, pcaReduce, SEURAT and dbscan")


### Appendix
ari.df <- ldply(ARI, data.frame)
adply(ari.df,.margins = 1,c(".id","X..i.."), grep("kumar",ari.df$.id))
sub("\\.", "", ari.df$.id)
strip_splits(ari.df$.id,"\\.")

dply(strsplit(ari.df$.id,"\\."),data.frame)

str_(ari.df$.id,"\\.")
ggplot(data=ari.df, aes(x=ari.df$.id, y=ari.df$X..i..))+geom_bar(stat = "identity")+theme_grey()

