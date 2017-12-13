#######################################################
# Stability analysis of methods using bootstrap
# for one dataset only
#######################################################

# load librareis
library(dplyr)
library(mclust)
# load helper files
source("skript/helper_files/Helper_functions.R")
#load methods for stability anlysis
source("skript/stability_analysis/interface_clusterboot.R")


#-------------------------------------------------------------

# load a data set and the labels thereof
source("FILES.R")
data <- load_data(files, DATA_DIR)
labels <- load_labels(data)[[1]] %>% as.factor %>% as.integer
data.normcts <- assay(data[[1]], "normcounts") %>% rbind(labels=labels)
data.raw <-  assay(data[[1]], "counts") %>% rbind(labels=labels)

#-------------------------------------------------------------

# sample by replacement 20 times
# function

sample_rep <- function(data) {
  data[, sample(1:ncol(data), ncol(data), replace = TRUE)]
}
#bootstrap
sample.norm <- NULL
for ( i in 1:20) {
sample.norm[[i]] <-  (sample_rep(data.normcts ) )
}
sample.raw <- NULL
for ( i in 1:20) {
  sample.raw[[i]] <-  (sample_rep(data.raw ) )
}
# extract labels
labels.norm <- lapply(sample.norm, "[",  "labels", 1:ncol(data.normcts)) %>%lapply(as.integer)

labels.raw <- lapply(sample.raw, "[",  "labels", 1:ncol(data.raw))%>%lapply(as.integer)
# extract sample
sample.norm <- lapply(sample.norm, function(x) x[ -grep("labels",rownames(x)), ] )
sample.raw <- lapply(sample.raw, function(x) x[ -grep("labels",rownames(x)), ] )
#-------------------------------------------------------------
# run methods over bootstrap samples

#boot.rtsnekmeans <- lapply(  sample.norm, rtsnekmeansCBI, k=3, perplexity=30 ) %>% lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.norm)

boot.cidr <- lapply( sample.norm, cidrCBI, par.k=3, par.nPC=4) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.norm) 

boot.linnorm <- lapply( sample.raw, linnormCBI, par.minNonZeroPortion=0.75, par.num_center=3 ,par.BE_strength=0.5) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.raw)

boot.simlr <- lapply( sample.norm, simlrCBI,par.c=3, normalize =TRUE ) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.norm)

#boot.sc3 <- lapply( sample.norm, sc3CBI , par.ks = NULL, par.k_estimator=FALSE , par.k =3, pct_dropout_max =90) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.norm)

boot.pcareduce <- lapply( sample.norm, pcareduceCBI, par.nbt=100, par.q=30 ,n.cluster=3) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.norm)

boot.seurat <- lapply( sample.raw,seuratCBI,par.resolution=0.6, k.param = 25 , par.dims.use=1:9 ) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.raw)

boot.tscan <- lapply( sample.raw,tbscanCBI,par.minexpr_percent=0.5  ,par.clusternum=3 ) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.raw)

boot.raceid <- lapply( sample.raw,raceidCBI,par.mintotal=3000, par.maxexpr=Inf ,do.gap=FALSE,cln=3 ) %>%lapply("[[", 1) %>% mapply(FUN=adjustedRandIndex, labels.raw)

#-------------------------------------------------------------
# save ari results
save(boot.cidr , file = "results/stability_analysis/boot.cidr.rda")
save(boot.linnorm , file = "results/stability_analysis/boot.linnorm.rda")
save(boot.simlr , file = "results/stability_analysis/boot.simlr.rda")
save(boot.pcareduce  , file = "results/stability_analysis/boot.pcareduce.rda")
save(boot.seurat  , file = "results/stability_analysis/boot.seurat.rda")
save(boot.tscan  , file = "results/stability_analysis/boot.tscan.rda")
save(boot.raceid , file = "results/stability_analysis/boot.raceid.rda")


write.table(boot.cidr , file = "results/stability_analysis/boot.cidr.txt", sep="\t")
write.table(boot.linnorm , file = "results/stability_analysis/boot.linnorm.txt", sep="\t")
write.table(boot.simlr , file = "results/stability_analysis/boot.simlr.txt", sep="\t")
write.table(boot.pcareduce  , file = "results/stability_analysis/boot.pcareduce.txt", sep="\t")
write.table(boot.seurat  , file = "results/stability_analysis/boot.seurat.txt", sep="\t")
write.table(boot.tscan  , file = "results/stability_analysis/boot.tscan.txt", sep="\t")
write.table(boot.raceid , file = "results/stability_analysis/boot.raceid.txt", sep="\t")





res.boot <- list(   rtsnekmeans = boot.rtsnekmeans, 
                    cidr = boot.cidr, 
                    linnorm = boot.linnorm , 
                    simlr = boot.simlr , 
                    sc3 =boot.sc3, 
                    pcareduce = boot.pcareduce, 
                    seurat =boot.seurat,
                    tscan = boot.tscan,
                    raceid=boot.raceid )


save(res.boot , file = file.path())

