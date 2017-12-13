
#####################
# RACEID
#####################
# RaceID is an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. 
# The method is based on transcript counts obtained with unique molecular identifies.


run_function_raceid <- function( data, labels, datatype ,par.mintotal, par.maxexpr,do.gap,cln ) {
  
  
  source("method_resources/RaceID/RaceID_class.R")
  source("skript/helper_files/Helper_functions.R")
  require(tsne)
  require(pheatmap)
  require(MASS)
  require(cluster)
  require(mclust)
  require(flexmix)
  require(lattice)
  require(fpc)
  require(amap)
  require(RColorBrewer) 
  require(locfit)

# create store vectors
list<- vector("list", length(data))
names(list) <- names(data)
list->sys.time->sc->res.cluster 
for (i in names(data)){
  print(i)
  sys.time [[i]] <- system.time({
    sc[[i]] <- SCseq(as.data.frame(counts(data[[i]]))) # extract the expression data
    sc[[i]] <- filterdata(sc[[i]], mintotal = par.mintotal[[i]], minexpr = 5, 
                     minnumber = 1, maxexpr = par.maxexpr[[i]] , downsample = FALSE, dsn = 1, rseed = 1234)
    
    res.cluster[[i]]<- clustexp(sc[[i]], metric="pearson", cln=cln[[i]], do.gap=do.gap[[i]], clustnr=20, B.gap=50, SE.method="Tibs2001SEmax", 
                   SE.factor=.25, bootnr=50, rseed=1234)@kmeans$kpart
  })
  
  
}
# save clusters

dir_cluster <- paste0("results/",datatype,"/raceid/raceid_clus_", names(res.cluster), ".txt")


save_clusters(res.cluster,dir_cluster)

# save systemtime

dir_systime <-  paste0("results/",datatype,"/raceid/raceid_systime_",names(sys.time),".txt")

save_systemtime(sys.time, dir_systime)

# save experiment labels

dir_labels <-  paste0("results/",datatype,"/raceid/raceid_labels_",names(labels), ".txt")
save_labels(labels, dir_labels )


###### Save Session Info
sink(file = "results/",datatype,"/raceid/session_info_raceid.txt")
sessionInfo()
sink()

}
