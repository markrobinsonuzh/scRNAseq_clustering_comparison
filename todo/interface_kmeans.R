############################################################
# tSNE kmeans Interface function for bootclust
#########################################################

# number of clusters in kmeans
rand.seed <- 1234
par.k <- list(
  kumar2015 = 3,
  trapnell2014 = 12,
  xue2013 = 8
)

par.perp <- list(
  kumar2015 = 20,
  trapnell2014 = 20,
  xue2013 = 5
)
# Run tSNE and kmeans
list <- vector("list", length(data))
names(list) <- names(data)
res.cluster <- sys.time<- list


for (i in names(data)){
  
  sys.time [[i]] <- system.time({
    data[[i]] <- plotTSNE(data[[i]], exprs_values= "counts",rand_seed = rand.seed, perplexity= par.perp[[i]],return_SCESet = TRUE, draw_plot= FALSE) # use Rtsne? function plotTSNE uses Rtsne anyway
    pData(data[[i]])$tSNE_kmeans <- as.character(kmeans(data[[i]]@reducedDimension, centers = par.k[[i]])$clust)
  })
  res.cluster[[i]] <- pData(data[[i]])$tSNE_kmean
  
}


function() {
  
  
  
}

iris %>%
  subset(Sepal.Length > mean(Sepal.Length)) %>%
  cor(Sepal.Length, Sepal.Width)
