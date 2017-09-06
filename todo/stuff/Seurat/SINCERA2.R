library(SINCERA)
library(SC3)
library(pheatmap)


# import data

load(file = "~/Desktop/masterthesis/data/sceset_final.Rda")

# create SCEset
cts <- assays(experiments(res)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(res)[["gene"]])[["TPM"]]
phn <- pData(res)
res<- newSCESet(
  countData = cts, 
  tpmData = tpms,
  phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn))
)

# pcaReduce : use genefilter
# use the same gene filter as in SC3


res <- calculateQCMetrics(res)
res <- sc3_prepare(res, ks = 2:5)
res<- sc3_estimate_k(res)
res@sc3$k_estimation
res <- sc3(res, ks = 2, biology = TRUE)

# extract expression data
input <- exprs(res[fData(res)$sc3_gene_filter, ])

# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
#dissimilarities
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)

hc <- hclust(dd, method = "average")
# !whats happenng here, slice clusters so that no singleton exists
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)

# store variable
pData(res)$SINCERA <- as.character(cutree(hc, k = i))
plotPCA(res, colour_by = "SINCERA")
#
# Let's now visualize the SINCERA results as a heatmap:
pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = 14,
  kmeans_k = 100,
  show_rownames = FALSE
)


