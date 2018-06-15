## Apply FlowSOM

suppressPackageStartupMessages({
  library(flowCore)
  library(FlowSOM)
})

apply_FlowSOM <- function(sce, params, k) {
  tryCatch({
    dat <- logcounts(sce)
    st <- system.time({
      pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
      pca <- pca$x[, seq_len(params$nPC), drop = FALSE]
      ff <- flowFrame(exprs = pca)
      fSOM <- FlowSOM::ReadInput(ff, compensate = FALSE, transform = FALSE, 
                                 scale = FALSE, silent = TRUE)
      fSOM <- FlowSOM::BuildSOM(fSOM, silent = TRUE, xdim = params$xdim, 
                                ydim = params$ydim)
      metaClustering <- metaClustering_consensus(fSOM$map$codes, k = k)
      cluster <- metaClustering[fSOM$map$mapping[, 1]]
      names(cluster) <- colnames(dat)
    })
    
    st <- c(user.self = st[["user.self"]], sys.self = st[["sys.self"]], 
            user.child = st[["user.child"]], sys.child = st[["sys.child"]],
            elapsed = st[["elapsed"]])
    list(st = st, cluster = cluster, est_k = NA)
  },
  error = function(e) {
    list(st = c(user.self = NA, sys.self = NA, user.child = NA, sys.child = NA,
                elapsed = NA),
         cluster = structure(rep(NA, ncol(sce)), names = colnames(sce)),
         est_k = NA)
  })
}
