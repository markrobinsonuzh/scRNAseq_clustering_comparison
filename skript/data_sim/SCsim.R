################################################
#ScriptforsimulateddatausingSCsim#
################################################
# libraries
library(dplyr)
library(scater)
library(cowplot)
#sourcethepackage
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/AllGenerics.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/Biais_methods.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/Classes.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/DEG_methods.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/Distribution_methods.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/Misc_methods.R')
source('~/Desktop/masterarbeit/scRNAseq_clustering_comparison/method_resources/SCsim-master/R/SCsimSet_methods.R')
#Parametersetting
nGenes=1000
nCells=300
nPop=5# 
pPop=c(30,20,5,10,35) #
seed=125

#Definebasalgeneexpression(default) , mean gene expression distribution 
distribution="gamma" # which distribution, gamma or NB
fstParam=2 # shpae or size pararmeter in gamma or NB
sndParam=0.75# rate or mean in gamma or NB

#Definebatcheffectinthedata
nbBatch=2# nr of batches
cellsPerBatch=c(30, 70) # cells per batch,  proportions
batchEffect=c(0.1,0.2) # mean shift in library size for each batch.

# Define differentially expressed genes proportion, Differentially expressed genes indexes / type / regulation in
#' the gene expression table to generate
pDEG=c(4,4,3,2,10)# fraction of  DE genes per population
pDE=c(20,40,30,40,10)# proportion of fixed up and down regulated genes seperated per poopulation
pDP=30 # proportion of genes on trajectory between supopulations
pDM=c(30,10,20,10,40)# proportion of genes , DEexpressed with propability 1-x?
pDC=20 # whats that?

pUp=c(70,50,40,60,70) # proportion up
pDown=c(30,50,60,40,30)# and down regulated

cellMixedDP="pseudo"
mixDP=10 # fraction of genes on trajectory that show DE
cellMixedDM="pseudo"
mixDM=25
popMixDP=NULL
trajectory=list(c(1,2,3),c(3,4,5))# which trajectory between cells
doublet=2 # how many?

distrUpFc="medium"
distrDownFc="medium"

dropoutPct=40 #


t<-newSCsimSet(nGenes=nGenes,nCells=nCells,nPop=nPop,pPop=pPop,seed=seed,distribution=distribution,fstParam=fstParam,sndParam=sndParam,
nbBatch=nbBatch,cellsPerBatch=cellsPerBatch,batchEffect=batchEffect,pDEG=pDEG,pDE=pDE,pDP=pDP,pDM=pDM,pDC=pDC,pUp=pUp,pDown=pDown,
cellMixedDP=cellMixedDP,mixDP=mixDP,cellMixedDM=cellMixedDM,mixDM=mixDM,popMixDP=popMixDP,trajectory=trajectory,doublet=doublet,distrUpFc=distrUpFc,distrDownFc=distrDownFc,dropoutPct=dropoutPct
)
# extract expression matrix
cts <- t@effectiveCounts%>%as.matrix
# extract the cell labels
phenoid <- t@batch_table$population

# create Single Cell Expression Set object

sceset <- SingleCellExperiment(assays=list(counts=cts))
colData(sceset)$phenoid <- phenoid
# plot 
p.meanexpr <- plotBaseMean(t) # mean expression of genes
p.batch <- plotBatch(t) # batch effect, distributions of batch effects per cells
p.libeffect <- plotCellBiais(t) # library effect , distr. of lib effects per cell population
p.pcell <- plotDatasetInfo(t) # fractions of subpopulations
p.cont <- plotCellsContent(t)  # 
p.deg <- plotDEG(t) # fractions of DEG per populations
p8 <- plotDEGDensity(t) # nort working
p.typedeg <- plotDEGTypes(t) # type of DEG, DC=?, DE=?, DM=?, DP=?
p.dropp <- plotDropoutProba(t)# dropout probabilites per cell and per gene, before after ?
p.drpdist <- plotDropoutQC(t) # dropout distr. gene and cell wise
p.libsize <- plotLibrarySize(t) # dist. of libsize per population
p.tsne <-  plotTSNE(sceset, exprs_values = "counts",colour_by = "phenoid" )
p.pca <-  plotPCA(sceset, exprs_values = "counts", colour_by = "phenoid" )
p.grid <- plot_grid( p.meanexpr,p.pcell, p.typedeg , p.deg ,p.cont$gtable, p.batch , p.libeffect, p.libsize, p.pca  , ncol = 3 )
save_plot(filename="results/QC_data/scsim.pdf",p.grid, base_height=10)
# save 
save(sceset, file = "data/scsim.rda" )



