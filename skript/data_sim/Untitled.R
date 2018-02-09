################################################
#ScriptforsimulateddatausingSCsim#
################################################

#sourcethepackage
source("method_resources/SCsim-master/RAllGenerics.R")
source("method_resources/SCsim-master/RBiais_methods.R")
source("method_resources/SCsim-master/RClasses.R")
source("method_resources/SCsim-master/RDistribution_methods.R")
source("method_resources/SCsim-master/RMisc_methods.R")
source("method_resources/SCsim-master/RSCsimSet_methods.R")
#Parametersetting
nGenes=5000
nCells=600
nPop=6
pPop=c(30,20,5,10,30,5)
seed=125

#Definebasalgeneexpression(default)
distribution="gamma"
fstParam=2
sndParam=0.75

#Definebatcheffectinthedata
nbBatch=1
cellsPerBatch=NULL
batchEffect=NULL

#Definedifferentiallyexpressedgenesproportion
pDEG=c(4,4,3,2,2,5)
pDE=c(20,40,30,40,10,60)
pDP=20
pDM=c(40,20,30,20,50,0)
pDC=20
pUp=c(70,50,40,60,70,30)
pDown=c(30,50,60,40,30,70)

cellMixedDP="pseudo"
mixDP=25
cellMixedDM="pseudo"
mixDM=25
popMixDP=NULL
trajectory=list(c(1,2,3,4),c(1,2,5))
doublet=2

distrUpFc="medium"
distrDownFc="medium"

dropoutPct=50


t<-newSCsimSet(nGenes=nGenes,nCells=nCells,nPop=nPop,pPop=pPop,seed=seed,distribution=distribution,fstParam=fstParam,sndParam=sndParam,
nbBatch=nbBatch,cellsPerBatch=cellsPerBatch,batchEffect=batchEffect,pDEG=pDEG,pDE=pDE,pDP=pDP,pDM=pDM,pDC=pDC,pUp=pUp,pDown=pDown,
cellMixedDP=cellMixedDP,mixDP=mixDP,cellMixedDM=cellMixedDM,mixDM=mixDM,popMixDP=popMixDP,trajectory=trajectory,doublet=doublet,distrUpFc=distrUpFc,distrDownFc=distrDownFc,dropoutPct=dropoutPct
)
# extract expression matrix
expr <- t@cell_table
# extract the cell labels
phenoid <- t@batch_table$population
# create Single Cell Expression Set object

# save 






