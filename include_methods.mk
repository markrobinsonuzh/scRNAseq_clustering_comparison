comma := ,
empty :=
space := $(empty) $(empty)

METHODS := PCAKmeans RtsneKmeans Seurat FlowSOM SC3 CIDR PCAHC SC3svm pcaReduce TSCAN ascend SAFE monocle RaceID2# SIMLR

METHODSc := $(subst $(space),$(comma),$(METHODS))
