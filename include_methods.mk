comma := ,
empty :=
space := $(empty) $(empty)

METHODS := PCAKmeans RtsneKmeans Seurat FlowSOM SC3 CIDR PCAHC SC3svm pcaReduce TSCAN ascend# RaceID2

METHODSc := $(subst $(space),$(comma),$(METHODS))