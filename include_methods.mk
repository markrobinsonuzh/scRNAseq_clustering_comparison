comma := ,
empty :=
space := $(empty) $(empty)

METHODS := PCAKmeans RtsneKmeans Seurat FlowSOM pcaReduce SIMLRlargescale SC3 SIMLR CIDR Linnorm RaceID SC3svm TSCAN
METHODSc := $(subst $(space),$(comma),$(METHODS))