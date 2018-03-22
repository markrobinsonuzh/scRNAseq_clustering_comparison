comma := ,
empty :=
space := $(empty) $(empty)

METHODSbig := PCAKmeans RtsneKmeans Seurat FlowSOM SIMLRlargescale SC3 CIDR RaceID PCAHC
METHODSsmall := PCAKmeans RtsneKmeans Seurat FlowSOM pcaReduce SIMLRlargescale SC3 SIMLR CIDR Linnorm RaceID SC3svm TSCAN PCAHC
METHODS := $(METHODSsmall)

METHODSbigc := $(subst $(space),$(comma),$(METHODSbig))
METHODSsmallc := $(subst $(space),$(comma),$(METHODSsmall))
METHODSc := $(subst $(space),$(comma),$(METHODS))