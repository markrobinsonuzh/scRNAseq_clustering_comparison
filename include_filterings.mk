comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := filtered# full
FILTERINGSc := $(subst $(space),$(comma),$(FILTERINGS))