comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := filteredExpr# full
FILTERINGSc := $(subst $(space),$(comma),$(FILTERINGS))