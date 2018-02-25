comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := full filteredExpr
FILTERINGSc := $(subst $(space),$(comma),$(FILTERINGS))