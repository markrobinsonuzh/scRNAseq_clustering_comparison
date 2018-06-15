comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := Expr M3Drop HVG
PCTKEEP := 10
ALLFILTERINGS := $(foreach F,$(FILTERINGS),$(foreach P,$(PCTKEEP),filtered$(F)$(P)))# full
ALLFILTERINGSc := $(subst $(space),$(comma),$(ALLFILTERINGS))
