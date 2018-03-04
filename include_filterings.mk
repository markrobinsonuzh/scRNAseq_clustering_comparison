comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := Expr HVG
PCTKEEP := 10 50
ALLFILTERINGS := $(foreach F,$(FILTERINGS),$(foreach P,$(PCTKEEP),filtered$(F)$(P)))# full
ALLFILTERINGSc := $(subst $(space),$(comma),$(ALLFILTERINGS))
