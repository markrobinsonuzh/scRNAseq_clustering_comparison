comma := ,
empty :=
space := $(empty) $(empty)

FILTERINGS := Expr
PCTS := 50
ALLFILTERINGS := full $(foreach F,$(FILTERINGS),$(foreach P,$(PCTS),$(F)$(P)))
ALLFILTERINGSc := $(subst $(space),$(comma),$(ALLFILTERINGS))
