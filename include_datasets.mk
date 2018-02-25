comma := ,
empty :=
space := $(empty) $(empty)

DATASETS := Kumar Trapnell Koh# Zhengmix# SimKumar
DATASETSc := $(subst $(space),$(comma),$(DATASETS))

