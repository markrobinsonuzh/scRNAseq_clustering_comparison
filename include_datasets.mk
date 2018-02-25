comma := ,
empty :=
space := $(empty) $(empty)

DATASETS := Kumar Trapnell Koh Zhengmix4eq Zhengmix4uneq Zhengmix8eq SimKumar4easy SimKumar4hard SimKumar8hard
DATASETSc := $(subst $(space),$(comma),$(DATASETS))

