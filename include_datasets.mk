comma := ,
empty :=
space := $(empty) $(empty)

DATASETS := Kumar Trapnell Koh Zhengmix4eq Zhengmix4uneq Zhengmix8eq SimKumarEasy SimKumarHard
DATASETSc := $(subst $(space),$(comma),$(DATASETS))

