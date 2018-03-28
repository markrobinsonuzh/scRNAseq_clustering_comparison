comma := ,
empty :=
space := $(empty) $(empty)

DATASETSbig := Zhengmix4eq Zhengmix4uneq Zhengmix8eq
DATASETSsmall := Kumar Trapnell Koh SimKumar4easy SimKumar4hard SimKumar8hard KohTCC TrapnellTCC KumarTCC
DATASETS := $(DATASETSbig) $(DATASETSsmall)

DATASETSbigc := $(subst $(space),$(comma),$(DATASETSbig))
DATASETSsmallc := $(subst $(space),$(comma),$(DATASETSsmall))
DATASETSc := $(subst $(space),$(comma),$(DATASETS))

