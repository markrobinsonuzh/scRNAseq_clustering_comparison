## Define the versions of R and the paths to the libraries
#R := R_LIBS=/home/Shared/Rlib/release-3.6-lib/ /usr/local/R/R-3.4.2/bin/R CMD BATCH --no-restore --no-save
#Rscript := /usr/local/R/R-3.4.2/src/unix/Rscript --no-restore --no-save
R := R CMD BATCH --no-restore --no-save
Rscript := Rscript --no-restore --no-save

DATASETS := GSE52529-GPL16791 GSE60749-GPL13112 SRP073808 Zhengmix SimKumar

.PHONY: all

all: 

## Set up directory structure and download raw data
setup: 
	mkdir -p data/data_raw
	mkdir -p data/sce_full
	mkdir -p data/sce_filtered
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE60749-GPL13112.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE52529-GPL16791.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/SRP073808.rds
	mkdir -p plots/qc_data

## Prepare data sets and generate QC plots
data/sce_filtered/sce_filtered_Trapnell.rds: Rscripts/import_datasets/import_QC_Trapnell.Rmd \
data/data_raw/GSE52529-GPL16791.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

data/sce_filtered/sce_filtered_Koh.rds: Rscripts/import_datasets/import_QC_Koh.Rmd \
data/data_raw/SRP073808.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

data/sce_filtered/sce_filtered_Zhengmix.rds: Rscripts/import_datasets/import_QC_Zhengmix.Rmd \
data/data_raw/zheng/b_cells_filtered/hg19/matrix.mtx data/data_raw/zheng/cd14_monocytes_filtered/hg19/matrix.mtx \
data/data_raw/zheng/naive_cytotoxic_filtered/hg19/matrix.mtx data/data_raw/zheng/regulatory_t_filtered/hg19/matrix.mtx
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

define qckumarrule
data/sce_filtered/sce_filtered_$(1).rds: Rscripts/import_datasets/import_QC_$(1).Rmd \
data/data_raw/GSE60749-GPL13112.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$$(<F)', clean = TRUE)"
endef
$(foreach d,Kumar SimKumar,$(eval $(call qckumarrule,$(d))))