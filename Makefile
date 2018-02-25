## Define the versions of R and the paths to the libraries
R := R_LIBS=/home/Shared/Rlib/release-3.6-lib/ /usr/local/R/R-3.4.2/bin/R CMD BATCH --no-restore --no-save
Rscript := /usr/local/R/R-3.4.2/src/unix/Rscript --no-restore --no-save

## Include lists of methods, data sets and gene filtering approaches to use
include include_methods.mk
include include_datasets.mk
include include_filterings.mk

.PHONY: all prepare_data cluster

## Define rules
all: cluster

## Prepare data
prepare_data: $(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach p,$(PCTKEEP),data/sce_filtered$(f)$(p)/sce_filtered$(f)$(p)_$(d).rds)))

#$(foreach d,$(DATASETS),plots/qc_data/$(d).rds)

cluster: $(foreach f,$(ALLFILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),results/sce_$(f)_$(d)_$(m).rds)))

plots: plots/runtime/runtime_by_k.rds plots/performance/performance_by_k.rds

## ------------------------------------------------------------------------------------ ##
## Setup
## Note! Run the setup rule only once, to set up the directory structure and download the data. 
## If it is rerun, it will update the raw data files and cause the results to be out of date.
## If additional parameters need to be defined or modified at a later stage, generate only 
## those parameter files manually.
## ------------------------------------------------------------------------------------ ##
## Set up directory structure, download raw data and define parameter settings
setup: directories zheng conquer
	$(R) Rscripts/parameter_settings/generate_parameter_settings.R Rout/generate_parameter_settings.Rout

directories:
	mkdir -p data/data_raw
	mkdir -p data/data_raw/zheng
	mkdir -p data/sce_full
	mkdir -p plots/qc_data
	mkdir -p plots/performance
	mkdir -p plots/runtime
	mkdir -p results
	mkdir -p Rout	

conquer: directories
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE60749-GPL13112.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE52529-GPL16791.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/SRP073808.rds

zheng: directories
	for D in b_cells naive_cytotoxic cd14_monocytes regulatory_t cd56_nk memory_t cd4_t_helper naive_t; do \
		mkdir -p data/data_raw/zheng/$${D}; \
		wget -P data/data_raw/zheng http://cf.10xgenomics.com/samples/cell-exp/1.1.0/$${D}/$${D}_filtered_gene_bc_matrices.tar.gz; \
		tar zxvf data/data_raw/zheng/$${D}_filtered_gene_bc_matrices.tar.gz -C data/data_raw/zheng/$${D}; \
		rm -f data/data_raw/zheng/$${D}_filtered_gene_bc_matrices.tar.gz; \
	done

## ------------------------------------------------------------------------------------ ##
## Prepare data sets
## ------------------------------------------------------------------------------------ ##
data/sce_full/sce_full_Trapnell.rds: Rscripts/import_datasets/import_QC_Trapnell.Rmd \
data/data_raw/GSE52529-GPL16791.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

data/sce_full/sce_full_Koh.rds: Rscripts/import_datasets/import_QC_Koh.Rmd \
data/data_raw/SRP073808.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

define zheng4rule
data/sce_full/sce_full_$(1).rds: Rscripts/import_datasets/import_QC_$(1).Rmd \
data/data_raw/zheng/b_cells/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/cd14_monocytes/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/naive_cytotoxic/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/regulatory_t/filtered_matrices_mex/hg19/matrix.mtx
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$$(<F)', clean = TRUE)"
endef
$(foreach d,Zhengmix4eq Zhengmix4uneq,$(eval $(call zheng4rule,$(d))))

define zheng8rule
data/sce_full/sce_full_$(1).rds: Rscripts/import_datasets/import_QC_$(1).Rmd \
data/data_raw/zheng/b_cells/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/cd14_monocytes/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/naive_cytotoxic/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/regulatory_t/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/cd4_t_helper/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/cd56_nk/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/memory_t/filtered_matrices_mex/hg19/matrix.mtx \
data/data_raw/zheng/naive_t/filtered_matrices_mex/hg19/matrix.mtx
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$$(<F)', clean = TRUE)"
endef
$(foreach d,Zhengmix8eq,$(eval $(call zheng8rule,$(d))))

define qckumarrule
data/sce_full/sce_full_$(1).rds: Rscripts/import_datasets/import_QC_$(1).Rmd \
data/data_raw/GSE60749-GPL13112.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$$(<F)', clean = TRUE)"
endef
$(foreach d,Kumar SimKumar4easy SimKumar4hard SimKumar8hard,$(eval $(call qckumarrule,$(d))))

## ------------------------------------------------------------------------------------ ##
## Generate filtered data sets
## ------------------------------------------------------------------------------------ ##
define filterrule
data/sce_filtered$(2)$(3)/sce_filtered$(2)$(3)_$(1).rds: data/sce_full/sce_full_$(1).rds Rscripts/filtering/filter$(2).R
	mkdir -p $$(@D)
	$(R) "--args scefile='$$<' method='$(2)' pctkeep=$(3) outrds='$$@'" Rscripts/filtering/filter_genes.R Rout/filter_genes_$(1)_$(2)$(3).Rout
endef
$(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach p,$(PCTKEEP),$(eval $(call filterrule,$(d),$(f),$(p))))))

## ------------------------------------------------------------------------------------ ##
## Generate QC plots
## ------------------------------------------------------------------------------------ ##
define qcrule
plots/qc_data/$(1).rds: data/sce_full/sce_full_$(1).rds Rscripts/evaluate_datasets/plot_dataset_characteristics.R
	$(R) "--args scefull='data/sce_full/sce_full_$(1).rds' scefiltered='data/sce_filteredExpr/sce_filteredExpr_$(1).rds' outrds='$$@'" Rscripts/evaluate_datasets/plot_dataset_characteristics.R Rout/plot_dataset_characteristics_$(1).Rout
endef
$(foreach d,$(DATASETS),$(eval $(call qcrule,$(d))))

## ------------------------------------------------------------------------------------ ##
## Apply clustering methods
## ------------------------------------------------------------------------------------ ##
define clusterrule ## $(1) - sce_full, sce_filteredExpr. $(2) - dataset. $(3) - clustering method
results/$(1)_$(2)_$(3).rds: data/sce_full/sce_full_$(2).rds Rscripts/clustering/run_clustering.R \
Rscripts/clustering/apply_$(3).R parameter_settings/$(1)_$(2).json parameter_settings/$(1)_$(2)_$(3).json \
parameter_settings/$(3).json
	$(R) "--args scefile='data/$(1)/$(1)_$(2).rds' method='$(3)' outrds='results/$(1)_$(2)_$(3).rds'" Rscripts/clustering/run_clustering.R Rout/run_clustering_$(1)_$(2)_$(3).Rout
endef
$(foreach f,$(FILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),$(eval $(call clusterrule,sce_$(f),$(d),$(m))))))

## ------------------------------------------------------------------------------------ ##
## Plot performance
## ------------------------------------------------------------------------------------ ##
plots/performance/performance_by_k.rds: Rscripts/evaluate_results/plot_performance_by_k.R \
$(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach m,$(METHODS),results/sce_$(f)_$(d)_$(m).rds)))
	$(R) "--args datasets='$(DATASETSc)' filterings='$(FILTERINGSc)' methods='$(METHODSc)' outrds='$@'" Rscripts/evaluate_results/plot_performance_by_k.R Rout/plot_performance_by_k.Rout

## ------------------------------------------------------------------------------------ ##
## Plot runtime
## ------------------------------------------------------------------------------------ ##
plots/runtime/runtime_by_k.rds: Rscripts/evaluate_results/plot_runtime_by_k.R \
$(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach m,$(METHODS),results/sce_$(f)_$(d)_$(m).rds)))
	$(R) "--args datasets='$(DATASETSc)' filterings='$(FILTERINGSc)' methods='$(METHODSc)' outrds='$@'" Rscripts/evaluate_results/plot_runtime_by_k.R Rout/plot_runtime_by_k.Rout




