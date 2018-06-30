## Define the versions of R and the paths to the libraries
#R := R_LIBS=/home/Shared/Rlib/release-3.6-lib/ /usr/local/R/R-3.4.2/bin/R CMD BATCH --no-restore --no-save
#Rd := R_LIBS=/home/Shared/Rlib/devel-lib/ /usr/local/R/R-devel/bin/R CMD BATCH --no-restore --no-save
Rscript := /usr/local/R/R-3.4.2/src/unix/Rscript --no-restore --no-save
R := R_LIBS=Rlibrary3.5 /usr/local/R/R-3.5.0/bin/R CMD BATCH --no-restore --no-save

## Include lists of methods, data sets and gene filtering approaches to use
include include_methods.mk
include include_datasets.mk
include include_filterings.mk

ncores := 12

.PHONY: all prepare_data cluster

## Define rules
all: cluster

## Prepare data
prepare_data: $(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach p,$(PCTKEEP),data/sce_filtered$(f)$(p)/sce_filtered$(f)$(p)_$(d).rds)))

cluster: $(foreach f,$(ALLFILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),results/sce_$(f)_$(d)_$(m).rds)))

cluster5: $(foreach f,$(ALLFILTERINGS),$(foreach m,RaceID2,$(foreach d,$(DATASETS),results/sce_$(f)_$(d)_$(m).rds)))

summarise: output/consensus/consensus.rds output/ensemble/ensemble.rds output/silhouettes/silhouettes.rds

figs: plots/manuscript/figure1.rds plots/manuscript/figure2.rds plots/manuscript/figure3.rds \
plots/performance/seurat_diagnostics.rds \
plots/performance/res_performance_cons.rds plots/ensemble/ensemble_vs_individual.rds \
plots/similarities_between_methods/similarities.rds plots/shared_genes_filterings/shared_genes_filterings.rds \
plots/facets_clustering/facets_clustering.rds

memoryusage: plots/memory_usage/memory_usage.rds

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
	mkdir -p output/clustering_summary
	mkdir -p output/parameter_range_investigation

conquer: directories
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE60749-GPL13112.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE52529-GPL16791.rds
	wget -P data/data_raw http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/SRP073808.rds
	wget -O data/data_raw/SRP073808TCC.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/data-tcc/SRP073808.rds
	wget -O data/data_raw/GSE52529-GPL16791TCC.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/data-tcc/GSE52529-GPL16791.rds
	wget -O data/data_raw/GSE60749-GPL13112TCC.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/data-tcc/GSE60749-GPL13112.rds

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

data/sce_full/sce_full_TrapnellTCC.rds: Rscripts/import_datasets/import_QC_TrapnellTCC.Rmd \
data/data_raw/GSE52529-GPL16791TCC.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

data/sce_full/sce_full_Koh.rds: Rscripts/import_datasets/import_QC_Koh.Rmd \
data/data_raw/SRP073808.rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$(<F)', clean = TRUE)"

data/sce_full/sce_full_KohTCC.rds: Rscripts/import_datasets/import_QC_KohTCC.Rmd \
data/data_raw/SRP073808TCC.rds
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
data/sce_full/sce_full_$(1)$(2).rds: Rscripts/import_datasets/import_QC_$(1)$(2).Rmd \
data/data_raw/GSE60749-GPL13112$(2).rds
	cd Rscripts/import_datasets && \
	$(Rscript) -e "rmarkdown::render('$$(<F)', clean = TRUE)"
endef
$(foreach d,Kumar SimKumar4easy SimKumar4hard SimKumar8hard,$(eval $(call qckumarrule,$(d),)))
$(foreach d,Kumar,$(eval $(call qckumarrule,$(d),TCC)))

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
#define qcrule
#plots/qc_data/$(1).rds: data/sce_full/sce_full_$(1).rds Rscripts/evaluate_datasets/plot_dataset_characteristics.R
#	$(R) "--args scefull='data/sce_full/sce_full_$(1).rds' scefiltered='data/sce_filteredExpr/sce_filteredExpr_$(1).rds' outrds='$$@'" Rscripts/evaluate_datasets/plot_dataset_characteristics.R Rout/plot_dataset_characteristics_$(1).Rout
#endef
#$(foreach d,$(DATASETS),$(eval $(call qcrule,$(d))))

## ------------------------------------------------------------------------------------ ##
## Apply clustering methods
## ------------------------------------------------------------------------------------ ##
define clusterrule ## $(1) - sce_full, sce_filteredExpr, sce_filteredHVG, sce_filteredM3Drop. $(2) - dataset. $(3) - clustering method
results/$(1)_$(2)_$(3).rds: data/$(1)/$(1)_$(2).rds parameter_settings/$(1)_$(2).json \
parameter_settings/$(1)_$(2)_$(3).json parameter_settings/$(3).json \
Rscripts/clustering/apply_$(3).R Rscripts/clustering/run_clustering.R
	$(4) "--args scefile='data/$(1)/$(1)_$(2).rds' method='$(3)' outrds='results/$(1)_$(2)_$(3).rds'" Rscripts/clustering/run_clustering.R Rout/run_clustering_$(1)_$(2)_$(3).Rout
endef
$(foreach f,$(ALLFILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),$(eval $(call clusterrule,sce_$(f),$(d),$(m),$(R))))))

## ------------------------------------------------------------------------------------ ##
## Summarize clustering performance
## ------------------------------------------------------------------------------------ ##
output/clustering_summary/clustering_summary.rds: Rscripts/evaluate_results/summarize_clustering_results.R \
$(foreach f,$(ALLFILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),results/sce_$(f)_$(d)_$(m).rds)))
	mkdir -p $(@D)
	$(R) "--args datasets='$(DATASETSc)' filterings='$(ALLFILTERINGSc)' methods='$(METHODSc)' outrds='$@'" Rscripts/evaluate_results/summarize_clustering_results.R Rout/summarize_clustering_results.Rout

## ------------------------------------------------------------------------------------ ##
## Compute consensus
## ------------------------------------------------------------------------------------ ##
output/consensus/consensus.rds: output/clustering_summary/clustering_summary.rds Rscripts/similarities_consensus/compute_consensus.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' ncores=$(ncores) outrds='$@'" Rscripts/similarities_consensus/compute_consensus.R Rout/compute_consensus.Rout

## ------------------------------------------------------------------------------------ ##
## Compute ensembles
## ------------------------------------------------------------------------------------ ##
output/ensemble/ensemble.rds: output/clustering_summary/clustering_summary.rds Rscripts/ensemble/compute_ensemble.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' ncores=$(ncores) outrds='$@'" Rscripts/ensemble/compute_ensemble.R Rout/compute_ensemble.Rout

## ------------------------------------------------------------------------------------ ##
## Compute silhouette widths for full data sets
## ------------------------------------------------------------------------------------ ##
output/silhouettes/silhouettes.rds: $(foreach d,$(DATASETS),data/sce_full/sce_full_$(d).rds) \
Rscripts/evaluate_datasets/compute_avg_silhouette_width.R
	mkdir -p $(@D)
	$(R) "--args datadir='data' datasets='$(DATASETSc)' ncores=$(ncores) outrds='$@'" Rscripts/evaluate_datasets/compute_avg_silhouette_width.R Rout/compute_avg_silhouette_width.Rout

## ------------------------------------------------------------------------------------ ##
## Plots
## ------------------------------------------------------------------------------------ ##
## performance
plots/performance/performance.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_performance.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_performance.R Rout/plot_performance.Rout

## timing
plots/runtime/runtime.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_timing.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_timing.R Rout/plot_timing.Rout

## stability
plots/performance/stability.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_stability.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_stability.R Rout/plot_stability.Rout

## entropy
plots/performance/entropy.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_entropy.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_entropy.R Rout/plot_entropy.Rout

## difference between estimated or optimal k and true k
plots/performance/difference_in_k.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_k_differences.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_k_differences.R Rout/plot_k_differences.Rout

## Seurat k vs resolution
plots/performance/seurat_diagnostics.rds: output/clustering_summary/clustering_summary.rds \
Rscripts/evaluate_results/plot_Seurat_k_resolution.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_Seurat_k_resolution.R Rout/plot_Seurat_k_resolution.Rout

## similarities between methods
plots/similarities_between_methods/similarities.rds: output/consensus/consensus.rds \
Rscripts/similarities_consensus/plot_similarities_between_methods2.R 
	mkdir -p $(@D)
	$(R) "--args consensusrds='$<' outrds='$@'" Rscripts/similarities_consensus/plot_similarities_between_methods2.R Rout/plot_similarities_between_methods2.Rout

## performance for consensus clusters
plots/performance/res_performance_cons.rds: output/consensus/consensus.rds \
Rscripts/similarities_consensus/plot_evaluate_performance_cons.R Rscripts/Colorscheme.R
	mkdir -p $(@D)
	$(R) "--args consensusrds='$<' outrds='$@'" Rscripts/similarities_consensus/plot_evaluate_performance_cons.R Rout/plot_evaluate_performance_cons.Rout

## ensembles vs individual methods
plots/ensemble/ensemble_vs_individual.rds: output/clustering_summary/clustering_summary.rds \
output/ensemble/ensemble.rds Rscripts/ensemble/plot_ensemble_vs_individual.R
	mkdir -p $(@D)
	$(R) "--args clusteringsummary='$<' ensemblerds='$(word 2,$^)' outrds='$@'" Rscripts/ensemble/plot_ensemble_vs_individual.R Rout/plot_ensemble_vs_individual.Rout

plots/shared_genes_filterings/shared_genes_filterings.rds: \
$(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach p,$(PCTKEEP),data/sce_filtered$(f)$(p)/sce_filtered$(f)$(p)_$(d).rds))) \
Rscripts/evaluate_datasets/plot_shared_genes_venn.R
	mkdir -p $(@D)
	$(R) "--args datadir='data' datasets='$(DATASETSc)' filterings='$(FILTERINGSc)' pctkeep=10 outrds='$@'" Rscripts/evaluate_datasets/plot_shared_genes_venn.R Rout/plot_shared_genes_venn.Rout

plots/facets_clustering/facets_clustering.rds: output/clustering_summary/clustering_summary.rds \
$(foreach d,$(DATASETS),$(foreach f,$(FILTERINGS),$(foreach p,$(PCTKEEP),data/sce_filtered$(f)$(p)/sce_filtered$(f)$(p)_$(d).rds))) \
Rscripts/evaluate_results/plot_facets_of_clustering_tSNE.R
	mkdir -p $(@D)
	$(R) "--args datadir='data' datasets='$(DATASETSc)' filterings='$(FILTERINGSc)' pctkeep=10 clusteringsummary='$<' outrds='$@'" Rscripts/evaluate_results/plot_facets_of_clustering_tSNE.R Rout/plot_facets_of_clustering_tSNE.Rout





## ------------------------------------------------------------------------------------ ##
## Manuscript figures
## ------------------------------------------------------------------------------------ ##
plots/manuscript/figure1.rds: plots/performance/performance.rds Rscripts/manuscript/plot_figure1.R
	mkdir -p $(@D)
	$(R) "--args performancerds='$<' outrds='$@'" Rscripts/manuscript/plot_figure1.R Rout/plot_figure1.Rout

plots/manuscript/figure2.rds: plots/performance/stability.rds plots/performance/entropy.rds \
plots/performance/difference_in_k.rds plots/runtime/runtime.rds Rscripts/manuscript/plot_figure2.R
	mkdir -p $(@D)
	$(R) "--args stabilityrds='$(word 1,$^)' entropyrds='$(word 2,$^)' diffrds='$(word 3,$^)' timerds='$(word 4,$^)' outrds='$@'" Rscripts/manuscript/plot_figure2.R Rout/plot_figure2.Rout

plots/manuscript/figure3.rds: plots/ensemble/ensemble_vs_individual.rds Rscripts/manuscript/plot_figure3.R
	mkdir -p $(@D)
	$(R) "--args ensemblerds='$(word 1,$^)' outrds='$@'" Rscripts/manuscript/plot_figure3.R Rout/plot_figure3.Rout


















## ------------------------------------------------------------------------------------ ##
## Plot memory usage
## ------------------------------------------------------------------------------------ ##
plots/memory_usage/memory_usage.rds: Rscripts/plot_memory_usage.R \
$(foreach f,$(ALLFILTERINGS),$(foreach m,$(METHODS),$(foreach d,$(DATASETS),results/sce_$(f)_$(d)_$(m).rds)))
	mkdir -p $(@D)
	grep "Ncells|Vcells" Rout/* > $(@D)/memory_usage.txt
	$(R) "--args memusetxt='memory_usage.txt' outrds='$@'" Rscripts/plot_memory_usage.R Rout/plot_memory_usage.Rout

