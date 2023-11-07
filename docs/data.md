# Example datasets

Descriptions and download links for the datasets used in the Hooke vignettes.  

#### Silica-induced pulmonary fibrosis [mouse]:

[Silicosis data (.RDS)](https://depts.washington.edu/trapnell-lab/software/hooke/silicosis_cds.rds){ .md-button }

Subclustered and finely annotated whole-lung single-nucleas RNA sequencing from a silica-induced pulmonary fibrosis mouse model. 

This studied was published in [Hasegawa, Franks, et al. _bioRxiv_, (2023)](https://www.biorxiv.org/content/10.1101/2023.02.17.528996v1)

##### Cell metadata breakdown 

Important columns in this data: 

* `ID` - individual sample ID. 
* `Timepoint` - timepoint the lung was sampled. 
* `Rep` - replicate
* `fine_annotation` - The finest cell type annotation level.
* `broad_annotation` - A broader level of cell type annotation than `cell_type_sub`, but still capturing all uniquely identified cell types.

#### Cranial sensory ganglia [zebrafish]:

Subclustered and finely annotated sensory cranial ganglia cells (plus Rohon-Beard neurons). This dataset contains just wild type and control-injected cells.  

[Cranial ganglia subset (.RDS)](https://depts.washington.edu/trapnell-lab/software/hooke/all-geno_sensory-cranial-ganglion_neuron_29k_cds.RDS){ .md-button }

Our findings are published in [Saunders, Srivatsan, et al. _Nature_, in press (2023)](https://www.biorxiv.org/content/10.1101/2022.08.04.502764v1). For more information about this dataset, see the [ZSCAPE website](https://cole-trapnell-lab.github.io/zscape/). 

##### Cell metadata breakdown 

Our cell metadata contains lots of information about time, perturbation, statistical metrics, and annotation. Here is a breakdown of those attributes according to our column names:

* `timepoint`: The developmental stage in hours post fertilization (hpf) of embryos. Embryos were staged according to key landmarks according to [(Kimmel, et al (1995))](https://zfin.org/zf_info/zfbook/stages/index.html). Options are 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 72, and 96 hpf.
* `expt`: Denotes unique library preparation instances.
* `cell_type_sub`: The finest cell type annotation level.
* `cell_type_broad`: A broader level of cell type annotation than `cell_type_sub`, but still capturing all uniquely identified cell types.
* `tissue`: The tissue type annotation that contains the cell types.
* `germ_layer`: The germ layer of origin for each cell type, if known.
* `gene_target`: The genetic perturbation target. Controls include injected (scrambled) and unjected. Sibling controls are listed as `ctrl-<target>`, for null mutants included in the study.
* `mean_nn_time`: The mean time points of the nearest 15 neighbor cells
* `embryo`: The individual embryo barcode.
* `temp`: The growth temperature of the embryos (standard conditions are 28C).
* `pseudostage`: Embryo-level staging prediction by cell composition.
