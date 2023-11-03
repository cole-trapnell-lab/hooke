# Overview
        
## Installing Hooke

Hooke depends on the PLNmodels package. 

Use the github install: 

`devtools::install_github("cole-trapnell-lab/PLNmodels)`

See the [PLN website](https://pln-team.github.io/PLNmodels/index.html) for more details. 

`devtools::install_github("cole-trapnell-lab/hooke)`

See our [Github repository](https://github.com/cole-trapnell-lab/hooke) for more details. 


## Using a Seurat object


Currently Hooke only supports Monocle3 cell data set objects. If using Seurat, please see the [Seurat documentation](https://satijalab.org/seurat/reference/as.celldataset) on how to convert a Seurat object to a Monocle3 object. 

## Example datasets

##### Processed data subsets:


[Silicosis data (.RDS)](https://depts.washington.edu/trapnell-lab/software/hooke/silicosis_cds.rds){ .md-button }


This studied was published [Hasegawa, Franks, et al. _bioRxiv_, (2023)](https://www.biorxiv.org/content/10.1101/2023.02.17.528996v1)




Our findings are published in [Saunders, Srivatsan, et al. _Nature_, in press (2023)](https://www.biorxiv.org/content/10.1101/2022.08.04.502764v1)

For more information about this dataset, see the [ZSCAPE website](https://cole-trapnell-lab.github.io/zscape/)

* Subclustered and finely annotated sensory cranial ganglia cells (plus Rohon-Beard neurons). This dataset contains just wild type and control-injected cells.  

[Cranial ganglia subset (.RDS)](https://depts.washington.edu/trapnell-lab/software/hooke/all-geno_sensory-cranial-ganglion_neuron_29k_cds.RDS){ .md-button }


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
