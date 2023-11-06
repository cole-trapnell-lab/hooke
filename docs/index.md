# Overview
        
Hooke is a new software package that uses Poisson-Lognormal models to perform differential analysis of cell abundances for perturbation experiments read out by single-cell RNA-seq. This versatile framework allows users to both 1) perform multivariate statistical regression to describe how perturbations alter the relative abundances of each cell state and 2) easily performed pseudobulked differential gene expression analysis. 

## Installation

Hooke runs in the [R statistical computing environment](https://www.r-project.org/). Hooke is currently only available for Github install. 

### Required software

Hooke builds on top of the [Monocle3 package](https://cole-trapnell-lab.github.io/monocle3/docs/installation/). 

```devtools::install_github("cole-trapnell-lab/monocle3, ref="develop")```

Hooke depends on the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html). Currently we are using a forked version of the PLNmodels package until our pull request is approved. 

Use the github install: 

```devtools::install_github("cole-trapnell-lab/PLNmodels)```

Finally, install the hooke package as follows: 

```devtools::install_github("cole-trapnell-lab/hooke)```

See our [Github repository](https://github.com/cole-trapnell-lab/hooke) for more details.

**_NOTE:_** Hooke is currently in the beta phase of its development. The documentation on this page is also still under construction. Not all features currently implemented have been completely documented. Please report any issues to your [github page](https://github.com/cole-trapnell-lab/hooke/issues). 


## Data type requirements

Hooke is built for experiments with multiple treatments and replicates, taking advantage of replicates in various groups or perturbations. As input, Hooke takes in a cell x gene matrix where cells are annotated according to type (or cluster) and by which sample or specimen they came from. It will aggregate cells according to type and by sample. This collapses the matrix into a new, smaller matrix where rows are cell types and the columns denote how many cells of that type were present in each sample. We refer to this as a *cell type abundance matrix*. 

![aggregation_cells](assets/aggregation_cells.png)


## Using a Seurat object


Currently Hooke only supports Monocle3 cell data set objects. If using Seurat, please see the [Seurat documentation](https://satijalab.org/seurat/reference/as.celldataset) on how to convert a Seurat object to a Monocle3 object. 

## Example datasets

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
