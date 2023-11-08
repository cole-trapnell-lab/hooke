# Overview
        
Hooke is a new software package that uses Poisson-Lognormal models to perform differential analysis of cell abundances for perturbation experiments read out by single-cell RNA-seq. This versatile framework allows users to both perform multivariate statistical regression to describe how perturbations alter the relative abundances of each cell state and visualize cell type abundance kinetics. 

![website_figure](assets/hooke_website_figure.png)


## Data type requirements

Hooke is built for experiments with multiple treatments and replicates, taking advantage of replicates in various groups or perturbations. As input, Hooke takes in a cell x gene matrix where cells are annotated according to type (or cluster) and by which sample or specimen they came from. It will aggregate cells according to type and by sample. This collapses the matrix into a new, smaller matrix where rows are cell types and the columns denote how many cells of that type were present in each sample. We refer to this as a *cell type abundance matrix*. 

![aggregation_cells](assets/aggregation_example_cells.png){width=75%}


## Using a Seurat object


Currently Hooke only supports Monocle3 cell data set objects. If using Seurat, please see the [Seurat documentation](https://satijalab.org/seurat/reference/as.celldataset) on how to convert a Seurat object to a Monocle3 object. 

