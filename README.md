# Hooke

Hooke is a new software package that uses Poisson-Lognormal models to perform differential analysis of cell abundances for perturbation experiments read out by single-cell RNA-seq. This versatile framework allows users to both perform multivariate statistical regression to describe how perturbations alter the relative abundances of each cell state and visualize cell type abundance kinetics.

See the [Hooke website](https://cole-trapnell-lab.github.io/hooke/) for a more comprehensive introduction. 

## Installation

Hooke runs in the [R statistical computing environment](https://www.r-project.org/). It requires R >= 3.5.0. Hooke is currently only available for Github install. 

### Required software

Hooke builds on top of the [Monocle3 package](https://cole-trapnell-lab.github.io/monocle3/docs/installation/). 

```r
devtools::install_github("cole-trapnell-lab/monocle3")
```

Hooke depends on the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html). 

Use the github install: 

```r 
remotes::install_github("pln-team/PLNmodels")
```

Finally, install the hooke package as follows: 

```r
devtools::install_github("cole-trapnell-lab/hooke")
```
