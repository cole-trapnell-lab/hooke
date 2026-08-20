# Installation

Hooke runs in the [R statistical computing environment](https://www.r-project.org/). It requires R >= 3.5.0. Hooke is currently only available for Github install. 

> **_NOTE:_** Hooke is currently in the beta phase of its development. The documentation on this page is also still under construction. Not all features currently implemented have been completely documented. Please report any issues to your [github page](https://github.com/cole-trapnell-lab/hooke/issues). 

## Required software

Hooke builds on top of the [Monocle3 package](https://cole-trapnell-lab.github.io/monocle3/docs/installation/). 

```r 
devtools::install_github("cole-trapnell-lab/monocle3")
```

Hooke depends on the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html), which is available on CRAN. 
You can install as follows: 

```r
install.packages("PLNmodels")
```

Hooke also imports several [Bioconductor](https://www.bioconductor.org/) packages (`Biobase`, `batchelor`, `BiocGenerics`, `S4Vectors`, `Rgraphviz`), so you will need `BiocManager` to install those if they aren't already present: 

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biobase", "batchelor", "BiocGenerics", "S4Vectors", "Rgraphviz"))
```

Finally, install the hooke package as follows: 

```r
devtools::install_github("cole-trapnell-lab/hooke")
```

If you wish to install the develop branch of hooke, execute: 

```r
devtools::install_github("cole-trapnell-lab/hooke", ref="develop")
```

See our [Github repository](https://github.com/cole-trapnell-lab/hooke) for more details.

