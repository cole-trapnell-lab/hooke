# Installation

Hooke runs in the [R statistical computing environment](https://www.r-project.org/). It requires R >= 3.5.0. Hooke is currently only available for Github install. 

## Required software

Hooke builds on top of the [Monocle3 package](https://cole-trapnell-lab.github.io/monocle3/docs/installation/). 

```r 
devtools::install_github("cole-trapnell-lab/monocle3")
```

Hooke depends on the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html). Currently we are using a forked version of the PLNmodels package until our pull request is approved. 

Use the github install: 

```r
devtools::install_github("cole-trapnell-lab/PLNmodels")
```

Finally, install the hooke package as follows: 

```r
devtools::install_github("cole-trapnell-lab/hooke")
```

See our [Github repository](https://github.com/cole-trapnell-lab/hooke) for more details.

> **_NOTE:_** Hooke is currently in the beta phase of its development. The documentation on this page is also still under construction. Not all features currently implemented have been completely documented. Please report any issues to your [github page](https://github.com/cole-trapnell-lab/hooke/issues). 