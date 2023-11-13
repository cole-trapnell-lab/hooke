<!-- toc -->

November 13, 2023

# DESCRIPTION

```
Package: hooke
Title: Differential analysis of cell counts in single-cell experiments
Version: 0.0.1
Authors@R: 
    person(given = "Madelein",
           family = "Duran",
           role = c("aut", "cre"),
           email = "duran@uw.edu")
    person(given = "Cole",
           family = "Trapnell",
           role = c("aut", "cre"),
           email = "coletrap@uw.edu")
Description: Hooke models abundances of different cell types and how they change in single-cell genomics experiments.
License: MIT + file LICENSE
Encoding: UTF-8
LazyData: true
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.1
LinkingTo: 
    Rcpp
Depends:
    Biobase, 
    monocle3,
    PLNmodels
Imports:
    assertthat (>= 0.2.1),
    batchelor,
    BiocGenerics (>= 0.28.0),
    dplyr (>= 0.8.0.1),
    ggplot2 (>= 3.1.1),
    ggrepel (>= 0.8.1),
    grr,
    igraph (>= 1.2.4),
    Matrix (>= 1.2-17),
    methods,
    purrr (>= 0.3.2),
    RColorBrewer,
    Rcpp (>= 1.0.1),
    S4Vectors,
    stringr (>= 1.4.0),
    tibble (>= 2.1.1),
    tidyr (>= 0.8.3),
    viridis (>= 0.5.1),
    maxmatching,
    Rgraphviz,
    ggnetwork,
    tidygraph,
    ggforce,
    ggnewscale
Suggests: 
    testthat (>= 2.1.0),
    pryr (>= 0.1.4),
    knitr,
    rmarkdown,
    spelling
VignetteBuilder: knitr
Language: en-US```


# `add_umap_coords`

adds umap coords to a data frame


## Description

adds umap coords to a data frame


## Usage

```r
add_umap_coords(df, umap_centers)
```


# `blacklist`

return the complementing edges


## Description

return the complementing edges


## Usage

```r
blacklist(edges)
```


# `calc_max_flow`

Calculate max flow between two points


## Description

Calculate max flow between two points


## Usage

```r
calc_max_flow(edges, source, target)
```


# `calc_mst`

Calculate a minimum spanning tree


## Description

Calculate a minimum spanning tree


## Usage

```r
calc_mst(edges, weight = "pcor")
```


# `calc_shortest_path`

Calculate the shortest path between two points


## Description

Calculate the shortest path between two points


## Usage

```r
calc_shortest_path(G, from, to)
```


## Arguments

Argument      |Description
------------- |----------------
`edges`     |     data frame of edges with edge weights


## Value

data frame containing the shortest path


# `cds`

Get the underlying cell_data_set object from a cell_count_model.


## Description

Get the underlying cell_data_set object from a cell_count_model.


## Usage

```r
cds(ccm)
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model object.


## Value

A cell_data_set object


## Examples

```r
cds(ccm)
```


# `cell_count_model`

The cell_count_model class


## Description

The main class used by Monocle3 to hold single-cell expression data.
 cell_count_model extends the Bioconductor SingleCellExperiment class.


## Details

This class is initialized from a matrix of expression values along with cell
 and feature metadata.


# `cell_count_set`

The cell_count_set class


## Description

The main class used by Hooke to hold cell abundances data.
 cell_count_set extends the Monocle's cell_data_set class.


## Details

This class is initialized from a matrix of expression values along with cell
 and feature metadata.


# `centroids`

Get the centroids of cell groups in UMAP/PCA space.


## Description

Get the centroids of cell groups in UMAP/PCA space.


## Usage

```r
centroids(ccs, reduction_method = "UMAP", switch_group = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`ccs`     |     A cell_count_set object.


## Value

A data frame of centroid coordinates


# `collect_pln_graph_edges`

Orient graph edges from a PLNnetwork using a contrast between conditions


## Description

Orient graph edges from a PLNnetwork using a contrast between conditions


## Usage

```r
collect_pln_graph_edges(
  ccm,
  cond_b_vs_a_tbl,
  log_abundance_thresh = 1 - 5,
  model_for_pcors = "reduced"
)
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model
`cond_b_vs_a_tbl`     |     A contrast between two conditions as returned by compare_abundances()


# `compare_abundances`

Compare two estimates of cell abundances from a Hooke model


## Description

Compare two estimates of cell abundances from a Hooke model


## Usage

```r
compare_abundances(ccm, cond_x, cond_y, method = "BH")
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model
`cond_x`     |     An estimate from estimate_abundances()
`cond_y`     |     An estimate from estimate_abundances()


## Value

A table contrasting cond_x and cond_y (interpret as Y/X)


# `estimate_abundances`

Predict cell type abundances given a PLN model and a set of inputs for its covariates


## Description

Predict cell type abundances given a PLN model and a set of inputs for its covariates


## Usage

```r
estimate_abundances(ccm, newdata, min_log_abund = -5)
```


## Arguments

Argument      |Description
------------- |----------------
`newdata`     |     needs to be suitable input to pln_model


# `get_paga_graph`

Extract a partitioned abstract graph from a Monocle cell_data_set object


## Description

QUESTION: Do the cells need to be grouped by monocle cluster? Or can they
 be grouped arbitrarily (e.g. by cell type? Would be good to generalize if
 possible. For now, only works when cells are grouped by Monocle cluster


## Usage

```r
get_paga_graph(cds, reduction_method = "UMAP")
```


## Arguments

Argument      |Description
------------- |----------------
`cds`     |     A cell_data_set object. cluster_cells() must have been called.
`reduction_method`     |     The coordinate space in which to build the graph


# `get_pcor_edges`

returns edges based on pcor values


## Description

returns edges based on pcor values


## Usage

```r
get_pcor_edges(ccm, selected_model = c("reduced", "full"))
```


# `hello`

Hello, World!


## Description

Prints 'Hello, world!'.


## Usage

```r
hello()
```


## Examples

```r
hello()
```


# `init_penalty_matrix`

Initialize the PLN network penalty matrix, accepting optional whitelists and
 blacklists of edges that are "free" or "off limits" between cell groups


## Description

Initialize the PLN network penalty matrix, accepting optional whitelists and
 blacklists of edges that are "free" or "off limits" between cell groups


## Usage

```r
init_penalty_matrix(
  ccs,
  whitelist = NULL,
  blacklist = NULL,
  base_penalty = 1,
  min_penalty = 0.01,
  max_penalty = 1e+06
)
```


## Arguments

Argument      |Description
------------- |----------------
`ccs`     |     A cell_count_set of aggregated cell counts
`whitelist`     |     a data frame with two columns corresponding to (undirected) edges that should receive no penalty
`blacklist`     |     a data frame with two columns corresponding to (undirected) edges that should receive very high penalty
`dist_fun`     |     A function that returns a penalty based given a distance between two clusters


# `model`

Get the selected PLNnetwork from a cell_count_model object.


## Description

Get the selected PLNnetwork from a cell_count_model object.


## Usage

```r
model(ccm, model_to_return = c("full", "reduced"))
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model object.


## Value

An updated cell_count_model object


## Examples

```r
model(ccm)
```


# `new_cell_count_model`

Create a new cell_count_model object.


## Description

Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a
 way of encoding a graph prior. Automatically selects sparsity parameter, but
 allows user to update it.


## Usage

```r
new_cell_count_model(
  ccs,
  main_model_formula_str,
  nuisance_model_formula_str = "1",
  penalty_matrix = NULL,
  whitelist = NULL,
  blacklist = NULL,
  sparsity_factor = 0.1,
  base_penalty = 1,
  min_penalty = 0.01,
  max_penalty = 1e+06,
  verbose = FALSE,
  pseudocount = 0,
  pln_min_ratio = 0.001,
  pln_num_penalties = 30,
  size_factors = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`ccs`     |     A Hooke cell_count_set object.
`main_model_formula_str`     |     A character string specifying the model of cell abundances across samples, where terms refer to columns in `colData(ccs)` . Put main effects here.
`nuisance_model_formula_str`     |     A character string specifying the model of cell abundances across samples. Put nuisance effects here.
`penalty_matrix`     |     A numeric NxN symmetric matrix specifying penalties for the PLN model, where N is the number of cell types. Entries must be positive. Use to specify an undirected graph prior for the PLN model.
`sparsity_factor`     |     A positive number to control how sparse the PLN network is. Larger values make the network more sparse.


## Value

a new cell_count_model object


# `new_cell_count_set`

Create a new cell_data_set object.


## Description

Create a new cell_data_set object.


## Usage

```r
new_cell_count_set(
  cds,
  sample_group,
  cell_group,
  sample_metadata = NULL,
  cell_metadata = NULL,
  lower_threshold = NULL,
  upper_threshold = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`cds`     |     A Monocle cell data set object.
`sample_group`     |     A column in colData(cds) that specifes how cells are grouped into samples
`cell_group`     |     A column in colData(cds) that specifies how cells are grouped into types or states (e.g. cluster)
`sample_metadata`     |     data frame containing attributes of individual samples, where
`cell_metadata`     |     data frame containing attributes of individual cell groups, where `row.names(cell_metadata)` are entries in `cell_group`


## Value

a new cell_data_set object


# `plot_contrast`

Plot a UMAP colored by how cells shift in a given contrast


## Description

Plot a UMAP colored by how cells shift in a given contrast


## Usage

```r
plot_contrast(
  ccm,
  cond_b_vs_a_tbl,
  log_abundance_thresh = -5,
  scale_shifts_by = c("receiver", "sender", "none"),
  edge_size = 2,
  cell_size = 1,
  q_value_thresh = 1,
  group_label_size = 2,
  plot_labels = c("significant", "all", "none"),
  fc_limits = c(-3, 3),
  sender_cell_groups = NULL,
  receiver_cell_groups = NULL,
  plot_edges = c("all", "directed", "undirected", "none"),
  label_cell_groups = list(),
  repel_labels = TRUE,
  model_for_pcors = "reduced",
  switch_label = NULL,
  sub_cds = NULL,
  alpha = 1,
  x = 1,
  y = 2
)
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model object.
`criterion`     |     a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".


## Value

an updated cell_count_model object


# `select_model`

Select the model a cell_count_model should use


## Description

Select the model a cell_count_model should use


## Usage

```r
select_model(
  ccm,
  criterion = "EBIC",
  sparsity_factor = 1,
  models_to_update = c("both", "full", "reduced")
)
```


## Arguments

Argument      |Description
------------- |----------------
`ccm`     |     A cell_count_model object.
`criterion`     |     a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".


## Value

an updated cell_count_model object


