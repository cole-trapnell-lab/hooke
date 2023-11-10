Function reference • hooke

Toggle navigation [hooke](reference/../index.md) 0.0.1

*   [Reference](reference/../reference/index.md)
*   [Articles](reference/#)
    *   [Hooke Tutorial](reference/../articles/hooke_tutorial.md)

Reference
=========

All functions[](reference/#all-functions)
-------------------------------

`[add_covariate()](reference/add_covariate.md)`

add metadata to pb\_cds from cds

`[build_interval_formula()](reference/build_interval_formula.md)`

Builds a model formula for time series models based on the range of the data. This is a utility function that puts the knots in reasonable positions based on the range of the data.


`[cds()](reference/cds.md)`

Get the underlying cell\_data\_set object from a cell\_count\_model.

`[cell_count_model](reference/cell_count_model.md)` `[cell_count_model-class](reference/cell_count_model.md)`

The cell\_count\_model class

`[cell_count_set](reference/cell_count_set.md)` `[cell_count_set-class](reference/cell_count_set.md)`

The cell\_count\_set class

`[centroids()](reference/centroids.md)`

Get the centroids of cell groups in UMAP/PCA space.

`[compare_abundances()](reference/compare_abundances.md)`

Compare two estimates of cell abundances from a Hooke model.

`[compare_ko_to_wt_at_timepoint()](reference/compare_ko_to_wt_at_timepoint.md)`

Helper function to plot kinetics

`[estimate_abundances()](reference/estimate_abundances.md)`

Predict cell type abundances given a PLN model and a set of inputs for its covariates

`[estimate_abundances_over_interval()](reference/estimate_abundances_over_interval.md)`

Predict cell type abundances given a PLN model over a range of time or other interval

`[get_norm_counts()](reference/get_norm_counts.md)`

return a matrix of normalized counts

`[get_paga_graph()](reference/get_paga_graph.md)`

Extract a partitioned abstract graph from a Monocle cell\_data\_set object

`[hooke_theme_opts()](reference/hooke_theme_opts.md)`

Default plotting options for ggplot2

`[initial_pcor_graph()](reference/initial_pcor_graph.md)`

Get an initial graph for use as a whitelist in fitting a cell count model

`[model()](reference/model.md)`

Get the selected PLNnetwork from a cell\_count\_model object.

`[new_cell_count_model()](reference/new_cell_count_model.md)`

Create a new cell\_count\_model object.

`[new_cell_count_set()](reference/new_cell_count_set.md)`

Create a new cell\_data\_set object.

`[plot_contrast()](reference/plot_contrast.md)`

Plot a UMAP colored by how cells shift in a given contrast

`[pseudobulk_ccs_for_states()](reference/pseudobulk_ccs_for_states.md)`

Compute a pseudobulk expression matrix for a ccs

`[select_model()](reference/select_model.md)`

Select the model a cell\_count\_model should use

`[subset_ccs()](reference/subset_ccs.md)`

subset ccs by cell groups

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](reference/https://pkgdown.r-lib.org/) 2.0.7.