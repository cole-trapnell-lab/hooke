Function reference • hooke

Toggle navigation [hooke`](../index) 0.0.1

*   [Reference`](../index)
*   [Articles`](#)
    *   [Hooke Tutorial`](../articles/hooke_tutorial)

Reference
=========

All functions[`](#all-functions)
-------------------------------

[`add_covariate()`](new_cell_count_set)

[`add_covariate()`](add_covariate)

add metadata to pb\_cds from cds

[`build_interval_formula()](build_interval_formula)

Builds a model formula for time series models based on the range of the data. This is a utility function that puts the knots in reasonable positions based on the range of the data.


[`cds()](cds)

Get the underlying cell\_data\_set object from a cell\_count\_model.

[`cell_count_model`](cell_count_model) [`cell_count_model-class`](cell_count_model)

The cell\_count\_model class

[`cell_count_set`](cell_count_set) [`cell_count_set-class`](cell_count_set)

The cell\_count\_set class

[`centroids()](centroids)

Get the centroids of cell groups in UMAP/PCA space.

[`compare_abundances()](compare_abundances)

Compare two estimates of cell abundances from a Hooke model.

[`compare_ko_to_wt_at_timepoint()](compare_ko_to_wt_at_timepoint)

Helper function to plot kinetics

[`estimate_abundances()](estimate_abundances)

Predict cell type abundances given a PLN model and a set of inputs for its covariates

[`estimate_abundances_over_interval()](estimate_abundances_over_interval)

Predict cell type abundances given a PLN model over a range of time or other interval

[`get_norm_counts()](get_norm_counts)

return a matrix of normalized counts

[`get_paga_graph()](get_paga_graph)

Extract a partitioned abstract graph from a Monocle cell\_data\_set object

[`hooke_theme_opts()](hooke_theme_opts)

Default plotting options for ggplot2

[`initial_pcor_graph()](initial_pcor_graph)

Get an initial graph for use as a whitelist in fitting a cell count model

[`model()](model)

Get the selected PLNnetwork from a cell\_count\_model object.

[`new_cell_count_model()](new_cell_count_model)

Create a new cell\_count\_model object.

[`new_cell_count_set()](new_cell_count_set)

Create a new cell\_data\_set object.

[`plot_contrast()](plot_contrast)

Plot a UMAP colored by how cells shift in a given contrast

[`pseudobulk_ccs_for_states()](pseudobulk_ccs_for_states)

Compute a pseudobulk expression matrix for a ccs

[`select_model()](select_model)

Select the model a cell\_count\_model should use

[`subset_ccs()](subset_ccs)

subset ccs by cell groups

Contents
--------

Developed by Cole Trapnell.
