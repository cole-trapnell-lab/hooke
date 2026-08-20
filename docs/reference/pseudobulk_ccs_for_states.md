Compute a pseudobulk expression matrix for a ccs

    pseudobulk_ccs_for_states(
      ccs,
      state_col = NULL,
      collapse_samples = FALSE,
      gene_group_df = NULL,
      agg_rowdata = NULL,
      cell_agg_fun = "mean",
      gene_agg_fun = "sum",
      norm_method = "size_only",
      scale_agg_values = FALSE,
      pseudocount = 0
    )

Arguments
---------

ccs

a cell count set object

state\_col

column to aggregate expression. Defaults using the cell\_group used in ccs construction.

collapse\_samples

boolean Whether to collapse sample groups into one.

gene\_group\_df

A dataframe in which the first column contains gene ids or short gene names and the second contains group membership for creating metagene. If NULL, genes are not grouped.

agg\_rowdata

A dataframe with metagene information, with the rows equal to the number of gene groupings. Must be specified if `gene_group_df` is not NULL. Otherwise this is ignored.
