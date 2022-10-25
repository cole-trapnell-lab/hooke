test_that("new_cell_count_model works" ,{

  cds <- monocle3::load_a549()
  colData(cds)$sample = NULL

  cds = preprocess_cds(cds)
  cds = reduce_dimension(cds)
  cds = cluster_cells(cds, resolution=1e-2)
  #plot_cells(cds)

  ccs = new_cell_count_set(cds,
                           sample_group = "replicate",
                           cell_group = "cluster")

  ccm  = new_cell_count_model(ccs,
                              model_formula_str = "~log_dose")

  expect_is(ccm, "cell_count_model")

  # expect_error(cds <- new_cell_data_set(
  #   expression_data = as.data.frame(as.matrix(small_a549_exprs)),
  #   cell_metadata = small_a549_colData_df,
  #   gene_metadata = small_a549_rowData_df),
  #   "Argument expression_data must be a matrix - either sparse from the Matrix package or dense")
  #
  # expect_warning(cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                                         cell_metadata = small_a549_colData_df))
  #
  # expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                                       cell_metadata = small_a549_colData_df[1:100,],
  #                                       gene_metadata = small_a549_rowData_df),
  #              "cell_metadata must be NULL or have the same number of rows as columns in expression_data")
  #
  # expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                                       cell_metadata = small_a549_colData_df,
  #                                       gene_metadata = small_a549_rowData_df[1:100,]),
  #              "gene_metadata must be NULL or have the same number of rows as rows in expression_data")
  # temp <- small_a549_colData_df
  # row.names(temp)[1] <- "HP"
  # expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                                       cell_metadata = temp,
  #                                       gene_metadata = small_a549_rowData_df),
  #              "row.names of cell_metadata must be equal to colnames of expression_data")
  #
  # temp <- small_a549_rowData_df
  # row.names(temp)[1] <- "HP"
  # expect_error(cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                                       cell_metadata = small_a549_colData_df,
  #                                       gene_metadata = temp),
  #              "row.names of gene_metadata must be equal to row.names of expression_data")
  # cds <- new_cell_data_set(expression_data = small_a549_exprs,
  #                          cell_metadata = small_a549_colData_df,
  #                          gene_metadata = small_a549_rowData_df)
  # expect_is(cds, "cell_data_set")
  #
})

