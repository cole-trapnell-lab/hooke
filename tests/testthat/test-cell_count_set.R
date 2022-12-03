context('test-cell_count_set')

small_a549_colData_df <- readRDS(system.file("extdata", "small_a549_dex_pdata.rda", package = "monocle3"))
small_a549_rowData_df <- readRDS(system.file("extdata", "small_a549_dex_fdata.rda", package = "monocle3"))
small_a549_exprs <- readRDS(system.file("extdata", "small_a549_dex_exprs.rda", package = "monocle3"))
small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]


test_that("new_cell_count_set works" ,{
  set.seed(2016)

  cds <- monocle3::load_a549()
  colData(cds)$sample <- NULL

  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds, resolution=1e-2)
  colData(cds)$cluster <- clusters(cds)

  #plot_cells(cds)

  ccs <- new_cell_count_set(cds,
                           sample_group = "replicate",
                           cell_group = "cluster")

  expect_is(ccs, "cell_count_set")
  expect_equal(ccs@metadata$cell_group_assignments$group_id[1:4], c('7', '123', '50', '74'))
  expect_equal(ccs@metadata$cell_group_assignments$sample[1:4], c('AC03', 'BD08', 'AD06', 'BC02'))
  expect_equal(as.integer(ccs@metadata$cell_group_assignments$cell_group[1:4]), c(1, 1, 1, 2))


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

