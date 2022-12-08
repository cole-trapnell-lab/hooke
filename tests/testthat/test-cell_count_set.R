context('test-cell_count_set')


test_that('new_cell_count_set works', {

  set.seed(2016)

  cds <- readRDS(system.file('extdata', 'pap_cds.trim_cell_type_sub3.2.rds', package = 'hooke'))
  cds <- detect_genes(cds)
  colData(cds)$cell_type <- colData(cds)$CW_assignedCellType
  colData(cds)$cluster <- monocle3::clusters(cds)

  colData(cds)$Size_Factor <- size_factors(cds)


  colData(cds)$experiment <- colData(cds)$sample
  colData(cds)$sample <- NULL

  expect_is(ccs <- new_cell_count_set(cds,
                                      sample_group = 'sampleName',
                                      cell_group = 'cell_type'), 'cell_data_set')

  expected_group_id <- c('19', '32', '59', '81', '22', '84', '27', '68')
  expected_sample <- c('15.pap.RBCKO.3_cds.RDS', '17.pap.GMKO.2_cds.RDS',
                       '22.pap.WTOLD.1_cds.RDS', '4.pap.Csf2.4_cds.RDS',
                       '15.pap.RBCKO.3_cds.RDS', '4.pap.Csf2.4_cds.RDS', 
                       '16.pap.GMKO.1_cds.RDS', '23.pap.WTOLD.2_cds.RDS')
  expected_cell_group <- c('hPAP alveolar macrophages', 'B cells',
                           'Healthy alveolar macrophages',
                           'hPAP alveolar macrophages',
                           'B cells', 'B cells', 'B cells', 'T cells')

  expect_equal(metadata(ccs)[['cell_group_assignments']][['group_id']][1:8], expected_group_id)
  expect_equal(metadata(ccs)[['cell_group_assignments']][['sample']][1:8], expected_sample)
  expect_equal(as.character(metadata(ccs)[['cell_group_assignments']][['cell_group']][1:8]), expected_cell_group)

  # sort(rowSums(counts(ccs)))
  # 153 207 230 422 493 495
  lower_threshold <- 207
  upper_threshold <- 422
  expect_is(ccs <- new_cell_count_set(cds,
                                      sample_group = 'sampleName',
                                      cell_group = 'cell_type',
                                      lower_threshold = lower_threshold,
                                      upper_threshold = upper_threshold),  'cell_data_set')

  expect_equal(as.vector(sort(Matrix::rowSums(counts(ccs)))), c(207, 230, 422))


  ccs <- new_cell_count_set(cds,
                            sample_group = 'sampleName',
                            cell_group = 'cell_type')

  # length(unique(metadata(ccs)$cell_group_assignments$sample))
  # [1] 20
  # length(unique(metadata(ccs)$cell_group_assignments$cell_group))
  # [1] 6

  # unique(colData(cds)[[cell_group]])   -> cell_metadata      rows of ccs count matrix
  # unique(colData(cds)[[sample_group]]) -> sample_metadata    columns of ccs count matrix
  cell_metadata <- data.frame(cell_metad_a = c('r_1', 'r_2', 'r_3', 'r_4', 'r_5', 'r_6'))
  rownames(cell_metadata) <- unique(metadata(ccs)$cell_group_assignments$cell_group)
  sample_metadata <- data.frame(sample = c('15.pap.RBCKO.3_cds.RDS', '17.pap.GMKO.2_cds.RDS', '22.pap.WTOLD.1_cds.RDS',
                                           '4.pap.Csf2.4_cds.RDS', '16.pap.GMKO.1_cds.RDS', '23.pap.WTOLD.2_cds.RDS',
                                           '14.pap.RBCKO.2_cds.RDS', '18.pap.WTYNG.3_cds.RDS', '6.pap.Csf2rb.4_cds.RDS',
                                           '7.pap.Csf2ra.4_cds.RDS', '24.pap.WTOLD.3_cds.RDS', '1.pap.Csf2ra.1_cds.RDS',
                                           '8.pap.WTOLD.4_cds.RDS', '21.pap.GMKO.3_cds.RDS', '3.pap.Csf2ra.2_cds.RDS', 
                                           '5.pap.Csf2ra.3_cds.RDS', '19.pap.WTYNG.1_cds.RDS', '20.pap.WTYNG.2_cds.RDS',
                                           '2.pap.WTYNG.4_cds.RDS', '13.pap.RBCKO.1_cds.RDS'),
                                sample_metad_a = c('c_1', 'c_2', 'c_3', 'c_4', 'c_5',
                                                   'c_6', 'c_7', 'c_8', 'c_9', 'c_10',
                                                   'c_11', 'c_12', 'c_13', 'c_14','c_15',
                                                   'c_16', 'c_17', 'c_18', 'c_19', 'c_20'))
  rownames(sample_metadata) <- unique(metadata(ccs)$cell_group_assignments$sample)

  expect_is(ccs <- new_cell_count_set(cds,
                                      sample_group = 'sampleName',
                                      cell_group = 'cell_type',
                                      sample_metadata = sample_metadata), 'cell_data_set')

  expect_is(ccs <- new_cell_count_set(cds,
                                      sample_group = 'sampleName',
                                      cell_group = 'cell_type',
                                      cell_metadata = cell_metadata), 'cell_data_set')

  lower_threshold <- 207
  upper_threshold <- 422
  expect_is(ccs <- new_cell_count_set(cds,
                                      sample_group = 'sampleName',
                                      cell_group = 'cell_type',
                                      sample_metadata = sample_metadata,
                                      cell_metadata = cell_metadata,
                                      lower_threshold = lower_threshold,
                                      upper_threshold = upper_threshold), 'cell_data_set')


})


test_that('new_cell_count_set problems', {

  
  small_a549_colData_df <- readRDS(system.file("extdata", "small_a549_dex_pdata.rda", package = "monocle3"))
  small_a549_rowData_df <- readRDS(system.file("extdata", "small_a549_dex_fdata.rda", package = "monocle3"))
  small_a549_exprs <- readRDS(system.file("extdata", "small_a549_dex_exprs.rda", package = "monocle3"))
  small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]
  small_a549_cds <- new_cell_data_set(small_a549_exprs,
                                      cell_metadata = small_a549_colData_df,
                                      gene_metadata = small_a549_rowData_df)

  sample_group_names_cds <- unique(colData(small_a549_cds)$sample)
  cell_group_names_cds <- unique(colData(small_a549_cds)$cell_type)

  # unique(colData(cds)[[cell_group]])   -> cell_metadata      rows of ccs count matrix
  # unique(colData(cds)[[sample_group]]) -> sample_metadata    columns of ccs count matrix
  expect_error(new_cell_count_set( cds = as.data.frame(as.matrix(small_a549_exprs)),
                                   sample_metadata = small_a549_colData_df,
                                   cell_metadata = small_a549_rowData_df),
               "Argument cds must be a cell_data_set")

  sample_group_test_1 <- data.frame(sample='sample_a', datum='big')
  rownames(sample_group_test_1) <- 'sample_a'
  sample_group_test_2 <- data.frame(sample=c('sample_a', 'sample_b'), datum=c('big', 'little'))
  rownames(sample_group_test_2) <- c('sample_a', 'sample_b')

  expect_error(new_cell_count_set(cds = small_a549_cds, sample_group='sample', cell_group='cell_type',
                                  sample_metadata = sample_group_test_1),
               'Argument sample_metadata must have sample group names in a column called \'sample\'.')

  expect_error(new_cell_count_set(cds = small_a549_cds, sample_group='sample', cell_group='cell_type',
                                  sample_metadata = sample_group_test_2),
               'Argument sample_metadata must have the same number of rows as there are distinct sample names in the cds column data.')


  cell_group_test_1 <- data.frame(cell_type='cell_a', datum='big')
  rownames(cell_group_test_1) <- 'cell_a'
  cell_group_test_1 <- data.frame(cell_type='A549', datum='big')
  cell_group_test_2 <- data.frame(cell_type=c('cell_a', 'cell_b'), datum=c('big', 'little'))
  rownames(cell_group_test_2) <- c('cell_a', 'cell_b')

  expect_error(new_cell_count_set(cds = small_a549_cds, sample_group='sample', cell_group='cell_type',
                                  cell_metadata = cell_group_test_1),
               'Argument cell_metadata row names must contain the cell_group names.')

  expect_error(new_cell_count_set(cds = small_a549_cds, sample_group='sample', cell_group='cell_type',
                                  cell_metadata = cell_group_test_2),
               'Argument cell_metadata must have the same number of rows as there are distinct cell_group names in the cds column data.')

})

