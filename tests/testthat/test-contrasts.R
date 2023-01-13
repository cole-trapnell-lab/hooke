context('test-contrasts')


test_that('estimate_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  abund <- estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))

  expect_equal(abund[['log_abund']][[1]], 2.61, tol=1.0e-2)
  expect_equal(abund[['log_abund_se']][[1]], 0.0648, tol=1.0e-3)
  expect_equal(abund[['log_abund_sd']][[1]], 1.38, tol=1.0e-2)

})


test_that('estimate_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  abund <- estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))

  expect_equal(abund[['log_abund']][[1]], 2.85, tol=1.0e-2)
  expect_equal(abund[['log_abund_se']][[1]], 0.0648, tol=1.0e-3)
  expect_equal(abund[['log_abund_sd']][[1]], 1.38, tol=1.0e-2)

})


test_that('estimate_abundances problems', {

  # On Github Actions
  # skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  ccm2 <- 5
  expect_error(estimate_abundances(ccm2, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1)))
  newdata <- data.frame(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1)
  expect_error(estimate_abundances(ccm, newdata))
  newdata <- tibble::tibble(Genotype="Csf5ra-/-", Age="12-13 weeks", batch="A", experiment=1)
  expect_error(estimate_abundances(ccm, newdata))
  newdata <- tibble::tibble(Phenotype="odd", Age="12-13 weeks", batch="A", experiment=1)
  expect_error(estimate_abundances(ccm, newdata))

})


test_that('compare_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
  cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
  comp_abund = compare_abundances(ccm, cond_wt, cond_csf2ra)

  expect_equal(comp_abund[['log_abund_x']][1], 2.5, tol=1.0e-1)
  expect_equal(comp_abund[['log_abund_se_x']][1], 0.041, tol=1.0e-3)
  expect_equal(comp_abund[['log_abund_sd_x']][1], 1.38, tol=1.0e-2)

  expect_equal(comp_abund[['log_abund_y']][1], 2.61, tol=1.0e-2) # 69
  expect_equal(comp_abund[['log_abund_se_y']][1], 0.065, tol=1.0e-3)
  expect_equal(comp_abund[['log_abund_sd_y']][1], 1.377, tol=1.0e-1)

  expect_equal(comp_abund[['delta_log_abund']][1], 0.192, tol=1.0e-3)
  expect_equal(comp_abund[['delta_p_value']][1], 0.025, tol=1.0e-3)
  expect_equal(comp_abund[['delta_q_value']][1], 0.025, tol=1.0e-3)

})


test_that('compare_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
  cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
  comp_abund = compare_abundances(ccm, cond_wt, cond_csf2ra)

  expect_equal(comp_abund[['log_abund_x']][1], 2.5, tol=1.0e-1)
  expect_equal(comp_abund[['log_abund_se_x']][1], 0.041, tol=1.0e-3)
  expect_equal(comp_abund[['log_abund_sd_x']][1], 1.377, tol=1.0e-3)

  expect_equal(comp_abund[['log_abund_y']][1], 2.852, tol=1.0e-3)
  expect_equal(comp_abund[['log_abund_se_y']][1], 0.065, tol=1.0e-3)
  expect_equal(comp_abund[['log_abund_sd_y']][1], 1.377, tol=1.0e-1)

  expect_equal(comp_abund[['delta_log_abund']][1], 0.387, tol=1.0e-3)
  expect_equal(comp_abund[['delta_p_value']][1], 1.756e-4, tol=1.0e-3)
  expect_equal(comp_abund[['delta_q_value']][1], 2.107e-4, tol=1.0e-3)

})


test_that('compare_abundances problems', {

  # Not on Github Actions
  # skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
                                      cell_group = 'cell_type'),
            'cell_data_set')
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch"),
            'cell_count_model')
  ccm <- select_model(ccm, criterion='StARS', sparsity_factor=6.0, models_to_update='both')

  cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
  cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))

  ccm2 <- 5
  expect_error(compare_abundances(ccm2, cond_wt, cond_csf2ra))
  cond_wt2 <- as.data.frame(cond_wt)
  expect_error(compare_abundances(ccm, cond_wt2, cond_csf2ra))
  cond_csf2ra2 <- as.data.frame(cond_wt)
  expect_error(compare_abundances(ccm, cond_wt, cond_csf2ra2))
  method <- 'my_method'
  expect_error(compare_abundances(ccm, cond_wt, cond_csf2ra, method=method))


})

