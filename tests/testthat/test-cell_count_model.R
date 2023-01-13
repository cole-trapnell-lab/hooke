context('test-cell_count_model')


test_that('new_cell_count_model works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

  # alt: /home/brent/work/hooke_data/macrophage/PAP/pap.al_2021-09-15.cds.RDS

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

  # Reduced cell_count_model
  model_list <- model(ccm, model_to_return='reduced')
  expect_equal(model_list$nb_param, 22)
  expect_equal(model_list$loglik, -405.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -438.0, tol=1.0e-1)
  expect_equal(model_list$ICL, -497.0, tol=1.0e-1)
  expect_equal(model_list$n_edges, 10)
  expect_equal(model_list$EBIC, -440.0, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -410.0, tol=1.0e-1)
  expect_equal(model_list$density, 0.556, tol=1.0e-3)

  expect_equal(model_list$latent[1,1], 1.858, tol=1.0e-1)
  expect_equal(attr(coef(model_list), 'vcov_variational')[1,1], 0.0106, tol=1.0e-3)


  # Full cell_count_model
  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 37)
  expect_equal(model_list$loglik, -384.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.6, tol=1.0e-1)
  expect_equal(model_list$ICL, -501.2, tol=1.0e-1)
  expect_equal(model_list$n_edges, 7)
  expect_equal(model_list$EBIC, -442.3, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -388.0, tol=1.0e-1)
  expect_equal(model_list$density, 0.389, tol=1.0e-3)

  expect_equal(model_list$latent[1,1], 1.859, tol=1.0e-3)
  expect_equal(attr(coef(model_list), 'vcov_variational')[1,1], 0.0178, tol=1.0e-3)

  # Penalty matrix.
  penalty_vector <- c(0.0020000, 0.2927265,  0.5065795,  0.3674560,  0.9472785, 1.0020000,
                      0.2927265, 0.00200000, 0.04575442, 0.05469007, 0.2746410, 0.4931481,
                      0.5065795, 0.04575442, 0.00200000, 0.02917695, 0.1012464, 0.2887256,
                      0.3674560, 0.05469007, 0.02917695, 0.00200000, 0.1374082, 0.2307490,
                      0.9472785, 0.2746410,  0.1012464 , 0.1374082,  0.0020000, 0.1076127,
                      1.0020000, 0.4931481,  0.2887256 , 0.2307490,  0.1076127, 0.0020000)
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names

  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -348.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -393.0, tol=1.0e-1)

  # Whitelist by cell_group names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  min_penalty <- 0.001
  whitelist <- data.frame(a=c('Aerocytes', 'Monocytes'), b=c('T cells', 'B cells'))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        whitelist = whitelist,
                                        min_penalty = min_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 31)
  expect_equal(model_list$loglik, -348.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -396.0, tol=1.0e-1)
 
  # Blacklist by cell_group names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  max_penalty <- 1.0e6
  blacklist <- data.frame(a=c('Aerocytes', 'Monocytes'), b=c('T cells', 'B cells'))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        blacklist = blacklist,
                                        max_penalty = max_penalty),
            'cell_count_model')

  # Whitelist by cell_group indices.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  min_penalty <- 0.001
  whitelist <- data.frame(a=c(1, 4), b=c(5, 6))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        whitelist = whitelist,
                                        min_penalty = min_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 31)
  expect_equal(model_list$loglik, -348.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -396.0, tol=1.0e-1)

  # Blacklist by cell_group indices.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  max_penalty <- 1.0e6
  blacklist <- data.frame(a=c(1, 4), b=c(5, 6))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        blacklist = blacklist,
                                        max_penalty = max_penalty),
            'cell_count_model')


  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -348.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -393.0, tol=1.0e-1)

  # Sparsity factor argument. 
  sparsity_factor <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        sparsity_factor = sparsity_factor),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 37)
  expect_equal(model_list$loglik, -384.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.6, tol=1.0e-1)

  # Base_penalty
  base_penalty <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        base_penalty = base_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 37)
  expect_equal(model_list$loglik, -384.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.4, tol=1.0e-1)

  # Pseudocount
  pseudocount <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        pseudocount = pseudocount),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 33)   
  expect_equal(model_list$loglik, -350.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -400.0, tol=1.0e-1)

  # Pln_min_ratio
  pln_min_ratio <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        pln_min_ratio = pln_min_ratio),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -396.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -441.0, tol=1.0e-1)


})


test_that('new_cell_count_model works', {

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

  # Reduced cell_count_model
  model_list <- model(ccm, model_to_return='reduced')
  expect_equal(model_list$nb_param, 22)
  expect_equal(model_list$loglik, -406.7, tol=1.0e-1)
  expect_equal(model_list$BIC, -438.2, tol=1.0e-1)
  expect_equal(model_list$ICL, -496.4, tol=1.0e-1)
  expect_equal(model_list$n_edges, 10)
  expect_equal(model_list$EBIC, -440.5, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -411.1, tol=1.0e-1)
  expect_equal(model_list$density, 0.5, tol=1.0e-1)

  expect_equal(model_list$latent[1,1], 1.9, tol=1.0e-1)
  expect_equal(attr(coef(model_list), 'vcov_variational')[1,1], 0.011, tol=1.0e-3)


  # Full cell_count_model
  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 37)
  expect_equal(model_list$loglik, -384.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.6, tol=1.0e-1)
  expect_equal(model_list$ICL, -501.2, tol=1.0e-1)
  expect_equal(model_list$n_edges, 7)
  expect_equal(model_list$EBIC, -442.3, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -389.0, tol=1.0e-1)
  expect_equal(model_list$density, 0.389, tol=1.0e-3)

  expect_equal(model_list$latent[1,1], 1.859, tol=1.0e-3)
  expect_equal(attr(coef(model_list), 'vcov_variational')[1,1], 0.0178, tol=1.0e-3)


  # Penalty matrix.
  penalty_vector <- c(0.0020000, 0.2927265,  0.5065795,  0.3674560,  0.9472785, 1.0020000,
                      0.2927265, 0.00200000, 0.04575442, 0.05469007, 0.2746410, 0.4931481,
                      0.5065795, 0.04575442, 0.00200000, 0.02917695, 0.1012464, 0.2887256,
                      0.3674560, 0.05469007, 0.02917695, 0.00200000, 0.1374082, 0.2307490,
                      0.9472785, 0.2746410,  0.1012464 , 0.1374082,  0.0020000, 0.1076127,
                      1.0020000, 0.4931481,  0.2887256 , 0.2307490,  0.1076127, 0.0020000)
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names

  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -349.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -394.0, tol=1.0e-1)

  # Whitelist by cell_group names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  min_penalty <- 0.001
  whitelist <- data.frame(a=c('Aerocytes', 'Monocytes'), b=c('T cells', 'B cells'))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        whitelist = whitelist,
                                        min_penalty = min_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 32)
  expect_equal(model_list$loglik, -385.8, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.7, tol=1.0e-1)
 
  # Blacklist by cell_group names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  max_penalty <- 1.0e6
  blacklist <- data.frame(a=c('Aerocytes', 'Monocytes'), b=c('T cells', 'B cells'))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        blacklist = blacklist,
                                        max_penalty = max_penalty),
            'cell_count_model')

  # Whitelist by cell_group indices.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  min_penalty <- 0.001
  whitelist <- data.frame(a=c(1, 4), b=c(5, 6))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        whitelist = whitelist,
                                        min_penalty = min_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 32)
  expect_equal(model_list$loglik, -385.8, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.7, tol=1.0e-1)

  # Blacklist by cell_group indices.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- cell_group_names
  max_penalty <- 1.0e6
  blacklist <- data.frame(a=c(1, 4), b=c(5, 6))
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        penalty_matrix = penalty_matrix,
                                        blacklist = blacklist,
                                        max_penalty = max_penalty),
            'cell_count_model')


  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -349.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -394.0, tol=1.0e-1)

  # Sparsity factor argument. 
  sparsity_factor <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        sparsity_factor = sparsity_factor),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 37)
  expect_equal(model_list$loglik, -384.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.6, tol=1.0e-1)

  # Base_penalty
  base_penalty <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        base_penalty = base_penalty),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 38)
  expect_equal(model_list$loglik, -384.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -439.4, tol=1.0e-1)

  # Pseudocount
  pseudocount <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        pseudocount = pseudocount),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 32)
  expect_equal(model_list$loglik, -382.5, tol=1.0e-1)
  expect_equal(model_list$BIC, -397.0, tol=1.0e-1)

  # Pln_min_ratio
  pln_min_ratio <- 0.1
  expect_is(ccm <- new_cell_count_model(ccs,
                                        main_model_formula_str = "Genotype",
                                        nuisance_model_formula_str = "batch",
                                        pln_min_ratio = pln_min_ratio),
            'cell_count_model')

  model_list <- model(ccm, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -396.0, tol=1.0e-1)
  expect_equal(model_list$BIC, -441.0, tol=1.0e-1)


})


test_that('new_cell_count_model problems', {

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

  #
  # Checks
  #   o  main_model_formula_str
  #   o  nuisance_model_formula_str
  #   o  penalty_matrix: no row or column names; dimensions (other?)
  #   o  whitelist/blacklist: dimensions; row/col names
  #   o  sparsity_factor: positive factor
  #   o  min_penalty and max_penalty
  #   o  norm_method
  #   o  vhat_method
  #   o  backend

  # main_model_formula_str
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Phenotype",
                                    nuisance_model_formula_str = "batch" ),
               '^object \'Phenotype\' not found$')

  # main_model_formula_str 
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "basket" ),
               '^object \'basket\' not found$')

  # Penalty matrix.
  penalty_vector <- c(0.0020000, 0.2927265,  0.5065795,  0.3674560,  0.9472785, 1.0020000,
                      0.2927265, 0.00200000, 0.04575442, 0.05469007, 0.2746410, 0.4931481,
                      0.5065795, 0.04575442, 0.00200000, 0.02917695, 0.1012464, 0.2887256,
                      0.3674560, 0.05469007, 0.02917695, 0.00200000, 0.1374082, 0.2307490,
                      0.9472785, 0.2746410,  0.1012464 , 0.1374082,  0.0020000, 0.1076127,
                      1.0020000, 0.4931481,  0.2887256 , 0.2307490,  0.1076127, 0.0020000)
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  cell_group_names <- c('Aerocytes', 'Healthy alveolar macrophages', 'hPAP alveolar macrophages', 'Monocytes', 'T cells', 'B cells')

  # No row and no column names.
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    penalty_matrix = penalty_matrix))

  # Bad row names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  rownames(penalty_matrix) <- c('1', '2', '3', '4', '5', '6')
  colnames(penalty_matrix) <- cell_group_names

  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch", 
                                    penalty_matrix = penalty_matrix))

  # Bad column names.
  penalty_matrix <- matrix(penalty_vector, nrow=6)
  rownames(penalty_matrix) <- cell_group_names
  colnames(penalty_matrix) <- c('1', '2', '3', '4', '5', '6')

  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    penalty_matrix = penalty_matrix))

  # Bad whitelist.
  whitelist <- data.frame(a=c('Aerocytes', 'parsley'), b=c('Monocytes', 'T cells'))
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    whitelist = whitelist))

  whitelist <- data.frame(a=c(1,8), b=c(2,3))
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    whitelist = whitelist))

  # Bad blacklist.
  blacklist <- data.frame(a=c('Aerocytes', 'parsley'), b=c('Monocytes', 'T cells'))
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    blacklist = blacklist))

  blacklist <- data.frame(a=c(1,8), b=c(2,3))
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    blacklist = blacklist))

  sparsity_factor <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    sparsity_factor = sparsity_factor))

  base_penalty <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    base_penalty = base_penalty))

  min_penalty <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    min_penalty = min_penalty))

  max_penalty <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    max_penalty = max_penalty))

  pseudocount <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    pseudocount = pseudocount))

  pln_min_ratio <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    pln_min_ratio = pln_min_ratio))

  pln_num_penalties <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    pln_num_penalties = pln_num_penalties))

  norm_method <- 'turnips'
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    norm_method = norm_method))

  vhat_method <- 'rutabagas'
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    vhat_method = vhat_method))

  size_factors <- c('jackfruit', 'parsnips')
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    size_factors = size_factors))

  num_bootstraps <- -1
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    num_bootstraps = num_bootstraps))

  backend <- 'rocquet'
  expect_error(new_cell_count_model(ccs,
                                    main_model_formula_str = "Genotype",
                                    nuisance_model_formula_str = "batch",
                                    backend = backend))
})


test_that('select_model works', {

  # On Github Actions
  #  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
  ccs <- new_cell_count_set(cds, sample_group = 'sampleName',
                            cell_group = 'cell_type')
  ccm <- new_cell_count_model(ccs, main_model_formula_str = "Genotype",
                              nuisance_model_formula_str = "batch")

  ccm2 <- select_model(ccm, criterion='BIC', sparsity_factor=2.0, models_to_update='both')


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='BIC',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)
  expect_equal(model_list$loglik, -413.7, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.2, tol=1.0e-1)
  expect_equal(model_list$ICL, -494.9, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -438.6, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -464.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -390.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$ICL, -497.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -448.8, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='EBIC',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)
  expect_equal(model_list$loglik, -413.7, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.2, tol=1.0e-1)
  expect_equal(model_list$ICL, -494.9, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -438.6, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -464.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -390.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$ICL, -497.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -448.8, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='StARS',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)    
  expect_equal(model_list$loglik, -414.9, tol=1.0e-1)
  expect_equal(model_list$BIC, -437.3, tol=1.0e-1)
  expect_equal(model_list$ICL, -496.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -439.8, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -474.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -391.1, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.0, tol=1.0e-1)
  expect_equal(model_list$ICL, -498.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -436.0, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -460.0, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)

})


test_that('select_model works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

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
  ccs <- new_cell_count_set(cds, sample_group = 'sampleName',
                            cell_group = 'cell_type')
  ccm <- new_cell_count_model(ccs, main_model_formula_str = "Genotype",
                              nuisance_model_formula_str = "batch")

  ccm2 <- select_model(ccm, criterion='BIC', sparsity_factor=2.0, models_to_update='both')


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='BIC',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)
  expect_equal(model_list$loglik, -413.7, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.2, tol=1.0e-1)
  expect_equal(model_list$ICL, -494.9, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -438.6, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -464.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -390.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$ICL, -497.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -448.8, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='EBIC',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)
  expect_equal(model_list$loglik, -413.7, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.2, tol=1.0e-1)
  expect_equal(model_list$ICL, -494.9, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -438.6, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -464.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -390.2, tol=1.0e-1)
  expect_equal(model_list$BIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$ICL, -497.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -435.1, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -448.8, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)


  expect_is(ccm2 <- select_model(ccm,
                                 criterion='StARS',
                                 sparsity_factor=2.0,
                                 models_to_update='both'),
            'cell_count_model')

  model_list <- model(ccm2, model_to_return='reduced')
  expect_equal(model_list$nb_param, 15)    
  expect_equal(model_list$loglik, -414.9, tol=1.0e-1)
  expect_equal(model_list$BIC, -437.3, tol=1.0e-1)
  expect_equal(model_list$ICL, -496.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 3)
  expect_equal(model_list$EBIC, -439.8, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -474.2, tol=1.0e-1)
  expect_equal(model_list$density, .167, tol=1.0e-3)

  model_list <- model(ccm2, model_to_return='full')
  expect_equal(model_list$nb_param, 30)
  expect_equal(model_list$loglik, -391.1, tol=1.0e-1)
  expect_equal(model_list$BIC, -436.0, tol=1.0e-1)
  expect_equal(model_list$ICL, -498.1, tol=1.0e-1)
  expect_equal(model_list$n_edges, 0)
  expect_equal(model_list$EBIC, -436.0, tol=1.0e-1)
  expect_equal(model_list$pen_loglik, -460.0, tol=1.0e-1)
  expect_equal(model_list$density, .0, tol=1.0e-1)

})

