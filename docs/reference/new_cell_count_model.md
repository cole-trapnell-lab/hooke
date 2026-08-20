Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a way of encoding a graph prior. Automatically selects sparsity parameter, but allows user to update it.
  
    new_cell_count_model(
      ccs,
      main_model_formula_str,
      nuisance_model_formula_str = "1",
      penalty_matrix = NULL,
      allowlist = NULL,
      denylist = NULL,
      sparsity_factor = 0.1,
      base_penalty = 1,
      min_penalty = 0.01,
      max_penalty = 1e+06,
      verbose = FALSE,
      pseudocount = 0,
      keep_ccs = TRUE,
      pln_min_ratio = 0.001,
      pln_num_penalties = 30,
      vhat_method = c("bootstrap", "variational_var", "sandwich_var", "jackknife"),
      covariance_type = c("spherical", "full", "diagonal"),
      num_bootstraps = 100,
      inception = NULL,
      backend = c("nlopt", "torch"),
      num_threads = 1,
      ftol_rel = 1e-06,
      penalize_by_distance = TRUE,
      penalty_scale_exponent = 2,
      reduction_method = "UMAP",
      random.seed = 42,
      ...
  )
  
  Arguments
  ---------
    
    ccs
  
  A Hooke cell\_count\_set object.
  
  main\_model\_formula\_str
  
  A character string specifying the model of cell abundances across samples, where terms refer to columns in`colData(ccs)`. Put main effects here.
  
  nuisance\_model\_formula\_str
  
  A character string specifying the model of cell abundances across samples. Put nuisance effects here.
  
  penalty\_matrix
  
  A numeric NxN symmetric matrix specifying penalties for the PLN model, where N is the number of cell types. Entries must be positive and the rows and columns must be named with the cell\_group names. Use to specify an undirected graph prior for the PLN model.
  
  allowlist
  
  data.frame Optional two-column data frame of undirected edges that should receive min\_penalty. The columns are either cell\_group names or integers that refer to cell\_groups in penalty\_matrix.
  
  denylist
  
  data.frame Optional two-column data frame of undirected edges that should receive max\_penalty. The columns are either cell\_group names or integers that refer to cell\_groups in penalty\_matrix.
  
  sparsity\_factor
  
  A positive number to control how sparse the PLN network is. Larger values make the network more sparse.
  
  base\_penalty
  
  numeric A factor that scales the penalty matrix.
  
  min\_penalty
  
  numeric A positive value that is assigned to penalty matrix elements in the allowlist, which over-write existing values.
  
  max\_penalty
  
  numeric A positive value that is assigned to penalty matrix elements in the denylist, which over-write existing values.
  
  verbose
  
  logical Whether to emit verbose output.
  
  pseudocount
  
  integer A value added to the elements of the initial cell\_count\_set matrix.
  
  keep\_ccs
  
  logical Whether to retain the original `ccs` object inside the fitted `cell_count_model`.
  
  pln\_min\_ratio
  
  numeric Used in the definition of the sparsity penalty grid.
  
  pln\_num\_penalties
  
  integer Number of penalty values for the internally generated penalty grid.
  
  vhat\_method
  
  string Method used to estimate the covariance matrix of fixed-effect coefficients. One of "bootstrap", "variational\_var", "sandwich\_var", or "jackknife".
  
  covariance\_type
  
  string Covariance structure for the full model. One of "spherical", "full", or "diagonal".
  
    num\_bootstraps
  
  positive integer Number of iterations used with the bootstrap vhat\_method.
  
  inception
  
  Not used.
  
  backend
  
  string Optimization backend passed to PLNmodels ("nlopt" or "torch").
  
  num\_threads
  
  positive integer Number of BLAS threads used while fitting.
  
  ftol\_rel
  
  numeric Relative tolerance passed to the optimizer.
  
  penalize\_by\_distance
  
  logical Whether to initialize penalties from between-centroid distances.
  
  penalty\_scale\_exponent
  
  numeric Exponent used when mapping distances to penalty weights.
  
  reduction\_method
  
  string Reduced-dimension embedding used for distance-based penalties (for example "UMAP").
  
  random.seed
  
  integer Random seed set before model fitting.
  
  Value
  -----
    
    a new cell\_count\_model object
  
