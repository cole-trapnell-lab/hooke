Create a new cell\_count\_model object. — new\_cell\_count\_model • hooke

Toggle navigation [hooke](../index.html) 0.0.1

*   [Reference](../reference/index.html)
*   [Articles](#)
  *   [Hooke Tutorial](../articles/hooke_tutorial.html)
  
  Create a new cell\_count\_model object.
  =======================================
    
    `new_cell_count_model.Rd`
  
  Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a way of encoding a graph prior. Automatically selects sparsity parameter, but allows user to update it.
  
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
    vhat_method = c("bootstrap", "variational_var", "jackknife"),
    covariance_type = c("spherical", "diagonal"),
    num_bootstraps = 10,
    inception = NULL,
    backend = c("nlopt", "torch"),
    num_threads = 1,
    ftol_rel = 1e-06,
    penalize_by_distance = TRUE,
    penalty_scale_exponent = 2,
    reduction_method = "UMAP",
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
  
  whitelist
  
  list A data frame with two columns corresponding to (undirected) edges that should receive min\_penalty. The columns are integers that refer to cell clusters.
  
  blacklist
  
  list A data frame with two columns corresponding to (undirected) edges that should receive max\_penalty. The columns are integers that refer to cell clusters.
  
  sparsity\_factor
  
  A positive number to control how sparse the PLN network is. Larger values make the network more sparse. edges that should receive min\_penalty. The columns are either cell\_group names or integers that refer to cell\_groups in penalty\_matrix.
  
  base\_penalty
  
  numeric A factor that scales the penalty matrix.
  
  min\_penalty
  
  numeric A positive value that is assigned to whitelisted penalty matrix elements, which over-write existing values.
  
  max\_penalty
  
  numeric A positive value that is assigned to blacklisted penalty matrix elements. which over-write existing values.
  
  verbose
  
  logical Whether to emit verbose output.
  
  pseudocount
  
  integer A value added to the elements of the initial cell\_count\_set matrix.
  
  pln\_min\_ratio
  
  numeric Used in the definition of the sparsity penalty grid.
  
  pln\_num\_penalties
  
  integer Number of penalty values for the internally generated penalty grid.
  
  vhat\_method
  
  string Method used to compute covariance matrix?
    
    num\_bootstraps
  
  positive integer Number of iterations used with the bootstrap vhat\_method.
  
  inception
  
  Not used.
  
  backend
  
  Method used to run bootstrap iterations.
  
  Value
  -----
    
    a new cell\_count\_model object
  
  Contents
  --------
    
    Developed by Cole Trapnell.
  
  Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.