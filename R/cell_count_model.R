setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))
setOldClass(c("PLNnetworkfamily"), prototype=structure(list(), class="PLNnetworkfamily"))
setOldClass(c("PLNnetworkfit"), prototype=structure(list(), class="PLNnetworkfit"))


#' The cell_count_set class
#'
#' The main class used by Hooke to hold cell abundances data.
#' cell_count_set extends the Monocle's cell_data_set class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field cds cell_data_set, the Monocle cell data set object that this class models.
#' @name cell_count_set
#' @rdname cell_count_set
#' @aliases cell_count_set-class
#' @exportClass cell_count_set
setClass("cell_count_set",
         contains = "cell_data_set",
         slots = c(cds = "cell_data_set",
                   info = "SimpleList")
)


#' The cell_count_model class
#'
#' The main class used by Monocle3 to hold single-cell expression data.
#' cell_count_model extends the Bioconductor SingleCellExperiment class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field ccs cell_count_set the underlying object of cell counts this class models.
#' @field model_formula specifies the model for the cell abundances.
#' @field model_family PLNnetworkfamily a family of PLNnetwork models for the cell count data.
#' @field best_full_model PLNnetworkfit the current best PLN network model for the cell count data. Both main effects and nuisance terms
#' @field best_reduced_model PLNnetworkfit the current best (reduced) PLN network model for the cell count data. Nuisance terms only
#' @field model_aux SimpleList auxiliary information from from PLN model construction
#' @name cell_count_model
#' @rdname cell_count_model
#' @aliases cell_count_model-class
#' @import PLNmodels
#' @exportClass cell_count_model
setClass("cell_count_model",
         slots = c(ccs = "cell_count_set",
                   full_model_formula = "formula",
                   full_model_family = "PLNnetworkfamily",
                   best_full_model = "PLNnetworkfit",
                   reduced_model_formula = "formula",
                   reduced_model_family = "PLNnetworkfamily",
                   best_reduced_model = "PLNnetworkfit",
                   sparsity = "numeric",
                   model_aux = "SimpleList",
                   bootstrapped_vhat = "matrix")
)





#' Create a new cell_data_set object.
#'
#' @param cds A Monocle cell data set object.
#' @param sample_group A column in colData(cds) that specifes how cells are grouped into samples
#' @param cell_group A column in colData(cds) that specifies how cells are grouped into types or states (e.g. cluster)
#' @param sample_metadata data frame containing attributes of individual samples, where
#' @param cell_metadata data frame containing attributes of individual cell groups, where
#'   \code{row.names(cell_metadata)} are entries in \code{cell_group}
#' @return a new cell_data_set object
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
#' @export
new_cell_count_set <- function(cds,
                               sample_group,
                               cell_group,
                               sample_metadata = NULL,
                               cell_metadata = NULL,
                               lower_threshold = NULL,
                               upper_threshold = NULL) {

  colData(cds)$sample = NULL

  # check if anything contains NAs in it
  # if so drop them
  num_sample_group_NAs = sum(is.na(colData(cds)[[sample_group]]))
  if (num_sample_group_NAs != 0) {
    message(paste(num_sample_group_NAs, "NAs found in sample group. Dropping NAs."))
    cds = cds[, !is.na(colData(cds)[[sample_group]])]
  }

  num_cell_group_NAs = sum(is.na(colData(cds)[[cell_group]]))
  if (num_cell_group_NAs != 0) {
    message(paste(num_cell_group_NAs, "NAs found in cell group. Dropping NAs."))
    cds = cds[, !is.na(colData(cds)[[cell_group]])]
  }



  coldata_df = colData(cds) %>% tibble::as_tibble()
  # current commented out bc mess w projection clusters
  # coldata_df$cluster = monocle3::clusters(cds)
  # coldata_df$partition = partitions(cds)

  coldata_df = coldata_df %>% dplyr::rename_("sample" = sample_group, "cell_group" = as.character(cell_group))
  #coldata_df$cell_group = factor(coldata_df$cell_group, levels=unique(colData(cds)[,cell_group]))

  coldata_df$group_id = coldata_df %>%
    dplyr::group_indices_("sample", "cell_group") %>% as.character

  # add to cds
  colData(cds)$group_id = coldata_df$group_id

  cds_summary = coldata_df %>%
    dplyr::group_by_("sample", "cell_group") %>%
    dplyr::summarize(cells = dplyr::n())

  cds_covariates_df = coldata_df %>%
    dplyr::select(-cell_group) %>%
    dplyr::group_by_("sample") %>%
    dplyr::summarize(across(where(is.numeric), mean),
                     across(where(is.factor), function(x) { tail(names(sort(table(x))), 1) }),
                     across(where(is.character), function(x) { tail(names(sort(table(x, useNA="ifany"))), 1) } ))

  if (is.null(sample_metadata) == FALSE){
    cds_covariates_df = left_join(cds_covariates_df, sample_metadata, by=c("sample"="sample"))
  }

  cds_covariates_df = cds_covariates_df %>% as.data.frame(cds_covariates_df, stringsAsFactors=FALSE)
  row.names(cds_covariates_df) = cds_covariates_df %>% dplyr::pull(sample)

  cell_counts_wide = tidyr::spread(cds_summary, sample, cells, fill=0)
  cell_states = as.character(cell_counts_wide %>% dplyr::pull("cell_group"))
  cell_counts_wide = as.matrix(cell_counts_wide[,2:ncol(cell_counts_wide)])

  row.names(cell_counts_wide) = cell_states

  # filter out cell groups based on counts
  if (is.null(lower_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) >= lower_threshold, ]
  }
  if (is.null(upper_threshold) == FALSE) {
    cell_counts_wide = cell_counts_wide[Matrix::rowSums(cell_counts_wide) <= upper_threshold, ]
  }

  # remove from cds
  removed_cell_states = setdiff(cell_states, rownames(cell_counts_wide))

  #cell_counts_wide = t(cell_counts_wide)

  cds_covariates_df = cds_covariates_df[colnames(cell_counts_wide),]

  # This is super confusing because of the way the arguments are named in new_cell_data_set.
  # We are making a matrix of dimension MxN, where M are cell types and N are samples (e.g. embryos, replicates, etc).
  # The "gene" metadata monocle normally expects will actually be used to hold
  ccs = methods::new("cell_count_set",
               monocle3::new_cell_data_set(cell_counts_wide,
                                           cell_metadata=cds_covariates_df,
                                           gene_metadata=cell_metadata),
               cds=cds[, !colData(cds)[[cell_group]] %in% removed_cell_states],
               info=SimpleList(sample_group=sample_group,
                               cell_group=cell_group))


  # assertthat::assert_that(class(expression_data) == "matrix" ||
  #                           is_sparse_matrix(expression_data),
  #                         msg = paste("Argument expression_data must be a",
  #                                     "matrix - either sparse from the",
  #                                     "Matrix package or dense"))
  # if (!is.null(cell_metadata)) {
  #   assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data),
  #                           msg = paste("cell_metadata must be NULL or have",
  #                                       "the same number of rows as columns",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(cell_metadata)) &
  #                             all(row.names(cell_metadata) == colnames(expression_data)),
  #                           msg = paste("row.names of cell_metadata must be equal to colnames of",
  #                                       "expression_data"))
  # }
  #
  # if (!is.null(gene_metadata)) {
  #   assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data),
  #                           msg = paste("gene_metadata must be NULL or have",
  #                                       "the same number of rows as rows",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(
  #     row.names(gene_metadata) == row.names(expression_data)),
  #     msg = paste("row.names of gene_metadata must be equal to row.names of",
  #                 "expression_data"))
  # }
  #
  # if (is.null(cell_metadata)) {
  #   cell_metadata <- data.frame(cell = colnames(expression_data),
  #                               row.names = colnames(expression_data))
  # }
  #
  # if(!('gene_short_name' %in% colnames(gene_metadata))) {
  #   warning(paste("Warning: gene_metadata must contain a column verbatim",
  #                 "named 'gene_short_name' for certain functions."))
  # }
  #
  # sce <- SingleCellExperiment(list(counts=methods::as(expression_data, "dgCMatrix")),
  #                             rowData = gene_metadata,
  #                             colData = cell_metadata)
  #
  # cds <- methods::new("cell_data_set",
  #                     assays = SummarizedExperiment::Assays(
  #                       list(counts=methods::as(expression_data, "dgCMatrix"))),
  #                     colData = colData(sce),
  #                     int_elementMetadata =int_elementMetadata(sce),
  #                     int_colData = int_colData(sce),
  #                     int_metadata = int_metadata(sce),
  #                     metadata = metadata(sce),
  #                     NAMES = NULL,
  #                     elementMetadata = elementMetadata(sce)[,0],
  #                     rowRanges = rowRanges(sce))
  #
  # metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  # clusters <- stats::setNames(SimpleList(), character(0))
  # cds <- estimate_size_factors(cds)
  # cds



  ccs@metadata[["cell_group_assignments"]] = coldata_df %>% dplyr::select(group_id, sample, cell_group) %>% as.data.frame
  row.names(ccs@metadata[["cell_group_assignments"]]) = colnames(cds)
  ccs@metadata[["cell_group_assignments"]] = ccs@metadata[["cell_group_assignments"]] %>% filter(!cell_group %in% removed_cell_states)

  return (ccs)
}

#' resamples the ccs counts using a multinomial distribution
#' @param ccs
#' @param random.seed
bootstrap_ccs = function(ccs, random.seed=NULL) {
  count_mat = counts(ccs)
  num_cols = dim(count_mat)[2]
  set.seed(random.seed)
  count_list = lapply(seq(1,num_cols), function(i){
    counts = count_mat[,i]
    size = sum(count_mat[,i])
    prob = counts/size
    sample_counts = rmultinom(1, size, prob)

  })
  new_count_mat = do.call(cbind,count_list)
  colnames(new_count_mat) = colnames(count_mat)
  counts(ccs) = new_count_mat
  return(ccs)
}

#'
#' @param ccs
#' @param full_model_formula_str
#' @param best_full_model
#' @param reduced_pln_model
#' @param pseudocount
#' @param initial_penalties
#' @param pln_min_ratio
#' @param pln_num_penalties
#' @param random.seed
bootstrap_model = function(ccs,
                           full_model_formula_str,
                           best_full_model,
                           reduced_pln_model,
                           pseudocount,
                           initial_penalties,
                           pln_min_ratio,
                           pln_num_penalties,
                           random.seed) {

  # resample the counts
  sub_ccs = bootstrap_ccs(ccs, random.seed = random.seed)

  # remake data from new ccs
  sub_pln_data <- PLNmodels::prepare_data(counts = counts(sub_ccs) + pseudocount,
                                          covariates = colData(sub_ccs) %>% as.data.frame,
                                          offset = monocle3::size_factors(sub_ccs))
  # rerun the model using the same initial parameters
  # as the original, non bootstrapped model
  sub_full_model = do.call(PLNmodels::PLNnetwork, args=list(full_model_formula_str,
                                                            data = sub_pln_data,
                                                            penalties = reduced_pln_model$penalties,
                                                            control_init=list(min.ratio=pln_min_ratio, nPenalties=pln_num_penalties),
                                                            control_main=list(penalty_weights=initial_penalties,
                                                                              trace = ifelse(FALSE, 2, 0),
                                                                              inception = best_full_model)))

  return(sub_full_model)

}



#' computes the avg vhat across n bootstraps
#' @param ccm
#' @param num_bootstraps
bootstrap_vhat = function(ccs,
                          full_model_formula_str,
                          best_full_model,
                          best_reduced_model,
                          reduced_pln_model,
                          pseudocount,
                          initial_penalties,
                          pln_min_ratio,
                          pln_num_penalties,
                          verbose,
                          num_bootstraps = 2) {
  # to do: parallelize
  bootstraps = lapply(seq(1, num_bootstraps), function(i) {
    bootstrapped_model = bootstrap_model(ccs,
                                         full_model_formula_str,
                                         best_full_model,
                                         reduced_pln_model,
                                         pseudocount,
                                         initial_penalties,
                                         pln_min_ratio,
                                         pln_num_penalties,
                                         random.seed = i)

    best_bootstrapped_model = PLNmodels::getModel(bootstrapped_model, var=best_reduced_model$penalty)
    coef(best_bootstrapped_model) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_group")
  })

  # compute the covariance of the parameters
  get_cov_mat = function(data, cell_group) {

    cov_matrix = cov(data)
    rownames(cov_matrix) = paste0(cell_group, "_", rownames(cov_matrix))
    colnames(cov_matrix) = paste0(cell_group, "_", colnames(cov_matrix))
    return(cov_matrix)

  }

  bootstrapped_df = do.call(rbind, bootstraps) %>%
    group_by(cell_group) %>%
    tidyr::nest() %>%
    mutate(cov_mat = purrr::map2(.f = get_cov_mat,
                                 .x = data,
                                 .y = cell_group))

  bootstrapped_vhat = Matrix::bdiag(bootstrapped_df$cov_mat) %>% as.matrix()
  names = lapply(bootstrapped_df$cov_mat, function(m){ colnames(m)}) %>% unlist()
  rownames(bootstrapped_vhat) = names
  colnames(bootstrapped_vhat) = names


  return(bootstrapped_vhat)
}


#' Create a new cell_count_model object.
#'
#' Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a
#' way of encoding a graph prior. Automatically selects sparsity parameter, but
#' allows user to update it.
#'
#' @param ccs A Hooke cell_count_set object.
#' @param main_model_formula_str A character string specifying the model of cell abundances across samples,
#' where terms refer to columns in\code{colData(ccs)}. Put main effects here.
#' @param nuisance_model_formula_str A character string specifying the model of cell abundances across samples. Put nuisance effects here.
#' @param penalty_matrix A numeric NxN symmetric matrix specifying penalties for
#'   the PLN model, where N is the number of cell types. Entries must be
#'   positive. Use to specify an undirected graph prior for the PLN model.
#' @param sparsity_factor A positive number to control how sparse the PLN network is. Larger values make the network more sparse.
#' @return a new cell_count_model object
#' @importFrom PLNmodels prepare_data
#' @importFrom PLNmodels PLNnetwork
#' @importFrom PLNmodels getBestModel
#' @importFrom PLNmodels getModel
#' @export
new_cell_count_model <- function(ccs,
                                 main_model_formula_str,
                                 nuisance_model_formula_str = "1",
                                 penalty_matrix = NULL,
                                 whitelist=NULL,
                                 blacklist=NULL,
                                 sparsity_factor=0.1,
                                 base_penalty = 1,
                                 min_penalty=0.01,
                                 max_penalty=1e6,
                                 verbose=FALSE,
                                 pseudocount=0,
                                 pln_min_ratio=0.001,
                                 pln_num_penalties=30,
                                 norm_method = c("size_factors","TSS", "CSS",
                                                 "RLE", "GMPR", "Wrench", "none"),
                                 size_factors = NULL,
                                 num_bootstraps = NULL,
                                 inception = NULL,
                                 ...) {

  norm_method <- match.arg(norm_method)
  if (norm_method == "size_factors") {
    if (!is.null(size_factors)) {

      assertthat::assert_that(
        tryCatch(expr = identical(sort(colnames(ccs)), sort(names(size_factors))),
                 error = function(e) FALSE),
        msg = "size factors don't match")

      pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                          covariates = colData(ccs) %>% as.data.frame,
                                          offset = size_factors)
    } else {
      pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                          covariates = colData(ccs) %>% as.data.frame,
                                          offset = monocle3::size_factors(ccs))
    }
  } else {
    pln_data <- PLNmodels::prepare_data(counts = counts(ccs) + pseudocount,
                                        covariates = colData(ccs) %>% as.data.frame,
                                        offset = norm_method)
  }



  main_model_formula_str = stringr::str_replace_all(main_model_formula_str, "~", "")
  nuisance_model_formula_str = stringr::str_replace_all(nuisance_model_formula_str, "~", "")

  full_model_formula_str = paste("Abundance~", main_model_formula_str, "+", nuisance_model_formula_str, " + offset(log(Offset))")
  full_model_formula = as.formula(full_model_formula_str)

  reduced_model_formula_str = paste("Abundance~", nuisance_model_formula_str, " + offset(log(Offset))")
  reduced_model_formula = as.formula(reduced_model_formula_str)
  #pln_data <- as.name(deparse(substitute(pln_data)))

  if (is.null(penalty_matrix)){
    initial_penalties = init_penalty_matrix(ccs, whitelist=whitelist, blacklist=blacklist, base_penalty=base_penalty,min_penalty=min_penalty, max_penalty=max_penalty, ...)
    initial_penalties = initial_penalties[colnames(pln_data$Abundance), colnames(pln_data$Abundance)]
  }else{
    # TODO: check and validate dimensions of the user-provided penaties
    initial_penalties = penalty_matrix
  }

  # FIXME: This might only actually work when grouping cells by clusters and cluster names are
  # integers. We should make sure this generalizes when making white/black lists of cell groups
  # by type or other groupings
  if (is.null(whitelist) == FALSE){
    initial_penalties[as.matrix(whitelist[,c(1,2)])] = min_penalty
    initial_penalties[as.matrix(whitelist[,c(2,1)])] = min_penalty

  }

  if (is.null(blacklist) == FALSE){
    initial_penalties[as.matrix(blacklist[,c(1,2)])] = max_penalty
    initial_penalties[as.matrix(blacklist[,c(2,1)])] = max_penalty
  }

  #penalty_matrix = penalty_matrix[row.names(counts(ccs)), row.names(counts(ccs))]
  initial_penalties = initial_penalties[colnames(pln_data$Abundance), colnames(pln_data$Abundance)]

  # INSANE R BULLSHIT ALERT: for reasons I do not understand,
  # calling the fit via do.call prevents a weird error with formula
  # created with as.formula (e.g. after pasting).

  reduced_pln_model <- do.call(PLNmodels::PLNnetwork, args=list(reduced_model_formula_str,
                                                                data=pln_data,
                                                                control_init=list(min.ratio=pln_min_ratio, nPenalties=pln_num_penalties),
                                                                control_main=list(penalty_weights=initial_penalties,
                                                                                  trace = ifelse(verbose, 2, 0),
                                                                                  inception = inception),
                                                                ...),)


  full_pln_model <- do.call(PLNmodels::PLNnetwork, args=list(full_model_formula_str,
                                                               data=pln_data,
                                                               penalties = reduced_pln_model$penalties,
                                                               control_init=list(min.ratio=pln_min_ratio, nPenalties=pln_num_penalties),
                                                               control_main=list(penalty_weights=initial_penalties,
                                                                                 trace = ifelse(verbose, 2, 0),
                                                                                 inception = inception),
                                                               ...),)

  model_frame = model.frame(full_model_formula[-2], pln_data)
  xlevels = .getXlevels(terms(model_frame), model_frame)

  # Choose a model that isn't very aggressively sparsified
  best_reduced_model <- PLNmodels::getBestModel(reduced_pln_model, "EBIC")
  best_reduced_model <- PLNmodels::getModel(reduced_pln_model, var=best_reduced_model$penalty * sparsity_factor)

  # Choose a model that isn't very aggressively sparsified
  #best_full_model <- PLNmodels::getBestModel(full_pln_model, "EBIC")
  best_full_model <- PLNmodels::getModel(full_pln_model, var=best_reduced_model$penalty)

  if (is.null(num_bootstraps)) {
    bootstrapped_vhat = matrix(, nrow = 1, ncol = 1)
  } else {
    bootstrapped_vhat = bootstrap_vhat(ccs,
                                       full_model_formula_str,
                                       best_full_model,
                                       best_reduced_model,
                                       reduced_pln_model,
                                       pseudocount,
                                       initial_penalties,
                                       pln_min_ratio,
                                       pln_num_penalties,
                                       verbose,
                                       num_bootstraps)
  }

  ccm <- methods::new("cell_count_model",
                      ccs = ccs,
                      full_model_formula = full_model_formula,
                      best_full_model = best_full_model,
                      full_model_family = full_pln_model,
                      reduced_model_formula = reduced_model_formula,
                      best_reduced_model = best_reduced_model,
                      reduced_model_family = reduced_pln_model,
                      sparsity = sparsity_factor,
                      model_aux = SimpleList(model_frame=model_frame, xlevels=xlevels),
                      bootstrapped_vhat = bootstrapped_vhat
                      )
  #
  # metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  # clusters <- stats::setNames(SimpleList(), character(0))
  # cds <- estimate_size_factors(cds)
  # cds
  #ccm@model_aux[["best_model"]] = best_model
  #ccm@model_aux[["model_family"]] = pln_model


  ccm
}

#' Select the model a cell_count_model should use
#'
#' @param ccm A cell_count_model object.
#' @param criterion a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".
#' @return an updated cell_count_model object
#' @importFrom PLNmodels getBestModel
#' @importFrom PLNmodels getModel
#' @export
select_model <- function(ccm, criterion = "EBIC", sparsity_factor=1.0, models_to_update = c("both", "full", "reduced"))
{
  models_to_update = match.arg(models_to_update)

  #if (models_to_update == "reduced" || models_to_update == "both"){
    base_reduced_model <- PLNmodels::getBestModel(ccm@reduced_model_family, criterion)
    best_reduced_model <- PLNmodels::getModel(ccm@reduced_model_family, var=base_reduced_model$penalty * sparsity_factor)
    ccm@best_reduced_model = best_reduced_model
  #}

  #if (models_to_update == "full" || models_to_update == "both"){
  #  best_full_model <- PLNmodels::getBestModel(ccm@full_model_family, criterion)
    best_full_model <- PLNmodels::getModel(ccm@full_model_family, var=base_reduced_model$penalty * sparsity_factor)
    ccm@best_full_model = best_full_model
  #}

  ccm@sparsity = sparsity_factor

  return(ccm)
}


#' Initialize the PLN network penalty matrix, accepting optional whitelists and
#' blacklists of edges that are "free" or "off limits" between cell groups
#' @param ccs A cell_count_set of aggregated cell counts
#' @param whitelist a data frame with two columns corresponding to (undirected) edges that should receive no penalty
#' @param blacklist a data frame with two columns corresponding to (undirected) edges that should receive very high penalty
#' @param dist_fun A function that returns a penalty based given a distance between two clusters
init_penalty_matrix = function(ccs, whitelist=NULL, blacklist=NULL, base_penalty = 1, min_penalty=0.01, max_penalty=1e6){
  cell_group_centroids = centroids(ccs)
  dist_matrix = as.matrix(dist(cell_group_centroids[,-1], method = "euclidean", upper=T, diag = T))

  row.names(dist_matrix) <- cell_group_centroids$cell_group
  colnames(dist_matrix) <- cell_group_centroids$cell_group

  # TODO: do I need this? Probably the caller can and should do this.
  #dist_matrix = dist_matrix[colnames(data$Abundance), colnames(data$Abundance)]

  get_rho_mat <- function(DM, distance_parameter = 1, s=1, xmin = NULL) {
    if (is.null(xmin)){
      xmin = min(DM[DM > 0]) / 2
    }
    #out <- (1-(xmin/DM)^s) * distance_parameter
    out = min_penalty + (DM / max(DM))^2
    #out =  1 + DM^s
    # penalties have to be > 0
    out[!is.finite(out)] <- min_penalty
    out[out < 0] <- min_penalty
    return(out)
  }

  penalty_matrix = base_penalty * (min_penalty + get_rho_mat(dist_matrix, distance_parameter=1, s=2))

  # TODO: add support for whitelisting and blacklisting
  #qplot(as.numeric(dist_matrix), as.numeric(out))
  return(penalty_matrix)
}



#' Builds a model formula for time series models based on the range of the data
#'
#' This is just a utility function that puts the knots in reasonable positions based on the range of the data
#' @export
build_interval_formula <- function(ccs, num_breaks, interval_var="timepoint", interval_start=NULL, interval_stop=NULL){
  if (is.null(interval_start)){
    interval_start = as.numeric(min(colData(ccs@cds)[,interval_var]))
  }
  if (is.null(interval_stop)){
    interval_stop = as.numeric(max(colData(ccs@cds)[,interval_var]))
  }

  interval_breakpoints = seq(interval_start, interval_stop, length.out=num_breaks)
  interval_breakpoints = interval_breakpoints[2:(length(interval_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
  interval_formula_str = paste("~ ns(", interval_var, ", knots=", paste("c(",paste(interval_breakpoints, collapse=","), ")", sep=""), ")")
  return(interval_formula_str)
}
#debug(build_interval_formula)

