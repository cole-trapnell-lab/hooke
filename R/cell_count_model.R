setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))
setOldClass(c("PLNnetworkfamily"), prototype=structure(list(), class="PLNnetworkfamily"))
setOldClass(c("PLNnetworkfit"), prototype=structure(list(), class="PLNnetworkfit"))
setOldClass(c("PLNfamily"), prototype=structure(list(), class="PLNfamily"))
setOldClass(c("PLNfit"), prototype=structure(list(), class="PLNfit"))

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
                   cds_coldata = "tbl_df",
                   cds_reduced_dims = "SimpleList",
                   info = "SimpleList")
)
setMethod("is.na", "cell_count_set", function(x) FALSE)


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
                   full_model_family = "ANY", # this is probably unsafe
                   best_full_model = "ANY", # this is probably unsafe
                   reduced_model_formula = "formula",
                   reduced_model_family = "PLNnetworkfamily",
                   best_reduced_model = "PLNnetworkfit",
                   sparsity = "numeric",
                   model_aux = "SimpleList",
                   vhat = "dgCMatrix",
                   vhat_method = "character")
)
setMethod("is.na", "cell_count_model", function(x) FALSE)

empty_sparse_matrix = function (nrow = 0L, ncol = 0L, format = "R", dtype = "d")
{
  if (NROW(format) != 1L || !(format %in% c("R", "C", "T")))
    stop("'format' must be one of 'R', 'C', 'T'.")
  if (NROW(dtype) != 1L || !(dtype %in% c("d", "l", "n")))
    stop("'dtype' must be one of 'd', 'l', 'n'.")
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  if (NROW(nrow) != 1L || is.na(nrow) || nrow < 0)
    stop("'nrow' must be a non-negative integer.")
  if (NROW(ncol) != 1L || is.na(ncol) || ncol < 0)
    stop("'ncol' must be a non-negative integer.")
  target_class <- sprintf("%sg%sMatrix", dtype, format)
  out <- new(target_class)
  out@Dim <- as.integer(c(nrow, ncol))
  if (format == "R") {
    out@p <- integer(nrow + 1L)
  }
  else if (format == "C") {
    out@p <- integer(ncol + 1L)
  }
  return(out)
}

#' Create a new cell_data_set object.
#'
#' @param cds A Monocle cell data set object.
#' @param sample_group A column in colData(cds) that specifes how
#'    cells are grouped into samples.
#' @param cell_group A column in colData(cds) that specifies how
#'    cells are grouped into types or states (e.g. cluster).
#' @param sample_metadata Data frame containing attributes of
#'    individual samples, where the column named 'sample' has
#'    entries in \code{sample_group}.
#' @param cell_metadata Data frame containing attributes of
#'    individual cell groups, where \code{row.names(cell_metadata)}
#'    are entries in \code{cell_group}
#' @param lower_threshold numeric Minimum number of cells in
#'    retained cell_groups.
#' @param upper_threshold numeric Maximum number of cells in
#'    retained cell_groups.
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
                               upper_threshold = NULL,
                               keep_cds=TRUE) {

  assertthat::assert_that(is(cds, 'cell_data_set'),
                          msg = paste('Argument cds must be a cell_data_set.'))

  assertthat::assert_that(sample_group %in% colnames(colData(cds)),
                          msg = paste('Argument sample_group value must be a column name in the cell_data_set.'))

  assertthat::assert_that(cell_group %in% colnames(colData(cds)),
                          msg = paste('Argument cell_group value must be a column name in the cell_data_set.'))

  assertthat::assert_that(is.null(sample_metadata) || is.data.frame(sample_metadata),
                          msg = paste('Argument sample_metadata must be a data frame.'))

  sample_group_names_cds <- unique(colData(cds)[[sample_group]])
  assertthat::assert_that(is.null(sample_metadata) || nrow(sample_metadata) == length(sample_group_names_cds),
                          msg = paste('Argument sample_metadata must have the same',
                                      'number of rows as there are distinct sample',
                                      'names in the cds column data.'))

  assertthat::assert_that(is.null(sample_metadata) || all(sample_group_names_cds %in% sample_metadata[['sample']]),
                          msg = paste('Argument sample_metadata must have sample group names in',
                                       'a column called \'sample\'.'))

  assertthat::assert_that(is.null(cell_metadata) || is.data.frame(cell_metadata),
                          msg = paste('Argument cell_metadata must be a data frame.'))

  cell_group_names_cds <- unique(colData(cds)[[cell_group]])
  assertthat::assert_that(is.null(cell_metadata) || nrow(cell_metadata) == length(cell_group_names_cds),
                          msg = paste('Argument cell_metadata must have the same',
                                      'number of rows as there are distinct cell_group',
                                      'names in the cds column data.'))

  assertthat::assert_that(is.null(cell_metadata) || all(cell_group_names_cds %in% row.names(cell_metadata)),
                          msg = paste('Argument cell_metadata row names must contain the cell_group',
                                       'names.'))

  assertthat::assert_that(is.null(lower_threshold) || is.numeric(lower_threshold),
                          msg = paste('Argument lower_threshold must be numeric.'))

  assertthat::assert_that(is.null(upper_threshold) || is.numeric(upper_threshold),
                          msg = paste('Argument upper_threshold must be numeric.'))

  if(sample_group != 'sample')
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

  coldata_df = coldata_df %>% dplyr::rename("sample" = sample_group, "cell_group" = as.character(cell_group))
  #coldata_df$cell_group = factor(coldata_df$cell_group, levels=unique(colData(cds)[,cell_group]))

  coldata_df$group_id = coldata_df %>%
    dplyr::group_by(sample, cell_group) %>%
    dplyr::group_indices() %>% as.character

  # add to cds
  colData(cds)$group_id = coldata_df$group_id

  cds_summary = coldata_df %>%
    dplyr::group_by(sample, cell_group) %>%
    dplyr::summarize(cells = dplyr::n())

  cds_covariates_df = coldata_df %>%
    dplyr::select(-cell_group) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(across(where(is.numeric), function(x){mean(x)}),
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

  # This is super confusing because of the way the arguments are
  # named in new_cell_data_set. We are making a matrix of
  # dimension MxN, where M are cell types and N are samples
  # (e.g. embryos, replicates, etc). The "gene" metadata monocle
  # normally expects will actually be used to hold cell group
  # metadata.

  ccs_cds = cds[, !colData(cds)[[cell_group]] %in% removed_cell_states]

  # TODO: We could probably avoid duplicating this info if keep_cds == TRUE, providing it
  # through accessor functions directly from the cds

  # FIXME: potentially we should be using the filtered one above? Potentially rename cell_group, sample, etc?
  cds_coldata = colData(ccs_cds) %>% as_tibble
  cds_reducedDims = reducedDims(ccs_cds)

  cell_metadata_subset <- cell_metadata[rownames(cell_counts_wide),,drop=FALSE]
  ccs = methods::new("cell_count_set",
               monocle3::new_cell_data_set(cell_counts_wide,
                                           cell_metadata=cds_covariates_df,
                                           gene_metadata=cell_metadata_subset),
               cds=ccs_cds,
               cds_coldata=cds_coldata,
               cds_reduced_dims=cds_reducedDims,
               info=SimpleList(sample_group=sample_group,
                               cell_group=cell_group))

  if (keep_cds == FALSE)
    ccs@cds = new_cell_data_set(empty_sparse_matrix(format="C"))

  # if (!is.null(cell_metadata)) {
  #   assertthat::assert_that(!is.null(row.names(cell_metadata)) &
  #                             all(row.names(cell_metadata) == colnames(expression_data)),
  #                           msg = paste("row.names of cell_metadata must be equal to colnames of",
  #                                       "expression_data"))
  # }
  #
  # if (!is.null(gene_metadata)) {
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

  # Notes:
  #   o  ccs_cds has the original column names whereas coldata_df has
  #      several renamed columns.
  #   o  coldata_df has all rows
  #   o  ccs_cds has rows filtered by thresholds
  ccs@metadata[["cell_group_assignments"]] = coldata_df %>% dplyr::select(group_id, sample, cell_group) %>% as.data.frame
  ccs@metadata[["cell_group_assignments"]] = ccs@metadata[["cell_group_assignments"]] %>% filter(!cell_group %in% removed_cell_states)
  row.names(ccs@metadata[["cell_group_assignments"]]) = colnames(ccs_cds)

  return (ccs)
}


#' resamples the ccs counts using a multinomial distribution
#' @param ccs
#' @param random.seed
#' @noRd
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
#' @noRd
bootstrap_model = function(ccs,
                           full_model_formula_str,
                           best_full_model,
                           reduced_pln_model,
                           pseudocount,
                           initial_penalties,
                           pln_min_ratio,
                           pln_num_penalties,
                           random.seed,
                           norm_method,
                           backend = c('nlopt', 'torch')) {

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(backend) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument backend must be one of "nlopt" or "torch".'))
  backend <- match.arg(backend)

  # resample the counts
  sub_ccs = bootstrap_ccs(ccs, random.seed = random.seed)

  if (norm_method == "size_factors") {
    norm_method = monocle3::size_factors(sub_ccs)
  }

  # remake data from new ccs
  sub_pln_data <- PLNmodels::prepare_data(counts = counts(sub_ccs) + pseudocount,
                                          covariates = colData(sub_ccs) %>% as.data.frame,
                                          offset = norm_method)
  # rerun the model using the same initial parameters
  # as the original, non bootstrapped model

  if (backend == "torch") {
# bge (20221227): notes:
#                   o I am trying to track the code in the PLNmodels master branch at Github
#                   o the PLNmodels::PLN() function versions that I see do not include
#                     the arguments min.ratio, nPenalties, vcov_est, and penalty_weights. I
#                     am confused...
#                   o I revert to the original because the PLNmodels changes break hooke.
    sub_full_model = do.call(PLNmodels::PLN, args=list(full_model_formula_str,
                                                       data = sub_pln_data,
                                                       # penalties = reduced_pln_model$penalties,
                                                       control = PLNmodels::PLN_param(backend = 'torch',
                                                                                      trace = ifelse(FALSE, 2, 0),
                                                                                      inception = best_full_model,
                                                                                      config_optim = list(maxevel  = 10000,
                                                                                                          ftol_rel = 1e-8,
                                                                                                          xtol_rel = 1e-6))))

# bge (20221227): notes
#                   o I restore the PLN call for the earlier PLNmodels version commit 022d59d
#                   o the earlier version of PLNmodels was 'PLNmodels    * 0.11.7-9600 2022-11-29 [1] Github (PLN-team/PLNmodels@022d59d)'
#     sub_full_model = do.call(PLNmodels::PLN, args=list(full_model_formula_str,
#                                                        data = sub_pln_data,
#                                                        # penalties = reduced_pln_model$penalties,
#                                                        control = list(min.ratio=pln_min_ratio,
#                                                                       nPenalties=pln_num_penalties,
#                                                                       vcov_est = "none",
#                                                                       inception = best_full_model,
#                                                                       backend = "torch",
#                                                                       penalty_weights=initial_penalties,
#                                                                       trace = ifelse(FALSE, 2, 0))))

  } else {
# bge (20221227): notes:
#                   o I am trying to track the code in the PLNmodels master branch at Github
#                   o I revert to the original because the PLNmodels changes break hooke.
    sub_full_model = do.call(PLNmodels::PLNnetwork, args=list(full_model_formula_str,
                                                              data = sub_pln_data,
                                                              penalties = reduced_pln_model$penalties,
                                                              control = PLNmodels::PLNnetwork_param(backend = 'nlopt',
                                                                                                    trace = ifelse(FALSE, 2, 0),
                                                                                                    n_penalties = pln_num_penalties,
                                                                                                    min_ratio = pln_min_ratio,
                                                                                                    penalty_weights = initial_penalties,
                                                                                                    config_optim = list(algorithm = 'CCSAQ',
                                                                                                                        maxeval = 10000,
                                                                                                                        ftol_rel = 1e-8,
                                                                                                                        xtol_rel = 1e-6,
                                                                                                                        ftol_out = 1e-6,
                                                                                                                        maxit_out = 50,
                                                                                                                        ftol_abs = 0.0,
                                                                                                                        xtol_abs = 0.0,
                                                                                                                        maxtime = -1))))

# bge (20221227): notes
#                   o I restore the PLN call for the earlier PLNmodels version commit 022d59d
#                   o the earlier version of PLNmodels was 'PLNmodels    * 0.11.7-9600 2022-11-29 [1] Github (PLN-team/PLNmodels@022d59d)'
#     sub_full_model = do.call(PLNmodels::PLNnetwork, args=list(full_model_formula_str,
#                                                               data = sub_pln_data,
#                                                               penalties = reduced_pln_model$penalties,
#                                                               control_init=list(min.ratio=pln_min_ratio,
#                                                                                 nPenalties=pln_num_penalties,
#                                                                                 penalty_weights=initial_penalties),
#                                                               control_main=list(trace = ifelse(FALSE, 2, 0))))
  }

  return(sub_full_model)
}


#
# Removed 20230104 bge PLNmodels v1.0.0 removed get_vcov_hat
#' @noRd
# compute_vhat = function(model, model_family, type) {
#
#     if (model$d > 0) {
#       ## self$fisher$mat : Fisher Information matrix I_n(\Theta) = n * I(\Theta)
#       ## safe inversion using Matrix::solve and Matrix::diag and error handling
#
#       X = model_family$responses
#       Y = model_family$covariates
#       model$get_vcov_hat(type, X, Y)
#
#       vhat = vcov(model)
#
#       # zero out everything not on block diagonal
#       if (type == "sandwich") {
#
#         num_coef = dim(coef(model))[2]
#         num_blocks =dim(vhat)[1]/num_coef
#
#         blocks = lapply(1:num_blocks, function(i) {
#           start = num_coef*(i-1) + 1
#           end = num_coef*i
#           vhat[start:end,start:end]
#         })
#
#         vhat = Matrix::bdiag(blocks) %>% as.matrix()
#       }
#
#       # vcov_mat = vcov(model)
#
#     #   vhat <- matrix(0, nrow = nrow(vcov_mat), ncol = ncol(vcov_mat))
#     #
#     #   #dimnames(vhat) <- dimnames(vcov_mat)
#     #   safe_rows = safe_cols = Matrix::rowSums(abs(vcov_mat)) > 0
#     #   vcov_mat = vcov_mat[safe_rows, safe_cols]
#     #
#     #   out <- tryCatch(Matrix::solve(vcov_mat),
#     #                   error = function(e) {e})
#     #   row.names(out) = colnames(out) = names(safe_rows[safe_rows])
#     #   if (is(out, "error")) {
#     #     warning(paste("Inversion of the Fisher information matrix failed with following error message:",
#     #                   out$message,
#     #                   "Returning NA",
#     #                   sep = "\n"))
#     #     vhat <- matrix(NA, nrow = model$p, ncol = model$d)
#     #   } else {
#     #     row.names(out) = colnames(out) = names(safe_rows[safe_rows])
#     #     row.names(vhat) = colnames(vhat) = row.names(vcov(model))
#     #     vhat[safe_rows, safe_cols] = as.numeric(out) #as.numeric(out) #%>% sqrt %>% matrix(nrow = self$d) %>% t()
#     #   }
#     #   #dimnames(vhat) <- dimnames(vcov_mat)
#     } else {
#       vhat <- NULL
#     }
#     vhat = methods::as(vhat, "dgCMatrix")
#     vhat
# }


#' computes the avg vhat across n bootstraps
#' @param ccm
#' @param num_bootstraps
#' @noRd
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
                          norm_method,
                          num_bootstraps,
                          backend) {
  # to do: parallelize

  get_bootstrap_coef = function(random.seed,
                                ccs,
                                full_model_formula_str,
                                best_full_model,
                                reduced_pln_model,
                                pseudocount,
                                initial_penalties,
                                pln_min_ratio,
                                pln_num_penalties,
                                norm_method,
                                backend) {

    bootstrapped_model = bootstrap_model(ccs,
                                         full_model_formula_str,
                                         best_full_model,
                                         reduced_pln_model,
                                         pseudocount,
                                         initial_penalties,
                                         pln_min_ratio,
                                         pln_num_penalties,
                                         random.seed = random.seed,
                                         norm_method = norm_method,
                                         backend = backend)

    if (backend == "torch") {
      coef(bootstrapped_model) %>%
        as.data.frame() %>% t() %>%
        tibble::rownames_to_column("cell_group")
    } else {
      best_bootstrapped_model = PLNmodels::getModel(bootstrapped_model, var=best_reduced_model$penalty)
      coef(best_bootstrapped_model) %>% t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("cell_group")
    }
  }

  coef_df = data.frame(seed = seq(1, num_bootstraps)) %>%
    mutate(coef = purrr::map(.f = purrr::possibly(get_bootstrap_coef, NA_character_),
                             .x = seed,
                             ccs,
                             full_model_formula_str,
                             best_full_model,
                             reduced_pln_model,
                             pseudocount,
                             initial_penalties,
                             pln_min_ratio,
                             pln_num_penalties,
                             norm_method,
                             backend))
  coef_df = coef_df %>% filter(!is.na(coef))

  # compute the covariance of the parameters
  get_cov_mat = function(data, cell_group) {

    cov_matrix = cov(data)
    rownames(cov_matrix) = paste0(cell_group, "_", rownames(cov_matrix))
    colnames(cov_matrix) = paste0(cell_group, "_", colnames(cov_matrix))
    return(cov_matrix)
  }

  bootstrapped_df = do.call(rbind, coef_df$coef) %>%
    group_by(cell_group) %>%
    tidyr::nest() %>%
    mutate(cov_mat = purrr::map2(.f = get_cov_mat,
                                 .x = data,
                                 .y = cell_group))

  bootstrapped_vhat = Matrix::bdiag(bootstrapped_df$cov_mat) %>% as.matrix()
  names = lapply(bootstrapped_df$cov_mat, function(m){ colnames(m)}) %>% unlist()
  rownames(bootstrapped_vhat) = names
  colnames(bootstrapped_vhat) = names

  bootstrapped_vhat = methods::as(bootstrapped_vhat, "dgCMatrix")
  return(bootstrapped_vhat)
}


#' Create a new cell_count_model object.
#'
#' Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a
#' way of encoding a graph prior. Automatically selects sparsity parameter, but
#' allows user to update it.
#'
#' @param ccs A Hooke cell_count_set object.
#' @param main_model_formula_str A character string specifying the model of cell
#'    abundances across samples,
#'   where terms refer to columns in\code{colData(ccs)}. Put main effects here.
#' @param nuisance_model_formula_str A character string specifying the model of cell
#'    abundances across samples. Put nuisance effects here.
#' @param penalty_matrix A numeric NxN symmetric matrix specifying penalties for
#'   the PLN model, where N is the number of cell types. Entries must be
#'   positive and the rows and columns must be named with the cell_group names.
#'   Use to specify an undirected graph prior for the PLN model.
#' @param whitelist list A data frame with two columns corresponding to (undirected)
#'    edges that should receive min_penalty. The columns are integers that refer to
#'    cell clusters.
#' @param blacklist list A data frame with two columns corresponding to (undirected)
#'    edges that should receive max_penalty. The columns are integers that refer to
#'    cell clusters.
#' @param sparsity_factor A positive number to control how sparse the PLN network is. Larger values make the network more sparse.
#'    edges that should receive min_penalty. The columns are either cell_group
#'    names or integers that refer to cell_groups in penalty_matrix.
#' @param base_penalty numeric A factor that scales the penalty matrix.
#' @param min_penalty numeric A positive value that is assigned to whitelisted
#'    penalty matrix elements, which over-write existing values.
#' @param max_penalty numeric A positive value that is assigned to blacklisted
#'    penalty matrix elements. which over-write existing values.
#' @param verbose logical Whether to emit verbose output.
#' @param pseudocount integer A value added to the elements of the initial
#'    cell_count_set matrix.
#' @param pln_min_ratio numeric Used in the definition of the sparsity penalty grid.
#' @param pln_num_penalties integer Number of penalty values for the internally
#'    generated penalty grid.
#' @param norm_method string Normalization method used to compute scaling
#'    factors used as offset during PLN inference.
#' @param vhat_method string Method used to compute covariance matrix?
#' @param size_factors numeric vector or matrix User supplied vector or matrix of
#'    offsets passed the PLNmodels::prepare_data() method.
#' @param num_bootstraps positive integer Number of iterations used with the
#'    bootstrap vhat_method.
#' @param inception Not used.
#' @param backend Method used to run bootstrap iterations.
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
                                 vhat_method = c("variational_var", "jackknife", "bootstrap", "my_bootstrap"),
                                 covariance_type = c("spherical", "diagonal"),
                                 size_factors = NULL,
                                 num_bootstraps = 10,
                                 inception = NULL,
                                 backend = c("nlopt", "torch"),
                                 num_threads=1,
                                 ftol_rel = 1e-6,
                                 penalize_by_distance=TRUE,
                                 ...) {

  assertthat::assert_that(is(ccs, 'cell_count_set'))

  assertthat::assert_that(assertthat::is.string(main_model_formula_str))
  assertthat::assert_that(assertthat::is.string(nuisance_model_formula_str))

  assertthat::assert_that(assertthat::is.number(sparsity_factor) && sparsity_factor >= 0.0)
  assertthat::assert_that(assertthat::is.number(base_penalty) && base_penalty >= 0.0)
  assertthat::assert_that(assertthat::is.number(min_penalty) && min_penalty >= 0.0)
  assertthat::assert_that(assertthat::is.number(max_penalty) && max_penalty > 0.0)
  assertthat::assert_that(min_penalty <= max_penalty)
  assertthat::assert_that(assertthat::is.number(pseudocount) && pseudocount >= 0)
  assertthat::assert_that(assertthat::is.number(pln_min_ratio) && pln_min_ratio > 0.0)
  assertthat::assert_that(assertthat::is.count(pln_num_penalties) && pln_num_penalties >= 0)
  assertthat::assert_that(assertthat::is.count(num_bootstraps) && num_bootstraps >= 0)

  # Check that penalty matrix dimensions are N x N where N is the number of cell types.
  # Check that the penalty matrix has row/columm names that belong to the cell_groups,
  # which may be used for the whitelist and blacklist, and for ordering the matrix rows
  # and columns. See the line:
  #   initial_penalties = initial_penalties[colnames(pln_data$Abundance), colnames(pln_data$Abundance)]
  #
  ccs_num_cell_group <- dim(counts(ccs))[[1]]
  ccs_cell_group_names <- rownames(counts(ccs))

  assertthat::assert_that(is.null(penalty_matrix) || ( is.matrix(penalty_matrix) && is.numeric(penalty_matrix[1][1])))
  assertthat::assert_that(is.null(penalty_matrix) || identical(dim(penalty_matrix), c(ccs_num_cell_group, ccs_num_cell_group)))
  assertthat::assert_that(is.null(penalty_matrix) || (!is.null(rownames(penalty_matrix)) && all(rownames(penalty_matrix) %in% ccs_cell_group_names)))
  assertthat::assert_that(is.null(penalty_matrix) || (!is.null(rownames(penalty_matrix)) && all(colnames(penalty_matrix) %in% ccs_cell_group_names)))

  # Check the whitelist and blacklist for expected values.
  assertthat::assert_that(is.null(whitelist) || (is.numeric(whitelist[[1]][[1]]) &&
                                                 range(whitelist[[1]])[[1]] >= 0 &&
                                                 range(whitelist[[1]])[[2]] <= ccs_num_cell_group) ||
                                                (is.character(whitelist[[1]][[1]]) &&
                                                 all(whitelist[[1]] %in% ccs_cell_group_names)))

  assertthat::assert_that(is.null(whitelist) || (is.numeric(whitelist[[2]][[1]]) &&
                                                 range(whitelist[[2]])[[1]] >= 0 &&
                                                 range(whitelist[[2]])[[2]] <= ccs_num_cell_group) ||
                                                (is.character(whitelist[[2]][[1]]) &&
                                                 all(whitelist[[2]] %in% ccs_cell_group_names)))

  assertthat::assert_that(is.null(blacklist) || (is.numeric(blacklist[[1]][[1]]) &&
                                                 range(blacklist[[1]])[[1]] >= 0 &&
                                                 range(blacklist[[1]])[[2]] <= ccs_num_cell_group) ||
                                                (is.character(blacklist[[1]][[1]]) &&
                                                 all(blacklist[[1]] %in% ccs_cell_group_names)))

  assertthat::assert_that(is.null(blacklist) || (is.numeric(blacklist[[2]][[1]]) &&
                                                 range(blacklist[[2]])[[1]] >= 0 &&
                                                 range(blacklist[[2]])[[2]] <= ccs_num_cell_group) ||
                                                (is.character(blacklist[[2]][[1]]) &&
                                                 all(blacklist[[2]] %in% ccs_cell_group_names)))

  assertthat::assert_that(assertthat::is.flag(verbose))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument norm_method must be one of "size_factors",',
                '"TSS", "CSS", "RLE", "GMPR", "Wrench", or "none".'))
  norm_method <- match.arg(norm_method)

  # TO DO: FIX THIS
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(vhat_method) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste( 'Argument vhat_method must be one of "variational_var",',
                 '"jackknife","bootstrap", or "my_bootstrap".'))
  vhat_method <- match.arg(vhat_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(backend) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste( 'Argument backend must be one of "nlopt" or "torch".'))
  backend <- match.arg(backend)
  
  covariance_type = match.arg(covariance_type)

  #
  # PLNmodels::prepare_data returns (1) a matrix of cell abundances,
  # which were calculate in new_cell_count_set() where rows are
  # sample groups and the columns are cell groups, (2) covariates,
  # where is a copy of colData(cds), and (3) offsets, which are
  # calculated by PLNmodels::prepare_data.
  if (norm_method == "size_factors") {
    if (!is.null(size_factors)) {

      assertthat::assert_that(
        tryCatch(expr = identical(sort(colnames(ccs)), sort(names(size_factors))),
                 error = function(e) FALSE),
        msg = "Argument size factor names must match ccs column names.")

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

    if (norm_method == "none") {
      pln_data$Offset = 1
    }
  }

  main_model_formula_str = stringr::str_replace_all(main_model_formula_str, "~", "")
  nuisance_model_formula_str = stringr::str_replace_all(nuisance_model_formula_str, "~", "")

  full_model_formula_str = paste("Abundance~", main_model_formula_str, "+", nuisance_model_formula_str, " + offset(log(Offset))")
  # full_model_formula = as.formula(full_model_formula_str)
  full_model_formula <- tryCatch(
                          {
                            as.formula(full_model_formula_str)
                          },
                          error = function(condition) {
                            message(paste('Bad full_model_formula string', full_model_formula_str), ': ', condition, '.')
                          },
                          warn = function(condition) {
                            message(paste('Bad full_model_formula string', full_model_formula_str), ': ', condition, '.')
                          })

  reduced_model_formula_str = paste("Abundance~", nuisance_model_formula_str, " + offset(log(Offset))")
#  reduced_model_formula = as.formula(reduced_model_formula_str)
  reduced_model_formula <- tryCatch(
                          {
                            as.formula(reduced_model_formula_str)
                          },
                          error = function(condition) {
                            message(paste('Bad reduced_model_formula string', reduced_model_formula_str), ': ', condition, '.')
                          },
                          warn = function(condition) {
                            message(paste('Bad reduced_model_formula string', reduced_model_formula_str), ': ', condition, '.')
                          })

  #pln_data <- as.name(deparse(substitute(pln_data)))

  # Note: the whitelist and blacklist are applied later in this function
  #       (new_cell_count_model), not in init_penalty_matrix() so the
  #       arguments in the call to init_penalty_matrix() are unused there.
  #       This is so that the whitelist and blacklist penalties are applied
  #       to the user supplied penalty matrix.
    if (is.null(penalty_matrix)){
      if (penalize_by_distance){
        initial_penalties = init_penalty_matrix(ccs, whitelist=whitelist, blacklist=blacklist, base_penalty=base_penalty,min_penalty=min_penalty, max_penalty=max_penalty, ...)
        initial_penalties = initial_penalties[colnames(pln_data$Abundance), colnames(pln_data$Abundance)]
      }else{
        initial_penalties = NULL
      }
    }else{
      initial_penalties = penalty_matrix
    }

  # FIXME: This might only actually work when grouping cells by clusters and cluster names are
  # integers. We should make sure this generalizes when making white/black lists of cell groups
  # by type or other groupings.
  # Note: this appears to work when the cell contents are character strings of cell.
  #       group names.
  if (is.null(initial_penalties) == FALSE){
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
  }

  # INSANE R BULLSHIT ALERT: for reasons I do not understand,
  # calling the fit via do.call prevents a weird error with formula
  # created with as.formula (e.g. after pasting).

  tryCatch({
    if (num_threads > 1){
      print (paste("fitting model with", num_threads, "threads"))
      RhpcBLASctl::omp_set_num_threads(1)
      RhpcBLASctl::blas_set_num_threads(num_threads)
    }else{
      print (paste("fitting model with", 1, "threads"))
    }

    if (backend == "torch") {
      control_optim_args = list(
        # algorithm = 'CCSAQ',
        maxeval = 10000,
        ftol_rel = ftol_rel,
        xtol_rel = 1e-6)
    } else {
      control_optim_args = list(
        algorithm = 'CCSAQ',
        maxeval = 10000,
        ftol_rel = ftol_rel,
        xtol_rel = 1e-6,
        ftol_out = 1e-6,
        maxit_out = 50,
        ftol_abs = 0.0,
        xtol_abs = 0.0,
        maxtime = -1)

    }

    variational_var = TRUE
    sandwich_var = FALSE
    jackknife = FALSE
    bootstrap = FALSE
    my_bootstrap = FALSE

    if (vhat_method == "variational_var" | vhat_method == "my_bootstrap") {
      variational_var = TRUE
    }else{
      variational_var = FALSE # Don't compute the variational variance unless we have to, because it sometimes throws exceptions
    }

    if (vhat_method == "sandwich_var") {
      sandwich_var = TRUE
    }

    if (vhat_method == "jackknife") {
      jackknife = TRUE
    }

    if (vhat_method == "bootstrap") {
      bootstrap = num_bootstraps
    }


# bge (20221227): notes:
#                   o I am trying to track the code in the PLNmodels master branch at Github
#                   o I revert to the original because the PLNmodels changes break hooke.
    reduced_pln_model <- do.call(PLNmodels::PLNnetwork, args=list(reduced_model_formula_str,
                                                                  data=pln_data,
                                                                  control = PLNmodels::PLNnetwork_param(backend = backend,
                                                                                                        trace = ifelse(verbose, 2, 0),
                                                                                                        covariance = covariance_type,
                                                                                                        n_penalties = pln_num_penalties,
                                                                                                        min_ratio = pln_min_ratio,
                                                                                                        penalty_weights = initial_penalties,
                                                                                                        config_post = list(jackknife = FALSE,  # never jackknife the reduced model
                                                                                                                           bootstrap = FALSE, # never bootstrap the reduced model
                                                                                                                           variational_var = FALSE, # never compute variational variances on the reduced model
                                                                                                                           sandwich_var = FALSE,  # never bootstrap the reduced model
                                                                                                                           rsquared = FALSE),
                                                                                                        config_optim = control_optim_args),
                                                                  ...),)

    full_pln_model <- do.call(PLNmodels::PLN, args=list(full_model_formula_str,
                                                               data=pln_data,
                                                               control = PLNmodels::PLN_param(backend = backend,
                                                                                              covariance = covariance_type,
                                                                                              trace = ifelse(verbose, 2, 0),
                                                                                              config_post = list(jackknife = jackknife,
                                                                                                                 bootstrap = bootstrap,
                                                                                                                 variational_var = variational_var,
                                                                                                                 sandwich_var = sandwich_var,
                                                                                                                 rsquared = FALSE),
                                                                                              config_optim = control_optim_args),
                                                               ...),)


# bge (20221227): notes:
#                   o the previous version of PLNmodels was PLNmodels    * 0.11.7-9600 2022-11-29 [1] Github (PLN-team/PLNmodels@022d59d)
#   full_pln_model <- do.call(PLNmodels::PLNnetwork, args=list(full_model_formula_str,
#                                                                data=pln_data,
#                                                                penalties = reduced_pln_model$penalties,
#                                                                control_init=list(min.ratio=pln_min_ratio,
#                                                                                  nPenalties=pln_num_penalties,
#                                                                                  penalty_weights=initial_penalties),
#                                                                control_main=list(trace = ifelse(verbose, 2, 0)),
#                                                                ...),)
  }, finally = {
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
  })


  model_frame = model.frame(full_model_formula[-2], pln_data)
  xlevels = .getXlevels(terms(model_frame), model_frame)

  # Choose a model that isn't very aggressively sparsified
  best_reduced_model <- PLNmodels::getBestModel(reduced_pln_model, "EBIC")
  best_reduced_model <- PLNmodels::getModel(reduced_pln_model, var=best_reduced_model$penalty * sparsity_factor)

  # Choose a model that isn't very aggressively sparsified
  #best_full_model <- PLNmodels::getBestModel(full_pln_model, "EBIC")
  # best_full_model <- PLNmodels::getModel(full_pln_model, var=best_reduced_model$penalty)
  best_full_model <- full_pln_model

  if (vhat_method == "my_bootstrap") {
    vhat = bootstrap_vhat(ccs,
                          full_model_formula_str,
                          best_full_model,
                          best_reduced_model,
                          reduced_pln_model,
                          pseudocount,
                          initial_penalties,
                          pln_min_ratio,
                          pln_num_penalties,
                          verbose,
                          norm_method,
                          num_bootstraps,
                          backend)

  } else if (vhat_method == "jackknife" | vhat_method == "bootstrap") {

    vhat_coef = coef(best_full_model, type = "main")
    var_jack_mat = attributes(vhat_coef)[[paste0("vcov_", vhat_method)]]
    # var_jack_mat = var_jack %>% as.data.frame %>%
    #   tibble::rownames_to_column("term") %>%
    #   tidyr::pivot_longer(-term, names_to = "cell_group", values_to = "var") %>%
    #   arrange(cell_group) %>%
    #   mutate(rowname = paste0(cell_group, "_", term)) %>%
    #   select(-term, -cell_group) %>%
    #   mutate(colname = rowname) %>%
    #   tidyr::pivot_wider(values_from = var, names_from = colname) %>%
    #   replace(is.na(.), 0) %>%
    #   tibble::column_to_rownames("rowname") %>%
    #   as.matrix()
    vhat <- methods::as(var_jack_mat, "dgCMatrix")

  } else {
    vhat <- vcov(best_full_model, type= "main")
    vhat <- methods::as(vhat, "dgCMatrix")
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
                      vhat = vhat,
                      vhat_method = vhat_method
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
#' @param sparsity_factor A positive number to control how sparse the PLN network
#'    is. Larger values make the network sparser.
#' @param models_to_update string The model to update. Must be one of "both", "full", or "reduced".
#' @return an updated cell_count_model object.
#' @importFrom PLNmodels getBestModel
#' @importFrom PLNmodels getModel
#' @export
select_model <- function(ccm, criterion = c("BIC", "EBIC", "StARS"), sparsity_factor=1.0, models_to_update = c("both", "full", "reduced"))
{

  assertthat::assert_that(is(ccm, 'cell_count_model'))
  assertthat::assert_that(assertthat::is.number(sparsity_factor) && sparsity_factor >= 0.0)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(criterion) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument criterion must be one of "BIC","EBIC", or "StARS".'))
  criterion <- match.arg(criterion)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(models_to_update) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument models_to_update must be one of "both","full", or "reduced".'))
  models_to_update = match.arg(models_to_update)

  #if (models_to_update == "reduced" || models_to_update == "both"){
    base_reduced_model <- PLNmodels::getBestModel(ccm@reduced_model_family, criterion)
    best_reduced_model <- PLNmodels::getModel(ccm@reduced_model_family, var=base_reduced_model$penalty * sparsity_factor)
    ccm@best_reduced_model = best_reduced_model
  #}

  #if (models_to_update == "full" || models_to_update == "both"){
  #  best_full_model <- PLNmodels::getBestModel(ccm@full_model_family, criterion)
  #  best_full_model <- PLNmodels::getModel(ccm@full_model_family, var=base_reduced_model$penalty * sparsity_factor)
  #  ccm@best_full_model = best_full_model
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
#' @noRd
init_penalty_matrix = function(ccs, whitelist=NULL, blacklist=NULL, base_penalty = 1, min_penalty=0.01, max_penalty=1e6){
  cell_group_centroids = centroids(ccs)
  dist_matrix = as.matrix(dist(cell_group_centroids[,-1], method = "euclidean", upper=T, diag = T))

  row.names(dist_matrix) <- cell_group_centroids$cell_group
  colnames(dist_matrix) <- cell_group_centroids$cell_group

  # TODO: do I need this? Probably the caller can and should do this.
  #dist_matrix = dist_matrix[colnames(data$Abundance), colnames(data$Abundance)]

  get_rho_mat <- function(DM, s=2, xmin = NULL) {
    if (is.null(xmin)){
      xmin = min(DM[DM > 0]) / 2
    }
    #out <- (1-(xmin/DM)^s) * distance_parameter
    out = min_penalty + (DM / max(DM))^s
    #out =  1 + DM^s
    # penalties have to be > 0
    out[!is.finite(out)] <- min_penalty
    out[out < 0] <- min_penalty
    return(out)
  }

  penalty_matrix = base_penalty * (min_penalty + get_rho_mat(dist_matrix, s=2))

  # TODO: add support for whitelisting and blacklisting
  #qplot(as.numeric(dist_matrix), as.numeric(out))
  return(penalty_matrix)
}



#' Builds a model formula for time series models based on the range of the data.
#' This is a utility function that puts the knots in reasonable positions based on the range of the data.
#' @param numeric num_breaks Number of interval points.
#' @param character interval_var
#' @param numeric interval_start Interval start value.
#' @param numeric interval_stop Interval stop value.
#' @return An interval model formula.
#' @export
build_interval_formula <- function(ccs, num_breaks, interval_var="timepoint", interval_start=NULL, interval_stop=NULL){

  assertthat::assert_that(is(ccs, 'cell_count_set'))
  assertthat::assert_that(is.numeric(num_breaks))
  assertthat::assert_that(is.numeric(interval_start))
  assertthat::assert_that(is.numeric(interval_stop))
  assertthat::assert_that(is.character(interval_var))

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

