
# NOTE: This ignores offsets. DO NOT USE TO RECOVER FITTED VALUES FROM THE ORIGINAL DATA
my_plnnetwork_predict <- function (ccm, newdata, type = c("link", "response"), envir = parent.frame())
{
  type = match.arg(type)
  X <- model.matrix(terms(ccm@model_aux[["model_frame"]]), newdata,
                         xlev = ccm@model_aux[["xlevels"]])
  #O <- model.offset(ccm@model_aux[["model_frame"]])
  EZ <- tcrossprod(X, model(ccm)$model_par$Theta)
  #if (!is.null(O))
  #  EZ <- EZ + O
  EZ <- sweep(EZ, 2, 0.5 * diag(model(ccm)$model_par$Sigma), "+")
  colnames(EZ) <- colnames(model(ccm)$model_par$Sigma)
  results <- switch(type, link = EZ, response = exp(EZ))
  attr(results, "type") <- type
  results
}


#  subsamples the ccs and computes a vhat
#' @param ccm
#' @param prop
#' @param bycol

# subsampled_ccm = function(ccm, prop = 0.8, bycol = TRUE, random.seed = 42) {
#
#   set.seed(random.seed)
#   counts(ccm@ccs) = scuttle::downsampleMatrix(counts(ccm@ccs), prop=prop)
#
#   nuisance_model_formula = stringr::str_remove(format(ccm@reduced_model_formula), "Abundance ")
#   main_model_formula = stringr::str_remove(format(ccm@full_model_formula), "Abundance ")
#   main_model_formula = substr(main_model_formula, 1, nchar(main_model_formula) - nchar(nuisance_model_formula))
#   nuisance_model_formula = substr(nuisance_model_formula, 1, nchar(nuisance_model_formula) - nchar("+ offset(log(Offset))"))
#
#   sample_ccm = new_cell_count_model(ccm@ccs,
#                                     main_model_formula_str = main_model_formula,
#                                     nuisance_model_formula_str = nuisance_model_formula,
#                                     inception = ccm)
#
#
#   return(sample_ccm)
# }

subsample_ccs = function(ccs, prop = 0.8, random.seed = 42) {
  set.seed(random.seed)
  counts(ccs) = scuttle::downsampleMatrix(counts(ccs), prop=prop)
  return(ccs)
}

bootstrap_model = function(ccs,
                           full_model_formula_str,
                           best_full_model,
                           reduced_pln_model,
                           pseudocount,
                           initial_penalties,
                           pln_min_ratio,
                           pln_num_penalties,
                           random.seed) {

  sub_ccs = subsample_ccs(ccs, random.seed = random.seed)
  sub_pln_data <- PLNmodels::prepare_data(counts = counts(sub_ccs) + pseudocount,
                                      covariates = colData(sub_ccs) %>% as.data.frame,
                                      offset = size_factors(sub_ccs))

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
      rownames_to_column("cell_group")
  })


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


compute_vhat = function(ccm) {
  if (model(ccm)$d > 0) {
    ## self$fisher$mat : Fisher Information matrix I_n(\Theta) = n * I(\Theta)
    ## safe inversion using Matrix::solve and Matrix::diag and error handling

    vcov_mat = vcov(model(ccm))

    vhat <- matrix(0, nrow = nrow(vcov_mat), ncol = ncol(vcov_mat))

    #dimnames(vhat) <- dimnames(vcov_mat)
    safe_rows = safe_cols = Matrix::rowSums(abs(vcov_mat)) > 0
    vcov_mat = vcov_mat[safe_rows, safe_cols]

    out <- tryCatch(Matrix::solve(vcov_mat),
                    error = function(e) {e})
    row.names(out) = colnames(out) = names(safe_rows[safe_rows])
    if (is(out, "error")) {
      warning(paste("Inversion of the Fisher information matrix failed with following error message:",
                    out$message,
                    "Returning NA",
                    sep = "\n"))
      vhat <- matrix(NA, nrow = model(ccm)$p, ncol = model(ccm)$d)
    } else {
      row.names(out) = colnames(out) = names(safe_rows[safe_rows])
      row.names(vhat) = colnames(vhat) = row.names(vcov(model(ccm)))
      vhat[safe_rows, safe_cols] = as.numeric(out) #as.numeric(out) #%>% sqrt %>% matrix(nrow = self$d) %>% t()
    }
    #dimnames(vhat) <- dimnames(vcov_mat)
  } else {
    vhat <- NULL
  }
  vhat
}

#' Predict cell type abundances given a PLN model and a set of inputs for its covariates
#'
#' @param newdata needs to be suitable input to pln_model
#' @importFrom tibble tibble
#' @export
estimate_abundances <- function(ccm, newdata, min_log_abund=-5){

  # check that all terms in new data have been specified
  missing_terms = setdiff(names(ccm@model_aux$xlevels), names(newdata))

  if (length(missing_terms) > 1) {
    missing_terms = paste(missing_terms,collapse = ", ")
  }

  assertthat::assert_that(
    tryCatch(expr = length(missing_terms) == 0,
             error = function(e) FALSE),
    msg = paste0(missing_terms, " missing from newdata columns"))



  #stopifnot(nrow(newdata) == 1)
  newdata$Offset = 1

  base_X <- model.matrix(terms(ccm@model_aux[["model_frame"]]), newdata,
                         xlev = ccm@model_aux[["xlevels"]])


  #base_X <- model.matrix(formula(ccm@model_formula_str)[-2], newdata,
  #                  xlev = ccm@model_aux[["xlevels"]])
  X = Matrix::bdiag(rep.int(list(base_X), model(ccm)$p))

  # if it doesn't exist, use orig computation
  if (is.na(ccm@bootstrapped_vhat)) {
    v_hat = compute_vhat(ccm) #vcov(ccm@best_model)
  } else {
    v_hat = ccm@bootstrapped_vhat
  }


  se_fit = sqrt(diag(as.matrix(X %*% v_hat %*% Matrix::t(X)))) #/ sqrt(model(ccm)$n)

  pred_out = my_plnnetwork_predict(ccm, newdata=newdata)
  #pred_out = max(pred_out, -5)
  #log_abund = pred_out[1,]
  log_abund = as.numeric(pred_out)
  log_abund_sd = sqrt(diag(coef(model(ccm), type="covariance")))
  log_abund_se = se_fit

  log_abund[log_abund < min_log_abund] = min_log_abund
  #max_log_abundances = log(matrixStats::colMaxs(pln_model$fitted))
  #min_log_abundances = log(matrixStats::colMins(pln_model$fitted))
  #percent_max = 100 * (exp(log_abund)/exp(max_log_abundances))
  #percent_range = 100 * (exp(log_abund) - exp(min_log_abundances)) / (exp(max_log_abundances) - exp(min_log_abundances))
  pred_out_tbl = tibble::tibble(cell_group=colnames(pred_out),
                        log_abund,
                        log_abund_se)
  pred_out_tbl = left_join(pred_out_tbl,  tibble::tibble(cell_group=names(log_abund_sd), log_abund_sd), by=c("cell_group"))
  #max_log_abundances,
  #min_log_abundances,
  #percent_max,
  #percent_range)
  newdata$Offset = NULL
  pred_out_tbl = cbind(newdata, pred_out_tbl)
  return(pred_out_tbl)
}


# To do : need better error message for when you are missing a column that needs to be specified for the model
#' Predict cell type abundances given a PLN model over a range of time or other interval
#'
#' @importFrom tibble tibble
#' @export
estimate_abundances_over_interval <- function(ccm, start, stop, interval_col="timepoint", interval_step=2, ...) {

  timepoint_pred_df = tibble(IV= seq(start, stop, interval_step), ...)
  colnames(timepoint_pred_df)[1] = interval_col

  time_interval_pred_helper = function(tp, ...){
    tp_tbl = tibble(IV=tp, ...)
    colnames(tp_tbl)[1] = interval_col
    estimate_abundances(ccm, tp_tbl)
  }

  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      .f = time_interval_pred_helper, NA_real_),
      .x = !!sym(interval_col),
      ...)) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  return(timepoint_pred_df)
}


#' Compare two estimates of cell abundances from a Hooke model
#' @param ccm A cell_count_model
#' @param cond_x An estimate from estimate_abundances()
#' @param cond_y An estimate from estimate_abundances()
#' @return A table contrasting cond_x and cond_y (interpret as Y/X)
#' @importFrom dplyr full_join
#' @export
compare_abundances <- function(ccm, cond_x, cond_y, method = "BH"){
  contrast_tbl = dplyr::full_join(cond_x, cond_y, suffix = c("_x", "_y"), by="cell_group")
  contrast_tbl = contrast_tbl %>% dplyr::mutate(delta_log_abund = log_abund_y - log_abund_x,
                                                delta_p_value = pnorm(abs(delta_log_abund), sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2), lower.tail=FALSE),
                                                delta_q_value = p.adjust(delta_p_value, method = method))
  return(contrast_tbl)
}

correlate_abundance_changes <- function(pln_model, cond_b_vs_a_tbl){

  cov_graph <- return_igraph(pln_model)
  cov_edges <- igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight != 0.00)
  change_corr_tbl = cov_edges %>% dplyr::select(from, to, weight) %>% dplyr::rename(pcor = weight)

  #corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
  change_corr_tbl = dplyr::left_join(change_corr_tbl, cond_b_vs_a_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
  #dplyr::rename(log_abund_x,
  #              to_delta_log_abund = delta_log_abund)
  change_corr_tbl = dplyr::left_join(change_corr_tbl, cond_b_vs_a_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%
  #  dplyr::rename(from_delta_log_abund = delta_log_abund)
  return(change_corr_tbl)
}
