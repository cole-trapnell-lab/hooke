
# NOTE: This ignores offsets. DO NOT USE TO RECOVER FITTED VALUES FROM THE ORIGINAL DATA
#' @noRd
my_plnnetwork_predict <- function (ccm, newdata, type = c("link", "response"), envir = parent.frame())
{
  type = match.arg(type)
  X <- model.matrix(terms(ccm@model_aux[["model_frame"]]), newdata,
                         xlev = ccm@model_aux[["xlevels"]])
  #O <- model.offset(ccm@model_aux[["model_frame"]])
  EZ <- tcrossprod(X, t(model(ccm)$model_par$B))
  #if (!is.null(O))
  #  EZ <- EZ + O
  EZ <- sweep(EZ, 2, 0.5 * diag(model(ccm)$model_par$Sigma), "+")
  colnames(EZ) <- colnames(model(ccm)$model_par$Sigma)
  results <- switch(type, link = EZ, response = exp(EZ))
  attr(results, "type") <- type
  results
}


#' Predict cell type abundances given a PLN model and a set of inputs for its covariates
#'
#' @param ccm A cell_count_model.
#' @param newdata tibble A tibble of variables used for the prediction.
#' @param min_log_abund numeric Minimum log abundance value.
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances <- function(ccm, newdata, min_log_abund=-5) {

  assertthat::assert_that(is(ccm, 'cell_count_model'))
  assertthat::assert_that(tibble::is_tibble(newdata))
  assertthat::assert_that(is.numeric(min_log_abund))

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

  model_terms = terms(ccm@model_aux[["model_frame"]])
  base_X <- Matrix::sparse.model.matrix(model_terms, newdata,
                         xlev = ccm@model_aux[["xlevels"]])

  #base_X <- model.matrix(formula(ccm@model_formula_str)[-2], newdata,
  #                  xlev = ccm@model_aux[["xlevels"]])
  X = Matrix::bdiag(rep.int(list(base_X), model(ccm)$p))

  # if it doesn't exist, use orig computation
  # if (is.na(ccm@bootstrapped_vhat)[[1]]) {
  #   v_hat = compute_vhat(ccm)
  # } else {
  #   v_hat = ccm@bootstrapped_vhat
  # }

  vhat_coef <- coef(model(ccm), type="main")

  vcov_type <- grep('vcov', names(attributes(vhat_coef)), value=TRUE)
  v_hat <- attr(vhat_coef, vcov_type)
  v_hat_method <- ccm@vhat_method   

  if (v_hat_method == "wald") {
    se_fit = sqrt(diag(as.matrix(X %*% v_hat %*% Matrix::t(X)))) / sqrt(model(ccm)$n)
  } else {
    se_fit = sqrt(diag(as.matrix(X %*% v_hat %*% Matrix::t(X))))
  }

  pred_out = my_plnnetwork_predict(ccm, newdata=newdata)
  #pred_out = max(pred_out, -5)
  #log_abund = pred_out[1,]
  log_abund = as.numeric(pred_out)
  log_abund_sd = sqrt(diag(coef(model(ccm), type="covariance")))
  log_abund_se = se_fit

  below_thresh = log_abund < min_log_abund
  log_abund[below_thresh] = min_log_abund
  log_abund_se[below_thresh] = 0

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
  pred_out_tbl <- tibble::tibble(pred_out_tbl)

  return(pred_out_tbl)
}


# To do : need better error message for when you are missing a column that needs to be specified for the model
#' Predict cell type abundances given a PLN model over a range of time or other interval
#' @param ccm A cell_count_model.
#' @param interval_start numeric Interval start value.
#' @param interval_stop numeric Interval stop value.
#' @param interval_var character Interval values are taken from the interval_var data. Default is "timepoint".
#' @param interval_step numeric Interval size. Default is 2.
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances_over_interval <- function(ccm, interval_start, interval_stop, interval_var="timepoint", interval_step=2, ...) {

  assertthat::assert_that(is(ccm, 'cell_count_model'))
  assertthat::assert_that(is.numeric(min_log_abund))
  assertthat::assert_that(is.numeric(interval_start))
  assertthat::assert_that(is.numeric(interval_stop))
  assertthat::assert_that(interval_stop >= interval_start)
  assertthat::assert_that(is.numeric(interval_step))

  assertthat::assert_that(interval_var %in% attr(terms(ccm@model_aux[['model_frame']]), 'term.labels'))

  timepoint_pred_df = tibble(IV= seq(interval_start, interval_stop, interval_step), ...)
  colnames(timepoint_pred_df)[1] = interval_var

  time_interval_pred_helper = function(tp, ...){
    tp_tbl = tibble(IV=tp, ...)
    colnames(tp_tbl)[1] = interval_var
    estimate_abundances(ccm, tp_tbl)
  }

  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      .f = time_interval_pred_helper, NA_real_),
      .x = !!sym(interval_var),
      ...)) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  return(timepoint_pred_df)
}


#' Compare two estimates of cell abundances from a Hooke model.
#'
#' @param ccm A cell_count_model.
#' @param cond_x tibble A cell type abundance estimate from estimate_abundances().
#' @param cond_y tibble A cell type abundance estimate from estimate from estimate_abundances().
#' @param method string A method for correcting P-value multiple comparisons.
#'    This can be "BH" (Benjamini & Hochberg), "bonferroni" (Bonferroni),
#'    "hochberg" (Hochberg), "hommel", (Hommel), or "BYH" (Benjamini & Yekutieli).
#' @return tibble A table contrasting cond_x and cond_y (interpret as Y/X).
#' @importFrom dplyr full_join
#' @export
compare_abundances <- function(ccm, cond_x, cond_y, method = c("BH","bonferroni", "hochberg", "hommel", "BY")){

  assertthat::assert_that(is(ccm, 'cell_count_model'))
  assertthat::assert_that(tibble::is_tibble(cond_x))
  assertthat::assert_that(tibble::is_tibble(cond_y))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument method must be one of "size_factors",',
                '"BH", "bonferroni", "hochberg", "hommel", or "BY".'))
  method <- match.arg(method)

  contrast_tbl = dplyr::full_join(cond_x, cond_y, suffix = c("_x", "_y"), by="cell_group")

  # num samples
  n = nrow(model(ccm)$fitted)
  # num parameters
  k = length(rownames(coef(ccm@best_full_model)))
  df.r = n - k - 1

  contrast_tbl = contrast_tbl %>% dplyr::mutate(delta_log_abund = log_abund_y - log_abund_x,
                                                tvalue = delta_log_abund/(sqrt(log_abund_se_y^2 + log_abund_se_x^2)),
                                                delta_p_value = 2 * pt(-abs(tvalue), df.r),
                                                # delta_p_value = pnorm(abs(delta_log_abund), sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2), lower.tail=FALSE),
                                                delta_q_value = p.adjust(delta_p_value, method = method)) %>% select(-tvalue)
  return(contrast_tbl)
}


#' @noRd
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
