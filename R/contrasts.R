
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
  #stopifnot(nrow(newdata) == 1)
  newdata$Offset = 1

  base_X <- model.matrix(terms(ccm@model_aux[["model_frame"]]), newdata,
                         xlev = ccm@model_aux[["xlevels"]])


  #base_X <- model.matrix(formula(ccm@model_formula_str)[-2], newdata,
  #                  xlev = ccm@model_aux[["xlevels"]])
  X = Matrix::bdiag(rep.int(list(base_X), model(ccm)$p))



  v_hat = compute_vhat(ccm) #vcov(ccm@best_model)
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
