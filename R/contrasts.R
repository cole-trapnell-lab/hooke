
#' Predict cell type abundances given a PLN model and a set of inputs for its covariates
#'
#' @param newdata needs to be suitable input to pln_model
#' @importFrom tibble tibble
#' @export
estimate_abundances <- function(ccm, newdata, min_log_abund=-5){
  #stopifnot(nrow(newdata) == 1)
  newdata$Offset = 1
  pred_out = predict(model(ccm), newdata=newdata)
  #pred_out = max(pred_out, -5)
  log_abund = pred_out[1,]
  log_abund[log_abund < min_log_abund] = min_log_abund
  #max_log_abundances = log(matrixStats::colMaxs(pln_model$fitted))
  #min_log_abundances = log(matrixStats::colMins(pln_model$fitted))
  #percent_max = 100 * (exp(log_abund)/exp(max_log_abundances))
  #percent_range = 100 * (exp(log_abund) - exp(min_log_abundances)) / (exp(max_log_abundances) - exp(min_log_abundances))
  pred_out_tbl = tibble::tibble(cell_group=colnames(pred_out),
                        log_abund)
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
compare_abundances <- function(ccm, cond_x, cond_y){
  contrast_tbl = dplyr::full_join(cond_x, cond_y, suffix = c("_x", "_y"), by="cell_group")
  contrast_tbl = contrast_tbl %>% dplyr::mutate(delta_log_abund = log_abund_y - log_abund_x)
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
