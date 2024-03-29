#' Compute a pseudobulk expression matrix for a ccs
#' @param ccs a cell count set object
#' @param state_col column to aggregate expression. Defaults using the cell_group used in ccs construction.
#' @param collapse_samples boolean Whether to collapse sample groups into one.
#' @export
pseudobulk_ccs_for_states <- function(ccs, state_col=NULL, collapse_samples=FALSE){

  if (is.null(state_col)){
    cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]])
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = cell_group)
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, "cell_group", sep="_")) %>% dplyr::select(rowname, pseudobulk_id, cell_group)
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, cell_group) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
    #%>% select(rowname, cell_group)
  }else{
    cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]])
    cds_group_df = colData(ccs@cds) %>%
      as.data.frame %>% tibble::rownames_to_column() %>% dplyr::select(rowname, !!sym(state_col))
    cell_group_df = left_join(cell_group_df, cds_group_df, by=c("rowname"))
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = !!sym(state_col))
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, !!sym(state_col), sep="_")) %>% dplyr::select(rowname, pseudobulk_id, !!sym(state_col))
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, !!sym(state_col)) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
  }

  agg_expr_mat = monocle3::aggregate_gene_expression(ccs@cds,
                                                     cell_group_df=cell_group_df,
                                                     norm_method="size_only",
                                                     scale_agg_values = FALSE,
                                                     pseudocount=0,
                                                     cell_agg_fun="mean")

  agg_expr_mat = agg_expr_mat[,agg_coldata$pseudobulk_id]

  row.names(agg_coldata) = agg_coldata$pseudobulk_id
  agg_coldata = agg_coldata[colnames(agg_expr_mat),]

  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  return(pseudobulk_cds)
}

#' add metadata to pb_cds from cds
#' @param ccs a cell count set object
#' @param pb_cds a pseudobulked cds
#' @param covariate the column name in colData(ccs) to add to the pb_cds
#' @export
add_covariate <- function(ccs, pb_cds, covariate) {

  assertthat::assert_that(
    tryCatch(expr = covariate %in% colnames(colData(ccs@cds)),
             error = function(e) FALSE),
    msg = paste0(covariate, " not in colnames"))

  assertthat::assert_that(
    tryCatch(expr = !covariate %in% colnames(colData(pb_cds)),
             error = function(e) FALSE),
    msg = paste0(covariate," already in colnames"))

  group_to_covariate = colData(ccs@cds) %>%
    as.data.frame %>%
    select(group_id, all_of(covariate)) %>%
    distinct()

  pb_coldata = colData(pb_cds) %>%
    as.data.frame %>%
    mutate(group_id = gsub("_cell_group", "", pseudobulk_id)) %>%
    left_join(group_to_covariate, by = "group_id")

  colData(pb_cds)[[covariate]] =  pb_coldata[[covariate]]

  return(pb_cds)
}

