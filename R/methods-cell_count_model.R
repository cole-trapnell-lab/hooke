#' Get the selected PLNnetwork from a cell_count_model object.
#'
#' @param ccm A cell_count_model object.
#' @return An updated cell_count_model object
#' @examples
#' \dontrun{
#' model(ccm)
#' }
#' @export
model <- function( ccm, model_to_return=c("full", "reduced") ) {
  model_to_return = match.arg(model_to_return)
  stopifnot( methods::is( ccm, "cell_count_model" ) )

  if (model_to_return == "full"){
    return(ccm@best_full_model)
  }
  if (model_to_return == "reduced"){
    return(ccm@best_reduced_model)
  }
}

#' Get the underlying cell_data_set object from a cell_count_model.
#'
#' @param ccm A cell_count_model object.
#' @return A cell_data_set object
#' @examples
#' \dontrun{
#' cds(ccm)
#' }
#' @export
cds <- function( ccm ) {
  stopifnot( methods::is( ccm, "cell_count_model" ) )
  ccm@ccs@cds
}

#' Get the centroids of cell groups in UMAP/PCA space.
#'
#' @param ccs A cell_count_set object.
#' @return A data frame of centroid coordinates
#' @export
centroids <- function(ccs, reduction_method="UMAP", switch_group = NULL) {
  # TODO: checks that reduction_method is valid, exists in cds, etc.
  coord_matrix = reducedDims(ccs@cds)[[reduction_method]] %>% as.data.frame

  if (is.null(switch_group)==FALSE) {
    grp_assign = ccs@cds@colData %>% as.data.frame %>% dplyr::select(!!sym(switch_group))
    colnames(grp_assign) = "cell_group"
  } else {
    grp_assign = ccs@metadata[["cell_group_assignments"]]
    grp_assign = grp_assign %>% dplyr::select(cell_group)
  }

  coord_matrix = cbind(grp_assign, coord_matrix[row.names(grp_assign),])
  centroid_coords = aggregate(.~cell_group, data=coord_matrix, FUN=mean)
  colnames(centroid_coords)[-1] = paste0(tolower(reduction_method), "_", 1:(length(colnames(centroid_coords))-1))
  return (centroid_coords)
}
