#' Get the selected PLNnetwork from a cell_count_model object.
#'
#' @param ccm A cell_count_model object.
#' @return An updated cell_count_model object
#' @export
#' @examples
#' \dontrun{
#' model(ccm)
#' }
model <- function( ccm ) {
  stopifnot( methods::is( ccm, "cell_count_model" ) )
  ccm@best_model
}

#' Get the underlying cell_data_set object from a cell_count_model.
#'
#' @param ccm A cell_count_model object.
#' @return A cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' cds(ccm)
#' }
cds <- function( ccm ) {
  stopifnot( methods::is( ccm, "cell_count_model" ) )
  ccm@ccs@cds
}

#' Get the centroids of cell groups in UMAP/PCA space.
#'
#' @param ccs A cell_count_set object.
#' @return A data frame of centroid coordinates
#' @export
centroids <- function(ccs, reduction_method="UMAP") {
  # TODO: checks that reduction_method is valid, exists in cds, etc.
  coord_matrix = reducedDims(ccs@cds)[[reduction_method]] %>% as.data.frame
  grp_assign = ccs@metadata[["cell_group_assignments"]]
  grp_assign = grp_assign %>% dplyr::select(cell_group)
  coord_matrix = cbind(grp_assign, coord_matrix[row.names(grp_assign),])
  centroid_coords = aggregate(.~cell_group, data=coord_matrix, FUN=mean) 
  colnames(centroid_coords)[-1] = paste0(tolower(reduction_method), "_", 1:(length(colnames(centroid_coords))-1))
  return (centroid_coords)
}
