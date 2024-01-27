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
  coord_matrix = ccs@cds_reduced_dims[[reduction_method]] %>% as.data.frame #reducedDim(ccs@cds, reduction_method) %>% as.data.frame

  if (is.null(switch_group)==FALSE) {
    grp_assign = ccs@cds_coldata %>% as.data.frame %>% dplyr::select(!!sym(switch_group))
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

#' return a matrix of normalized counts
#' @param ccs cell_count_set
#' @param round TRUE if counts are rounded
#' @export
get_norm_counts <- function(ccs, round = FALSE) {
  if (round) {
    norm_counts = round(t((t(as.matrix(counts(ccs)))/size_factors(ccs))))
  } else {
    norm_counts = t((t(as.matrix(counts(ccs)))/size_factors(ccs)))
  }

  return(norm_counts)

}

get_count_df <- function(ccs, round=F, norm=F) {

  if (norm) {
    count_mat = get_norm_counts(ccs, round=round)
  } else {
    count_mat = counts(ccs)
  }

  count_df = count_mat %>%
    as.matrix %>%
    as.data.frame %>%
    rownames_to_column("cell_group") %>%
    pivot_longer(-cell_group, names_to = "sample", values_to = "count")

}

#' subset ccs by cell groups
#' #' @param ccs
#' #' @param cell_groups
#' #' @export
subset_ccs = function(ccs, ...) {

  filtered_cds_coldata = ccs@cds_coldata %>% as.data.frame %>% filter(...)
  group_ids = filtered_cds_coldata %>% pull(group_id)

  cell_group_assignments = ccs@metadata$cell_group_assignments[ccs@metadata$cell_group_assignments$group_id %in% group_ids,]
  samples = intersect( colnames(ccs), unique(cell_group_assignments$sample))
  cell_groups = intersect(rownames(ccs), unique(cell_group_assignments$cell_group))
  sub_ccs = ccs[cell_groups, samples]
  sub_ccs@metadata$cell_group_assignments = cell_group_assignments
  sub_ccs@cds = sub_ccs@cds[, colData(sub_ccs@cds)$cell %in% filtered_cds_coldata$cell]

  return(sub_ccs)

}



#' #' subset ccs by cell groups
#' #' @param ccs
#' #' @param cell_groups
#' #' @export
#' subset_ccs = function(ccs, cell_groups) {
#'
#'   # make sure that cell groups are in the ccs
#'
#'   assertthat::assert_that(
#'     tryCatch(expr = all(cell_groups %in% rownames(ccs)),
#'              error = function(e) FALSE),
#'     msg = paste0(setdiff(cell_groups, rownames(ccs)) %>% paste0(collapse = ", "),
#'                  " are not found in the ccs"))
#'
#'   sub_ccs = ccs[cell_groups, ]
#'
#'   sub_ccs@metadata$cell_group_assignments = sub_ccs@metadata$cell_group_assignments[sub_ccs@metadata$cell_group_assignments$cell_group %in% cell_groups,]
#'
#'   if (is.null(nrow(colnames(ccs@cds@colData))) == FALSE) {
#'     sub_ccs@cds = sub_ccs@cds[,rownames(sub_ccs@metadata$cell_group_assignments)]
#'   }
#'
#'   return(sub_ccs)
#'
#' }
