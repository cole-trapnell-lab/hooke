#' adds umap coords to a data frame
#' @import dplyr
add_umap_coords <- function(df, umap_centers) {

  from_df = left_join(df, umap_centers, by = c("from"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_from_", .x))

  to_df = left_join(df, umap_centers, by = c("to"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_to_", .x))

  return(full_join(from_df, to_df))

}

#' return the complementing edges
#' @importFrom igraph graph_from_data_frame
blacklist <- function(edges) {
  igraph::graph_from_data_frame(edges) %>%
    igraph::complementer() %>%
    igraph::as_data_frame()
}

#' returns edges based on pcor values
get_pcor_edges <- function(ccm) {
  ccm@best_model$latent_network() %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-c("from"), names_to = "to") %>%
    filter(from!=to, value!=0)
}

#' returns pairwise distances
get_distances <- function(ccs, method="euclidean", matrix=T) {
  umap_centers = centroids(ccs)
  row.names(umap_centers) <- umap_centers$cell_group
  
  dist_matrix = dist(umap_centers[,-1], method = method, upper=T, diag = T)
  
  if (matrix) {
    return(dist_matrix)
  } else {
    dist_df =  dist_matrix %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("from") %>%
      pivot_longer(-from, names_to = "to", values_to = "dist")
    return(dist_df)
  }
}




