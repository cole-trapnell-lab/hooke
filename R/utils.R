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


get_neg_edges <- function(ccm, cond_b_vs_a_tbl,  p_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(a549_ccm, cond_0_vs_10000_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold) 
  
}

get_pos_edges <- function(ccm, cond_b_vs_a_tbl,  p_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(a549_ccm, cond_0_vs_10000_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold)
  
}

add_covariate <- function(ccs, pb_cds, covariate) {
  assertthat::assert_that(
    tryCatch(expr = covariate %in% colnames(colData(ccs@cds)), 
             error = function(e) FALSE), 
    msg = "covariate not in colnames")
  
  assertthat::assert_that(
    tryCatch(expr = !covariate %in% colnames(colData(pb_cds)), 
             error = function(e) FALSE), 
    msg = "covariate already in colnames")
  
  group_to_covariate = colData(ccs@cds) %>% 
    as.data.frame %>%
    select(group_id, all_of(covariate)) %>% 
    distinct()
  
  pb_coldata = colData(pb_cds) %>% 
    as.data.frame %>% 
    left_join(group_to_covariate, by = "group_id")
  
  colData(pb_cds)[[covariate]] =  pb_coldata[[covariate]]
  
  return(pb_cds)
}




