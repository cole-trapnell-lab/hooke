#' Calculate the shortest path between two points
#' @param edges data frame of edges with edge weights
#' @param from
#' @param to
#' @return data frame containing the shortest path
#' @export
calc_shortest_path <- function(G, from, to) {

  mf = igraph::shortest_paths(G, from = from, to = to, weights = igraph::E(G)$weight, output="epath")

  if (is.null(mf$epath) | length(mf$epath[[1]]) == 0) {
    stop(paste("No path from", from, " to ", to))
    #return(data.frame("from"=from, "to"=to, weight = -1, distance_from_root = -1))
  }

  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph)

  return(new_pos_edge_coords_df)
}


#' @noRd
get_neg_dir_edges <- function(ccm, cond_b_vs_a_tbl, q_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_q_value < q_value_threshold &
             from_delta_q_value < q_value_threshold)

}

#' @noRd
get_positive_edges <- function(ccm, cond_b_vs_a_tbl, q_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
             to_delta_q_value < q_value_threshold &
             from_delta_q_value < q_value_threshold)

}
