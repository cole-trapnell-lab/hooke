#' Calculate the shortest path between two points
#' @param edges data frame of edges with edge weights
#' @param from
#' @param to
#' @return data frame containing the shortest path
#' @export
calc_shortest_path <- function(G, from, to) {
  #edges = edges %>% dplyr::select(from, to, weight) %>% distinct()

  # if from or to is not in edges, return NA
  #if (!from %in% edges$from | !to %in% edges$to) {
  #  return(data.frame("from"=from, "to"=to, weight = -1, distance_from_root = -1))
  #
  # }

  #G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))

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


#'
#' @param edges
#' @import tidygraph
#'
#' @noRd
distance_to_root <- function(edges) {

  g = edges %>% select(from,to) %>% distinct() %>% igraph::graph_from_data_frame()
  distances = igraph::distances(g) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("to") %>%
    tidyr::pivot_longer(-c("to"), names_to = "from", values_to = "distance_from_root")

  node_info = tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    mutate(is_leaf = tidygraph::node_is_leaf(),
           is_root = tidygraph::node_is_root()) %>%
    as.data.frame()

  distance_df = inner_join(distances,
                           filter(node_info, is_root),
                           by = c("from"="name")) %>%
    select(-c(is_leaf,is_root, from))

  left_join(edges, distance_df, by = c("to"))
}


#' weights edges based on partial correlation,
#' change in abundance, and transcriptomic distance
#' @param ccm
#' @param edges
#' @param alpha
#' @param beta
#' @param gamma
#' @param sum_weights
#' @noRd
get_weighted_edges <- function(ccm,
                               edges,
                               alpha= 1,
                               beta = 1,
                               gamma = 1,
                               sum_weights = T) {

  # get weighted path
  dist_df = get_distances(ccm@ccs, matrix = F)

  weighted_edge_df = left_join(edges, dist_df, by=c("from", "to")) %>%
    mutate(x = abs(from_delta_log_abund-to_delta_log_abund),
           y = 1/pcor,
           z = dist) %>%
    mutate(x_z = (x-min(x))/(max(x)-min(x)),
           y_z = (y-min(y))/(max(y)-min(y)),
           z_z = (z-min(z))/(max(z)-min(z)))

  if (sum_weights) {
    weighted_edge_df = weighted_edge_df %>% mutate(weight = alpha*x_z + beta*y_z + gamma*z_z)
  } else {
    weighted_edge_df = weighted_edge_df %>% mutate(weight = alpha*x_z * beta*y_z * gamma*z_z)
  }
  return(weighted_edge_df)

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
