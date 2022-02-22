#' Calculate max flow between two points
calc_max_flow <- function(edges, source, target) {

  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  mf = max_flow(G, source = source, target = target, capacity = igraph::E(G)$pcor)
  igraph::E(G)$flow <- mf$flow
  igraph::V(G)$pass <- igraph::strength(G,mode="in",weights=mf$flow)

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(G)

  # if flow values are neg, reverse them
  switch_dir_df = new_pos_edge_coords_df %>% filter(flow < 0) %>%
    dplyr::rename("from"="to",
                  "umap_from_1"="umap_to_1",
                  "umap_from_2"="umap_to_2",
                  "to"="from",
                  "umap_to_1"="umap_from_1",
                  "umap_to_2"="umap_from_2") %>%
    mutate(flow = -flow)
  max_flow_edge_df = rbind(filter(new_pos_edge_coords_df, flow >=0),
                           switch_dir_df) %>%
    filter(flow > 0)
  return(max_flow_edge_df)

}

#' Calculate a minimum spanning tree
calc_mst <- function(edges, weight = "pcor") {
  G <- igraph::graph_from_data_frame(edges, directed=FALSE)

  if (weight == "pcor") {
    G_mst <- mst(G, weights=igraph::E(G)$pcor)
  }
  if (weight == "weight") {
    G_mst <- mst(G, weights=igraph::E(G)$weight)
  }

  mst_df <- igraph::as_data_frame(G_mst)
  return(mst_df)
}


#' Calculate the shortest path between two points
#' @param edges data frame of edges with edge weights
#' @param from
#' @param to
#' @import igraph
#' @timport
#' @return data frame containing the shortest path
#' @export
calc_shortest_path <- function(edges, from, to) {
  edges = edges %>% dplyr::select(from, to, x, y, z, x_z, y_z, z_z, weight)

  # if from or to is not in edges, return NA
  if (!from %in% edges$from | !to %in% edges$to) {
    return(data.frame("from"=from, "to"=to, weight = -1, distance_from_root = -1))

  }

  G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))

  mf = igraph::shortest_paths(G, from = from, to = to, weights = igraph::E(G)$weight, output="epath")

  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph)

  return(new_pos_edge_coords_df)
}


#'
#' @param edges
#' @import tidygraph
#'
distance_to_root <- function(edges) {

  g = edges %>% select(from,to) %>% igraph::graph_from_data_frame()
  distances = igraph::distances(g) %>%
    as.data.frame() %>%
    rownames_to_column("to") %>%
    pivot_longer(-c("to"), names_to = "from", values_to = "distance_from_root")

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

#' calculates the shortest weighted path between two points
#' @param weighted_edge_df
#' @param source
#' @param target
#'
get_shortest_path <- function(from, to, weighted_edges) {
  # print(paste0(from, "-",to))
  shortest_path_df = calc_shortest_path(weighted_edges, from, to) %>%
    select(from, to, weight) %>%
    distance_to_root() %>%
    arrange(distance_from_root)


  return(shortest_path_df)
}

#'
#' @param shortest_path
get_path_order <- function(shortest_path) {
  state_order = shortest_path %>%
    mutate(order = paste(from, to, sep="_")) %>%
    pull(order) %>%
    stringr::str_split("_") %>%
    unlist() %>% unique
  return(state_order)
}

# select adjacent states
select_states <- function(ordered_path, start , n = 3) {
  i = which(ordered_path == start)
  j = i + n - 1
  return(ordered_path[i:j])
}

#'
#' @param ccm
#' @param cond_a_vs_tbl
#' @param p_value_threshold
#'
get_path <- function(ccm, cond_b_vs_a_tbl, p_value_threshold = 1.0) {

  pos_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
           to_delta_p_value < p_value_threshold &
           from_delta_p_value < p_value_threshold)

  assertthat::assert_that(
    tryCatch(expr = nrow(pos_edges) != 0,
             error = function(e) FALSE),
    msg = "no significant positive edges found")

  weighted_edges = get_weighted_edges(ccm, pos_edges)

  neg_rec_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold)

  assertthat::assert_that(
    tryCatch(expr = nrow(neg_rec_edges) != 0,
             error = function(e) FALSE),
    msg = "no significant negative reciprocal edges found")

  edge_path = neg_rec_edges %>%
    dplyr::mutate(shortest_path = purrr::map2(.f =
                                                purrr::possibly(get_shortest_path, NA_real_),
                                              .x = from, .y = to,
                                              weighted_edges)) %>%
    select(shortest_path) %>%
    tidyr::unnest(shortest_path) %>%
    select(-weight) %>%
    distinct()

  return(edge_path)
}

#' Find degs either along an entire path or between adjacent pairs
#' @param ccm a cell count model object
#' @param cond_b_vs_a_tbl
#' @param p_value_threshold p_value threshold for inclusion of edges
#' @param adj_pairs If TRUE, does the DEG testing between adjacent pairs. If FALSE, includes all states in the deg testing.
find_deg_path = function(ccm,
                         cond_b_vs_a_tbl,
                         p_value_threshold = 1.0,
                         adj_pairs = FALSE,
                         ...) {

  # 1. find all the negative reciprocal edges
  neg_rec_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold)

  assertthat::assert_that(
    tryCatch(expr = nrow(neg_rec_edges) != 0,
             error = function(e) FALSE),
    msg = "no significant negative reciprocal edges found")

  pos_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold)

  assertthat::assert_that(
    tryCatch(expr = nrow(pos_edges) != 0,
             error = function(e) FALSE),
    msg = "no significant positive edges found")

  weighted_edges = get_weighted_edges(ccm, pos_edges)

  # 2. find the shortest path between those edges
  neg_rec_edges = neg_rec_edges %>%
    dplyr::mutate(shortest_path = purrr::map2(.f = purrr::possibly(get_shortest_path, NA_real_),
                                              .x = from, .y = to,
                                              weighted_edges)) %>%
    dplyr::rename("neg_from" = from, "neg_to" = to) %>%
    select(neg_from, neg_to, shortest_path)


  # two options
  if (adj_pairs) {

    # 3. find adjacent pairs along the path
    neg_rec_edges = neg_rec_edges %>%
      select(shortest_path) %>%
      tidyr::unnest(shortest_path) %>%
      select(-weight) %>%
      distinct() %>%
      dplyr::mutate(states = purrr::map(.f = collect_transition_states,
                                        .x = from,
                                        .y = to))
  } else {

    # 3. collect the states along each path

    neg_rec_edges = neg_rec_edges %>%
      dplyr::mutate(states = purrr::map(.f = collect_between_transition_states,
                                        .x = shortest_path))

  }

  # 4. create pseudobulk cds and run DEG analysis
  pb_cds = pseudobulk(ccm@ccs, ...)

  neg_rec_edges = neg_rec_edges  %>%
    dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_),
                                    .x = states,
                                    pb_cds))

  return(neg_rec_edges)
}


get_neg_dir_edges <- function(ccm, cond_b_vs_a_tbl, q_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_q_value < q_value_threshold &
             from_delta_q_value < q_value_threshold)

}

get_positive_edges <- function(ccm, cond_b_vs_a_tbl, q_value_threshold = 1.0) {
  hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
             to_delta_q_value < q_value_threshold &
             from_delta_q_value < q_value_threshold)

}

