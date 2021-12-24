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
calc_mst <- function(edges) {
  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  G_mst <- mst(G, weights=igraph::E(G)$pcor)
  mst_df <- igraph::as_data_frame(G_mst)
  return(mst_df)
}


#' Calculate the shortest path between two points
#' @param edges data frame of edges with edge weights
#' @param source
#' @param target
#' @return data frame containing the shortest path
#' @export
calc_shortest_path <- function(edges, from, to) {
  edges = edges %>% dplyr::select(from, to, x, y, z, x_z, y_z, z_z, weight)
  G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))
  
  mf = igraph::shortest_paths(G, from = from, to = to, weights = igraph::E(G)$weight, output="epath")
  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])
  
  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph)
  
  return(new_pos_edge_coords_df)
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
get_shortest_path <- function(from, to, weighted_edges) {
  
  shortest_path_df = calc_shortest_path(weighted_edges, from, to) %>%
    select(from, to, weight)
  
  return(shortest_path_df)
}


         