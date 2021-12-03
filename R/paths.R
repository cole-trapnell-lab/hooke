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
calc_shortest_path <- function(edges, source, target) {
  edges = edges %>% dplyr::select(from, to, x, y, z, x_z, y_z, z_z, weight)
  G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))

  mf = igraph::shortest_paths(G, from = source, to = target, weights = igraph::E(G)$weight, output="epath")
  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph)

  return(new_pos_edge_coords_df)
}


#' @param umap_centers
#' @param corr_edge_coords_umap_delta_abund
#' @param alpha coefficient to weight
#' @param beta coefficient to weight pcor
#' @param gamma coefficient to weight distance
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
get_pos_pcor_edges <- function(umap_centers,
                               corr_edge_coords_umap_delta_abund,
                               alpha= 1,
                               beta = 1,
                               gamma = 1,
                               sum_weights = F) {

  row.names(umap_centers) <- umap_centers$cell_group

  dist_df = dist(umap_centers[,-1], method = "euclidean", upper=T, diag = T) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "z")

  pos_edge_df = corr_edge_coords_umap_delta_abund %>%
    filter(pcor > 0)  %>%
    left_join(dist_df, by=c("from", "to")) %>%
    mutate(x = abs(from_delta_log_abund-to_delta_log_abund),
           y = 1/pcor) %>%
    mutate(x_z = (x-min(x))/(max(x)-min(x)),
           y_z = (y-min(y))/(max(y)-min(y)),
           z_z = (z-min(z))/(max(z)-min(z)))

  if (sum_weights) {
    pos_edge_df = pos_edge_df %>% mutate(weight = alpha*x_z + beta*y_z + gamma*z_z)
  } else {
    pos_edge_df = pos_edge_df %>% mutate(weight = alpha*x_z * beta*y_z * gamma*z_z)
  }

  return(pos_edge_df)

}




get_shortest_path <- function(ccm,
                              cond_a_vs_b_tbl,
                              source,
                              target,
                              alpha=1,
                              beta=1,
                              gamma=1,
                              log_abundance_thresh = -5,
                              sum_weights = T) {

  umap_centers = centroids(ccm@ccs)
  corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                              cond_a_vs_b_tbl,
                                                              log_abundance_thresh)

  target_edge_df = filter(corr_edge_coords_umap_delta_abund,
                          from == source,
                          to == target)

  pos_edge_df = get_pos_pcor_edges(umap_centers,
                                   corr_edge_coords_umap_delta_abund,
                                   alpha = alpha,
                                   beta = beta,
                                   gamma = gamma,
                                   sum_weights=sum_weights)

  path_df = calc_shortest_path(pos_edge_df, source, target) %>%
            select(from,to,weight)

  return(path_df)
}


collate_edges <- function(corr_edge_coords_umap_delta_abund, umap_centers) {

  neg_edges = corr_edge_coords_umap_delta_abund %>% filter(edge_type!="undirected")

  if (nrow(neg_edges) > 0) {
    green_edges = lapply(1:nrow(neg_edges), function(i) {
      source = neg_edges[i,]$from
      target = neg_edges[i,]$to
      get_pos_pcor_edges(umap_centers,
                         corr_edge_coords_umap_delta_abund,
                         alpha = 1,
                         beta = 1,
                         gamma = 1,
                         sum_weights=T) %>%
        calc_shortest_path(source, target) %>% select(from,to,weight)
    })

    green_edge_df = do.call(rbind,green_edges) %>%
      group_by(from,to) %>%
      summarise(weight=sum(weight), n=n()) %>%
      add_umap_coords(umap_centers)
    return(green_edge_df)
  }

}


collate_edges_all_combos <- function(ccm,
                                     timepoints) {

  timepoint_combos = combn(timepoints, 2) %>%
      t() %>%
      as.data.frame() %>%
      rename("t1"=V1, "t2"=V2)

  greenedges = lapply(1:nrow(timepoint_combos), function(i) {
    t1 = timepoint_combos[i,]$t1
    t2 = timepoint_combos[i,]$t2

    time_t1 = estimate_abundances(ccm, data.frame(timepoint=t1))
    time_t2 = estimate_abundances(ccm, data.frame(timepoint=t2))

    cond_t1_vs_t2_tbl = compare_abundances(ccm, time_t1, time_t2)

    corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                                cond_t1_vs_t2_tbl,
                                                                log_abundance_thresh=-5)



    green_edge_df = collate_green_edges(corr_edge_coords_umap_delta_abund, umap_centers) %>%
                    mutate(timepoint_x=t1, timepoint_y=t2)

    green_edge_df

  })

  greenedges_time_df = do.call(rbind, greenedges)

  total_green_edges = greenedges_time_df %>%
    group_by(from,to) %>%
    summarise(total_weight = sum(weight), total_n = sum(n))

  return(total_green_edges)
}


update_cofficients <- function(ccm,
                               edges,
                               covariate = "gene_target",
                               group = "tbx16-msgn1",
                               wt_ccm = NULL,
                               update = T) {

  model_coef = coef(ccm@best_model)
  coef_df = as.data.frame(model_coef)[[paste0(covariate, group)]]
  names(coef_df) = rownames(as.data.frame(model_coef))

  if (update == F) {
    return(coef_df)
  }

  coef_updated = coef_df

  for (i in 1:nrow(edges)) {

    x = edges[i,]$from
    y = edges[i,]$to

    if (is.null(wt_ccm)) {
      coef_updated[y] = coef_df[y] - coef_df[x]
    } else {
      coef_updated[y] = coef_df[y] - coef_df[x]*wt_pcor_matrix[x,y]
    }
  }
  return(coef_updated)
}


edges_across_conditions <- function(ccm, gene_targets) {

  greenedges = lapply(gene_targets, function(gt) {

    umap_centers = centroids(ccm@ccs)

    time_18 = estimate_abundances(ccm, data.frame(timepoint="18", gene_target = gt))
    time_24 = estimate_abundances(ccm, data.frame(timepoint="24", gene_target = gt))

    cond_a_vs_b_tbl = compare_abundances(ccm, time_18, time_24)
    corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                                cond_a_vs_b_tbl,
                                                                log_abundance_thresh=-5)

    green_edge_df = collate_edges(corr_edge_coords_umap_delta_abund,
                                  umap_centers)

  })
  greenedges_time_df = do.call(rbind, greenedges)

  total_green_edges = greenedges_time_df %>%
    group_by(from,to) %>%
    summarise(total_weight = sum(weight), total_n = sum(n))
  return(total_green_edges)
}

