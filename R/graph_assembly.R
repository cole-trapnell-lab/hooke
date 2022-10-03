#' Orient graph edges from a PLNnetwork using a contrast between conditions
#'
#'@param ccm A cell_count_model
#'@param cond_b_vs_a_tbl A contrast between two conditions as returned by compare_abundances()
collect_pln_graph_edges <- function(ccm,
                                    cond_b_vs_a_tbl,
                                    log_abundance_thresh = 1-5,
                                    model_for_pcors = "reduced"){
  pln_model = model(ccm, model_for_pcors)

  abundance_corr_tbl = correlate_abundance_changes(pln_model, cond_b_vs_a_tbl)

  abundance_corr_tbl = abundance_corr_tbl %>% dplyr::filter(
    (to_log_abund_x > log_abundance_thresh | to_log_abund_y > log_abundance_thresh) &  # Keep if the "to" node is above abundance thresh in at least one condition
      (from_log_abund_x > log_abundance_thresh | from_log_abund_y > log_abundance_thresh) # Keep if the "to" node is above abundance thresh in at least one condition
  )

  corr_edge_delta_abund = abundance_corr_tbl

  corr_edge_delta_abund = corr_edge_delta_abund %>% dplyr::mutate(scaled_weight = -pcor)

  corr_edge_delta_abund = corr_edge_delta_abund %>%
    dplyr::mutate(edge_type = dplyr::case_when(
      from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor < 0 ~ "directed_to_from",
      from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor < 0 ~ "directed_from_to",
      from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor < 0 ~ "undirected",
      TRUE ~ "hidden"
    ))

  # Fix the edges so all directed edges are directed_from_to:
  backwards_edges = corr_edge_delta_abund %>% dplyr::filter(edge_type == "directed_to_from")
  ok_edges = corr_edge_delta_abund %>% dplyr::filter(edge_type != "directed_to_from")

  be_colnames = colnames(backwards_edges)
  be_colnames = unlist(lapply(be_colnames, stringr::str_replace, "^to", "from"))
  from_cols = grepl("^from", colnames(backwards_edges))
  be_colnames[from_cols] = unlist(lapply(be_colnames[from_cols], stringr::str_replace, "^from", "to"))
  colnames(backwards_edges) = be_colnames
  if (nrow(backwards_edges) != 0) {
    backwards_edges$edge_type = "directed_from_to"
    backwards_edges = backwards_edges[colnames(ok_edges)]
  }

  corr_edge_coords_umap_delta_abund = rbind(ok_edges, backwards_edges)
  return (corr_edge_coords_umap_delta_abund)
}

return_igraph <- function(model, type = "partial_cor", remove.isolated=FALSE,
                          edge.color = c("#F8766D", "#00BFC4"),
                          output = "igraph"){

  net <- model$latent_network(type = type)

  if (output == "igraph") {

    G <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::V(G)$label <- colnames(net)

    ## Nice nodes
    V.deg <- igraph::degree(G)/sum(igraph::degree(G))
    igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
    igraph::V(G)$size <- V.deg * 100
    igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
    igraph::V(G)$frame.color <- NA

    ## Nice edges
    igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
    if (type == "support")
      igraph::E(G)$width <- abs(igraph::E(G)$weight)
    else
      igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)


    if (remove.isolated) {
      G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
    }
    G$layout <- igraph::layout_in_circle
  }
  if (output == "corrplot") {

    if (ncol(net) > 100)
      colnames(net) <- rownames(net) <- rep(" ", ncol(net))
    G <- net
    diag(net) <- 0
    corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
  }
  G
}


#' Extract a partitioned abstract graph from a Monocle cell_data_set object
#'
#' QUESTION: Do the cells need to be grouped by monocle cluster? Or can they
#' be grouped arbitrarily (e.g. by cell type? Would be good to generalize if
#' possible. For now, only works when cells are grouped by Monocle cluster
#'
#' @param cds A cell_data_set object. cluster_cells() must have been called.
#' @param reduction_method The coordinate space in which to build the graph
#' @export
get_paga_graph <- function(cds, reduction_method = "UMAP", partition_q_value=0.05) {

  cluster_result <- cds@clusters[[reduction_method]]$cluster_result
  #clusters <- ccm@ccs@cds@clusters[[reduction_method]]$clusters
  #partitions <- ccm@ccs@cds@clusters[[reduction_method]]$partitions

  #colData(cds)$cluster = clusters
  #colData(cds)$partition = partitions

  #coldata = colData(ccm@ccs@cds) %>% as_tibble() %>% dplyr::rename(embryo = Oligo) %>%
  #  mutate(timepoint=as.numeric(stringr::str_remove(timepoint, "hpf")))
  #coldata = ccm@ccs@metadata[["cell_group_assignments"]]

  cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                     cluster_result$optim_res,
                                                     partition_q_value, FALSE)


  #cluster_time = coldata %>% group_by(cluster) %>% summarise(avg_time = mean(as.numeric(timepoint)))
  #cluster_info = merge(cluster_maj, cluster_time, by="cluster") %>%
  #  mutate("label" = paste0(cell_type, ".", round(avg_time), ".", cluster)) %>%
  #  arrange(as.numeric(cluster))

  cluster_g <- cluster_graph_res$cluster_g
  cluster_g <- igraph::set_vertex_attr(cluster_g, "name", value = stringr::str_replace(igraph::V(cluster_g)$name, "cell_membership", ""))
  cluster_g
}

#' Get an initial graph for use as a whitelist in fitting a cell count model
#'
#' @param ccs A cell_count_model object.
#'
#' @export
initial_pcor_graph = function(ccs) {
  paga_edges = get_paga_graph(ccs@cds) %>% igraph::as_data_frame() %>% as_tibble()

  # filter out values that aren't in the cds anymore
  cell_groups = unique(colData(ccs@cds)[[ccs@info$cell_group]])

  paga_edges = paga_edges %>% filter(from %in% cell_groups, to %in% cell_groups)

  return(paga_edges)
}



#' Return a dataframe that describes which cell types are present in a time interval
#'
#' @export
get_extant_cell_types <- function(ccm,
                                  start,
                                  stop,
                                  interval_col="timepoint",
                                  interval_step = 2,
                                  log_abund_detection_thresh = -5,
                                  percent_max_threshold=0.0,
                                  pct_dynamic_range = 0.25,
                                  min_cell_range = 2,
                                  ...){

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step=interval_step, ...)

  if (is.null(log_abund_detection_thresh)) {
    abund_range = range(timepoint_pred_df$log_abund)
    dynamic_range = abund_range[2]-abund_range[1]
    log_abund_detection_thresh = abund_range[1] + pct_dynamic_range*dynamic_range
  }

  timepoint_pred_df = timepoint_pred_df %>%
    group_by(cell_group) %>%
    mutate(max_abundance = max(exp(log_abund)),
           percent_max_abund = exp(log_abund) / max_abundance,
           cell_type_range = max(log_abund)-(min(log_abund)),
           above_log_abund_thresh = (log_abund - 2*log_abund_se > log_abund_detection_thresh) | cell_type_range < min_cell_range,
           above_percent_max_thresh = percent_max_abund > percent_max_threshold,
           present_flag = ifelse(above_log_abund_thresh & above_percent_max_thresh, TRUE, NA)) %>%
    ungroup()

  longest_present_interval <- function(tps_df){
    tryCatch(
      {
        delta_t = as.numeric(tps_df[2,1] - tps_df[1,1])
        ts_la = ts(tps_df$present_flag,
                   start=min(tps_df[,1]),
                   #end=max(tps_df[,1]),
                   deltat=delta_t)
        longest_contig = na.contiguous(ts_la)
        return(tibble(longest_contig_start = start(longest_contig)[1], longest_contig_end = end(longest_contig)[1]))
      }, error = function(e) {
        return (tibble(longest_contig_start = NA, longest_contig_end = NA))
      }
    )
  }

  # undebug(longest_present_interval)
  nested_timepoints_df = timepoint_pred_df %>%
    select(cell_group, !!sym(interval_col), present_flag) %>%
    group_by(cell_group)

  nested_timepoints_df = nested_timepoints_df %>%
    group_modify(~ longest_present_interval(.x))
  #mutate(cg_ts = purrr:::map2(.f=purrr::possibly(longest_present_interval, NA_real_),
  #                            .x=!!sym(interval_col),
  #                            .y=present_flag))

  timepoint_pred_df = left_join(timepoint_pred_df, nested_timepoints_df)

  timepoint_pred_df = timepoint_pred_df %>%
    mutate(present_above_thresh = !!sym(interval_col) >= longest_contig_start & !!sym(interval_col) <= longest_contig_end,
           present_above_thresh = ifelse(is.na(present_above_thresh), FALSE, present_above_thresh))

  extant_cell_type_df = timepoint_pred_df %>%
    select(!!sym(interval_col),
           cell_group,
           log_abund,
           max_abundance,
           percent_max_abund,
           longest_contig_start,
           longest_contig_end,
           present_above_thresh)
  return(extant_cell_type_df)
}
#undebug(get_extant_cell_types)
#get_extant_cell_types(wt_ccm_wl, 72, 96) %>% filter(cell_group == "21")


weigh_edges_by_umap_dist <- function(ccm, edges) {

  # get weighted path
  dist_df = hooke:::get_distances(ccm@ccs, matrix = F)

  weighted_edge_df = left_join(edges, dist_df, by=c("from", "to")) %>%
    mutate(weight = dist)


  return(weighted_edge_df)

}

weigh_edges <- function(ccm, edges) {
  # get weighted path
  dist_df = hooke:::get_distances(ccm@ccs, matrix = F)


}

# Helper function
get_destinations_that_need_origins <- function(neg_rec_edge_df,
                                               current_time,
                                               extant_cell_type_df,
                                               interval_col){
  expanding_destinations = neg_rec_edge_df %>% select(cell_group=to, to_delta_log_abund) %>% filter(to_delta_log_abund > 0)
  extant_cell_type_df = extant_cell_type_df %>% filter(!!sym(interval_col) == current_time & present_above_thresh)
  inner_join(expanding_destinations, extant_cell_type_df, by="cell_group") %>% pull(cell_group) %>% unique()
}

# Helper function
get_neg_rec_edges_for_cells_that_need_origins <- function(neg_rec_edge_df, need_origins){
  return(neg_rec_edge_df %>% filter(to %in% need_origins))
}

add_cross_component_pathfinding_links = function(ccm, pathfinding_graph){

  pcor_graph = tibble(cell_group = row.names(counts(ccm@ccs)))
  pcor_graph = pcor_graph %>% select(cell_group) %>% tidyr::expand(cell_group, cell_group)
  colnames(pcor_graph) = c("from", "to")
  pcor_graph = pcor_graph %>% filter(from != to)

  nz_pcors = igraph::graph_from_adjacency_matrix(model(ccm, "reduced")$latent_network(type = "partial_cor"), mode = "undirected", weighted = TRUE, diag = FALSE) %>% igraph::as.directed()
  nz_pcor_graph_edges = igraph::as_data_frame(nz_pcors, what="edges") %>% dplyr::rename(pcor=weight) %>% select(from, to, pcor)#%>% dplyr::filter(pcor != 0.00)

  pcor_graph = dplyr::left_join(pcor_graph, nz_pcor_graph_edges)

  #pcor_graph = pcor_graph %>% tidyr::replace_na(list(pcor = 0))

  #pcor_graph = pcor_graph %>% left_join()
  pcor_graph = hooke:::weigh_edges_by_umap_dist(ccm, pcor_graph) %>% dplyr::rename(umap_dist=weight)

  #pcor_graph = pcor_graph %>% filter(pcor != 0)
  pcor_graph = pcor_graph %>% mutate(abs_pcor = abs(pcor))
  #pcor_vs_dist_model = VGAM::vglm(pcor ~ I(1/umap_dist), data=pcor_graph, family=VGAM::gamma2(zero=NULL, lmu="identitylink", lshape="loglink"))
  #mod_predict = predict(pcor_vs_dist_model, newdata=pcor_graph) %>% as.data.frame()
  #pcor_graph$model_fit = mod_predict$mu
  #pcor_graph$model_sd = sqrt((mod_predict$mu)^2 / exp(mod_predict$`loglink(shape)`))

  pcor_vs_dist_model = VGAM::vglm(pcor ~ umap_dist, data=pcor_graph, family=VGAM::gaussianff(zero=NULL, lmean="identitylink", lsd="loglink"))
  mod_predict = predict(pcor_vs_dist_model, newdata=pcor_graph) %>% as.data.frame()
  pcor_graph$model_fit = mod_predict$mean
  #pcor_graph$model_sd = sqrt((mod_predict$mean)^2 / exp(mod_predict$`loglink(sd)`))
  #pcor_graph$model_sd = sqrt((mod_predict$mean)^2 / mod_predict$sd)
  pcor_graph$model_sd = mod_predict$sd
  pcor_graph$model_sd = exp(mod_predict$`loglink(sd)`)


  pcor_graph = pcor_graph %>%
    mutate(model_fit_lower = model_fit - 2 * model_sd,
           #model_fit_lower = ifelse(model_fit_lower > 0, model_fit_lower, min(abs_pcor)),
           model_fit_upper = model_fit + 2 * model_sd)

  unexpectedly_strong_pcor_edges = pcor_graph %>%
    filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
    select(from, to, weight=umap_dist)
  uspe_graph = igraph::graph_from_data_frame(unexpectedly_strong_pcor_edges, directed=TRUE, vertices=data.frame(id=row.names(counts(ccm@ccs))))
  updated_pathfinding_graph = igraph::union(pathfinding_graph, uspe_graph)

  return(updated_pathfinding_graph)
}

#' Initialize a graph over which cells can transition
#'
#' @param ccm A cell count model
#' @param extant_cell_type_df A data frame describing when cell types are present, generated by get_extant_cell_types()
#' @param allow_links_between_components Whether cells in separate partitions of the UMAP can be linked by paths
#' @export
init_pathfinding_graph <- function(ccm, extant_cell_type_df, allow_links_between_components=FALSE, weigh_by_pcor=F){

  # There are a number of different ways we could set up this "pathfinding graph" but for now
  # Let's just use the PAGA (weighed by distance in UMAP space), subtracting edges between which
  # there is zero partial correlation in the nuisance cell count model.

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = data.frame(id=cell_groups)

  paga_graph = initial_pcor_graph(ccm@ccs) %>% igraph::graph_from_data_frame(directed = FALSE, vertices=node_metadata) %>% igraph::as.directed()
  cov_graph = hooke:::return_igraph(model(ccm, "reduced"))
  cov_graph_edges = igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::rename(pcor=weight) %>% dplyr::filter(pcor != 0.00)

  # FIXME: While we should consider this more sophisticated weighting function, let's test just using UMAP distance for now:
  #weighted_edges = hooke:::get_weighted_edges(ccm, all_edges)
  weighted_edges = hooke:::weigh_edges_by_umap_dist(ccm, cov_graph_edges)

  # if (weigh_by_pcor) {
  #   weighted_edges = weighted_edges %>%
  #     mutate(weight = 1/pcor)
  #     # mutate(y = 1/pcor, z=dist,
  #     #        y_z = (y-min(y))/(max(y)-min(y)),
  #     #        z_z = (z-min(z))/(max(z)-min(z)),
  #     #        weight = y_z + z_z) %>% select(-c(y,z,y_z,z_z))
  # }

  paga_components = igraph::components(paga_graph)
  same_partition_mat = outer(paga_components$membership, paga_components$membership, FUN="==")
  weighted_edges = weighted_edges %>% group_by(from, to) %>%
    mutate(adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to),
           same_partition = same_partition_mat[from, to]) %>% ungroup()


  weighted_edges = weighted_edges %>% filter(adjacent_in_paga)

  pathfinding_graph = weighted_edges %>% select(from, to, weight) %>%
    igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
    igraph::as.directed()

  if (allow_links_between_components){
    pathfinding_graph = add_cross_component_pathfinding_links(ccm, pathfinding_graph)
    pathfinding_graph = igraph::as_data_frame(pathfinding_graph) %>%  select(from, to)
    pathfinding_graph = hooke:::weigh_edges_by_umap_dist(ccm, cov_graph_edges) %>%
      igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
      igraph::as.directed()
  }

  # edges_between_concurrent_states = left_join(pathfinding_graph %>% igraph::as_data_frame(what="edges"),
  #                                             extant_cell_type_df, by=c("from"="cell_group")) %>%
  #   select(from, to, from_start=longest_contig_start, from_end=longest_contig_end)
  # edges_between_concurrent_states = left_join(edges_between_concurrent_states, extant_cell_type_df, by=c("to"="cell_group")) %>%
  #   select(from, to, from_start, from_end, to_start=longest_contig_start, to_end=longest_contig_end)
  # edges_between_concurrent_states = edges_between_concurrent_states %>%
  #   filter(from_start <= to_start & from_end >= to_start) %>% select(from, to) %>% distinct()

  #if (require_presence_at_all_timepoints){
  #pathfinding_graph = igraph::intersection(pathfinding_graph,
  #                                         edges_between_concurrent_states %>%
  #                                         igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata))
  #}
  return (pathfinding_graph)
}

get_perturbation_paths <- function(perturbation_ccm,
                                                  time_window,
                                                  pathfinding_graph,
                                                  interval_step,
                                                  interval_col,
                                                  log_abund_detection_thresh,
                                                  q_val,
                                   model_for_pcors="reduced",
                                   min_pcor_for_edges=0.05,
                                                  ...){

  #print ("getting perturbation paths")
  #print (time_window)

  start_time = min(time_window$start_time)
  stop_time = min(time_window$stop_time)

  #print (start_time)
  #print (stop_time)

  earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
                                                   start_time=start_time,
                                                   stop_time=stop_time,
                                                   interval_step = interval_step,
                                                   interval_col=interval_col,
                                                   log_abund_detection_thresh=log_abund_detection_thresh,
                                                   q_val = q_val,
                                                   ...)
  #print (earliest_loss_tbl)


  pln_model = model(perturbation_ccm, model_for_pcors)

  #print (pln_model)

  #cov_graph <- return_igraph(pln_model)
  #cov_edges <- igraph::as_data_frame(cov_graph, what="edges")
  #print (cov_edges)
  #cov_edges = cov_edges %>%
  #  dplyr::filter(weight > min_pcor_for_edges)

  #
  #change_corr_tbl = cov_edges %>% dplyr::select(from, to, weight) %>% dplyr::rename(pcor = weight)
  change_corr_tbl = earliest_loss_tbl %>% filter (is.finite(earliest_time)) %>% select(cell_group)
  change_corr_tbl = change_corr_tbl %>% select(cell_group) %>% tidyr::expand(cell_group, cell_group)
  colnames(change_corr_tbl) = c("from", "to")
  change_corr_tbl = change_corr_tbl %>% filter(from != to)

  lost_cell_groups = change_corr_tbl %>% pull(from) %>% unique

  lost_cell_group_pathfinding_subgraph = igraph::subgraph(pathfinding_graph, lost_cell_groups)

  #print ("change corr tbl")
  #print (change_corr_tbl)

  #corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
  change_corr_tbl = dplyr::left_join(change_corr_tbl, earliest_loss_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
  #dplyr::rename(log_abund_x,
  #              to_delta_log_abund = delta_log_abund)
  change_corr_tbl = dplyr::left_join(change_corr_tbl, earliest_loss_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%

  change_corr_tbl = change_corr_tbl %>% mutate_if(is.numeric, tidyr::replace_na, replace = -Inf)

  #print (change_corr_tbl)

  #print ("possible loss pairs")
  #print (change_corr_tbl)

  #print ("debug pairs")
  #change_corr_tbl %>% filter(from == "24") %>% print
  concordant_fwd_loss_pairs = change_corr_tbl %>%
    #dplyr::filter(pcor < 0) %>% # do we just want negative again?
    dplyr::filter((is.finite(from_peak_loss_time) & is.finite(to_peak_loss_time) & from_peak_loss_time <= to_peak_loss_time) |
                  (is.finite(from_largest_loss_time) & is.finite(to_largest_loss_time) & from_largest_loss_time <= to_largest_loss_time) |
                  (is.finite(from_earliest_time) & is.finite(to_earliest_time) & from_earliest_time < to_earliest_time))

  #print ("loss pairs selected based on ordering")
  #print (concordant_fwd_loss_pairs)

  # Compute the shortest paths between each of those node pairs in the pathfinding graph
  concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>% select(from, to) %>% distinct()
  paths_between_concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>%
    mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                              .x = from, .y = to,
                              lost_cell_group_pathfinding_subgraph))

  #print ("loss pair paths")
  #print (paths_between_concordant_fwd_loss_pairs)

  #print ("debug pair paths")
  #paths_between_concordant_fwd_loss_pairs %>% filter(is.na(path) == FALSE) %>% rename(origin=from, destination=to) %>% tidyr::unnest(path) %>% print(n=1000)
  #print (paths_between_concordant_fwd_loss_pairs)#  %>% tidyr::unnest(path) %>% filter(from == "24") %>% print

  paths_between_concordant_fwd_loss_pairs = score_paths_for_perturbations(perturbation_ccm, paths_between_concordant_fwd_loss_pairs, perturbation_col="knockout", interval_col)
  return (paths_between_concordant_fwd_loss_pairs)
}


assess_graph_for_perturbation_support <- function(perturbation_ccm,
                                                  time_window,
                                                         pathfinding_graph,
                                                         interval_step,
                                                         interval_col,
                                                         log_abund_detection_thresh,
                                                         q_val,
                                                         ...){

  start_time = min(time_window$start_time)
  stop_time = min(time_window$stop_time)

  earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
                                                   start_time=start_time,
                                                   stop_time=stop_time,
                                                   interval_step = interval_step,
                                                   interval_col=interval_col,
                                                   log_abund_detection_thresh=log_abund_detection_thresh,
                                                   q_val = q_val,
                                                   ...)

  pathfinding_graph_as_df = igraph::as_data_frame(pathfinding_graph)

  #dplyr::rename(log_abund_x,
  #              to_delta_log_abund = delta_log_abund)
  pathfinding_graph_as_df = dplyr::left_join(pathfinding_graph_as_df, earliest_loss_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%
  pathfinding_graph_as_df = dplyr::left_join(pathfinding_graph_as_df, earliest_loss_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%

  #  dplyr::rename(from_delta_log_abund = delta_log_abund)
  pathfinding_graph_as_df = pathfinding_graph_as_df %>% mutate_if(is.numeric, coalesce, -Inf)

  print ("Here's the key pathfinding info")
  #print (pathfinding_graph_as_df)

  #pathfinding_graph_as_df %>% select (from, to, from_earliest_time, from_latest_time, to_earliest_time, to_latest_time) %>% print(n=1000)
  #

  # Include edges (from, to) where to is not lost before from:
  #pathfinding_graph_as_df = pathfinding_graph_as_df %>%
    #filter(!(to_earliest_time < from_earliest_time))


  # We DO WANT "from lost before to" - this is consistent with the edge direction
  # We DO NOT WANT "to lost before from" - this violates the edge direction
  pathfinding_graph_as_df = pathfinding_graph_as_df %>%
    mutate(to_lost = is.finite(to_earliest_time),
           from_lost = is.finite(from_earliest_time),
           to_lost_within_window = to_lost & to_earliest_time > start_time,
           from_lost_within_window = from_lost & from_earliest_time > start_time,
           from_lost_before_to = from_lost & to_lost & from_earliest_time < to_earliest_time,
           to_lost_before_from = to_lost & to_earliest_time < from_earliest_time) # this case also covers when to is lost but from is not. Also bad.

  # pathfinding_graph_as_df = pathfinding_graph_as_df %>%
  #   mutate(conflicts_with_perturb = to_lost_before_from,
  #          supported_by_perturb = from_lost_before_to)

  pathfinding_graph_as_df = pathfinding_graph_as_df %>%
    mutate(conflicts_with_perturb = to_lost_before_from,
           supported_by_perturb = from_lost_before_to)

  #pathfinding_graph_as_df %>% print (n=1000)

  # If for edge (from, to) from is lost, include the edge if from is lost at or after that :
  #pathfinding_graph_as_df = pathfinding_graph_as_df %>%
  #  filter(from_lost_before_to & !to_lost_before_from)


  print (pathfinding_graph_as_df)

  pathfinding_graph_revised = pathfinding_graph_as_df %>%
    igraph::graph_from_data_frame(directed=TRUE,  vertices=data.frame(id=igraph::V(pathfinding_graph)$name))
  #plot_state_graph_annotations(tbx16_ccm, xxx_pathfinding_graph_revised, group_nodes_by="cell_type", edge_size=2, hide_unlinked_nodes=F)



  #if (require_presence_at_all_timepoints){
  #pathfinding_graph = igraph::intersection(pathfinding_graph,
  #                                         edges_between_concurrent_states %>%
  #                                         igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata))
  #}
  return (pathfinding_graph_revised)
}


# 2) Make sure this graph is acyclic by deleting problematic edges
# from https://github.com/sachsmc/causaloptim/blob/master/R/graph-utilities.R

#' Find cycles in a graph
#'
#' @param g an igraph object
#' @return A list of vectors of integers, indicating the vertex sequences for the cycles found in the graph
#' @export
find_cycles = function(g) {
  Cycles = NULL
  for(v1 in igraph::V(g)) {
    if(igraph::degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = igraph::neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(igraph::all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 2)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}


#' score a path based on fitting a linear model of time ~ geodeisic distance
#' @param ccs
#' @param path_df
measure_time_delta_along_path <- function(path_df, ccs, interval_col="timepoint") {


  #cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df = normalized_counts(ccs, "size_only", pseudocount = 0)[vertices,] %>%
    as.matrix() %>%
    Matrix::t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(!matches("rowname")) %>%
    rename(sample=rowname, cell_group=name, num_cells=value)

  cells_along_path_df = cells_along_path_df %>%
    left_join(colData(ccs) %>% as.data.frame %>%
                select(sample, !!sym(interval_col)) %>%
                as_tibble(), by = c("sample" = "sample"))

  cells_along_path_df = cells_along_path_df %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  cells_along_path_df = cells_along_path_df %>% filter(num_cells > 0)
  cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$y = cells_along_path_df[[interval_col]]

  path_model = lm(y ~ geodesic_dist, data=cells_along_path_df, weights=cells_along_path_df$num_cells)

  path_model_tidied = broom::tidy(path_model)
  path_model_glanced = broom::glance(path_model)
  pm_stats = tibble(time_dist_effect = unlist(path_model_tidied[2, "estimate"]),
                    time_dist_effect_pval = unlist(path_model_tidied[2, "p.value"]),
                    time_dist_model_adj_rsq = unlist(path_model_glanced[1, "adj.r.squared"]),
                    time_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]),
                    path_length = max(path_df$geodesic_dist))
  return(pm_stats)
  #return(coef(path_model)[["distance_from_root"]])
}
#debug(cells_along_path)

#' score a path based on fitting a linear model of perturbation ~ geodesic distance
#' @param ccs
#' @param path_df
measure_perturbation_freq_along_path <- function(path_df, ccs, perturbation_col="knockout", interval_col="timepoint") {


  #cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df = normalized_counts(ccs, "size_only", pseudocount = 0)[vertices,] %>%
    as.matrix() %>%
    Matrix::t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(!matches("rowname")) %>%
    rename(sample=rowname, cell_group=name, num_cells=value)

  cells_along_path_df = cells_along_path_df %>%
    left_join(colData(ccs) %>% as.data.frame %>%
                select(sample, !!sym(perturbation_col), !!(interval_col)) %>%
                as_tibble(), by = c("sample" = "sample"))

  cells_along_path_df = cells_along_path_df %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  #print (cells_along_path_df)
  cells_along_path_df = cells_along_path_df %>% mutate_if(is.numeric, tidyr::replace_na, replace = 0)
  cells_along_path_df = cells_along_path_df %>% dplyr::filter (num_cells > 0)
  #cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  #cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$perturb = cells_along_path_df[[perturbation_col]]
  cells_along_path_df$timepoint = cells_along_path_df[[interval_col]]

  path_model = glm(#perturb ~ timepoint + geodesic_dist,
    perturb ~ geodesic_dist,
                   data=cells_along_path_df,
                   weights=cells_along_path_df$num_cells,
                   family=binomial(),
                   singular.ok = TRUE)

  path_model_tidied = broom::tidy(path_model)
  path_model_glanced = broom::glance(path_model)
  path_model_glanced$pseudo.r.squared = 1 - (path_model_glanced$deviance / path_model_glanced$null.deviance)

  pm_stats = tibble(perturb_dist_effect = unlist(path_model_tidied[2, "estimate"]),
                    perturb_dist_effect_pval = unlist(path_model_tidied[2, "p.value"]),
                    perturb_dist_model_adj_rsq = unlist(path_model_glanced[1, "pseudo.r.squared"]),
                    perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]))
  # if (unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) { # should not happen
  #   print (cells_along_path_df %>% as.data.frame)
  #   print (summary(path_model))
  # }

  if (path_model$converged == FALSE | path_model$boundary == TRUE | unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0)
  {
    return(tibble(perturb_dist_effect = 0,
                  perturb_dist_effect_pval = 1,
                  perturb_dist_model_adj_rsq = 0,
                  perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"])))
  }
  return(pm_stats)
  #return(coef(path_model)[["distance_from_root"]])
}
#debug(cells_along_path)




#' Identify the possible origins for each destination
#'
build_timeseries_transition_dag <- function(ccm,
                                            extant_cell_type_df,
                                            pathfinding_graph,
                                            q_val=0.01,
                                            start_time = NULL,
                                            stop_time = NULL,
                                            interval_col="timepoint",
                                            interval_step = 2,
                                            min_interval = 4,
                                            max_interval = 24,
                                            min_dist_vs_time_r_sq = 0,
                                            ...) {

  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start_time)){
    start_time = min(colData(ccm@ccs)[,interval_col])
  }
  if (is.null(stop_time)){
    stop_time = max(colData(ccm@ccs)[,interval_col])
  }

  timepoints = seq(start_time, stop_time, interval_step)

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start_time, stop_time, interval_col=interval_col, interval_step=interval_step, ...)

  select_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
    cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
    cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval)

  relevant_comparisons = time_contrasts %>%
    mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    interval_col=interval_col,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                  .x = comp_abund,
                                  ccm = ccm))
  relevant_comparisons = relevant_comparisons %>%
    tidyr::unnest(rec_edges) %>%
    #dplyr::filter(pcor < 0) %>% # do we just want negative again?
    dplyr::filter((from_delta_log_abund > 0 & to_delta_log_abund < 0) |
                    (to_delta_log_abund > 0 & from_delta_log_abund < 0)) %>%
    dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val)

  edge_union = relevant_comparisons %>% select(from, to) %>% distinct()


  #print ("finding shortest paths...")
  #print(head(relevant_comparisons))
  paths_for_relevant_edges = edge_union %>%
    mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                              .x = from, .y = to,
                              pathfinding_graph))


  #print ("scoring paths for time flow...")
  #print(head(paths_for_relevant_edges))
  paths_for_relevant_edges = paths_for_relevant_edges %>%
    filter(!is.na(path)) %>%
    mutate(time_vs_distance_model_stats = purrr::map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                     .x = path,
                                                     ccs=ccm@ccs,
                                                     interval_col=interval_col)) %>%
    tidyr::unnest(time_vs_distance_model_stats)
  #cells_along_path(ccs, paths_for_relevant_edges$path[[1]] %>% dplyr::select(from, to))


  #print ("tabulating model stats...")
  #print(head(paths_for_relevant_edges))
  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(#time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq,
    time_dist_model_score = time_dist_model_ncells * time_dist_model_adj_rsq,
    path_score = time_dist_model_score)
  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(time_dist_effect_q_val = p.adjust(time_dist_effect_pval, method="BH"))

  #print (paths_for_relevant_edges)
  selected_paths = paths_for_relevant_edges %>% filter (time_dist_effect > 0 &
                                                          time_dist_effect_pval < q_val &
                                                          time_dist_model_adj_rsq > min_dist_vs_time_r_sq) %>%
    group_by(to) %>%
    #slice_max(dist_model_score, n=3) %>%
    ungroup() %>% arrange(desc(time_dist_model_score)) %>%
    dplyr::mutate(path_contrast="time")

  #selected_paths %>% select(origin=from, destination=to, path, path_score) %>% tidyr::unnest(path) %>% arrange(origin, destination) %>% filter(from %in% c("20", "38") & to %in% c("20", "38"))  %>% print(n=1000)

  #G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = FALSE)

  G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = FALSE)

  covered_G = compute_min_path_cover(ccm, G, weight_attribute="support")

  edge_support = igraph::as_data_frame(covered_G) %>% select(from, to)

  edge_support = left_join(edge_support,
                           selected_paths %>% select(-from, -to) %>% tidyr::unnest(path))

  #print (edge_support)
  edge_support = edge_support %>% group_by(from, to) %>% slice_max(time_dist_model_adj_rsq, n=1)

  #print (edge_support)
  covered_G = igraph::graph_from_data_frame(edge_support, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))

  # missing_pathfinding_edges = dplyr::setdiff(as_tibble(igraph::as_edgelist(pathfinding_graph)), as_tibble(igraph::as_edgelist(G)))
  # colnames(missing_pathfinding_edges) = c("from", "to")
  # G = G + igraph::graph_from_edgelist(as.matrix(missing_pathfinding_edges))
  # #
  #
  # #G = break_cycles_in_state_transition_graph(G)
  #
  # cycle_breaking_scores = igraph::edge_attr(G, "support")
  # cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0
  # igraph::edge_attr(G, "support") = cycle_breaking_scores
  #
  # print (igraph::edge_attr(G, "support"))
  #
  # print ("breaking cycles...")
  # G = break_cycles_in_state_transition_graph(G)
  # print ("done")

  return(covered_G)

}

get_paths_between_recip_time_nodes <- function(ccm,
                                               pathfinding_graph,
                                               timepoint_pred_df,
                                               timepoints,
                                               min_interval,
                                               max_interval,
                                               q_val,
                                               interval_col){
  #timepoints = sort(unique(timepoint_pred_df[[interval_col]]))

  compare_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
    cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
    cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval)

  # Find the nodes that undergo reciprocal fold changes between successive time points
  recip_time_node_pairs = time_contrasts %>%
    mutate(comp_abund = purrr::map2(.f = compare_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    interval_col=interval_col,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                  .x = comp_abund,
                                  ccm = ccm))
  recip_time_node_pairs = recip_time_node_pairs %>%
    tidyr::unnest(rec_edges)


  recip_time_node_pairs = recip_time_node_pairs %>%
    #dplyr::filter(pcor < 0) %>% # do we just want negative again?
    dplyr::filter((from_delta_log_abund > 0 & to_delta_log_abund < 0) |
                    (to_delta_log_abund > 0 & from_delta_log_abund < 0)) %>%
    dplyr::filter(from_delta_p_value < q_val & to_delta_p_value < q_val)

  # recip_time_node_pairs %>% filter(from == "15") %>% select(from,
  #                                                           to,
  #                                                           from_adjusted_timepoint_x,
  #                                                            from_adjusted_timepoint_y,
  #                                                            from_delta_log_abund,
  #                                                            from_delta_p_value,
  #                                                            to_delta_log_abund,
  #                                                            to_delta_p_value) %>%
  #                                                            as.data.frame %>% print

  # Compute the shortest paths between each of those node pairs in the pathfinding graph
  recip_time_node_pairs = recip_time_node_pairs %>% select(from, to) %>% distinct()
  paths_between_recip_time_nodes = recip_time_node_pairs %>%
    mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                              .x = from, .y = to,
                              pathfinding_graph))

  #paths_between_recip_time_nodes %>% filter (from == "15") %>% print

  # Score each shortest path for correlation between time and distance along it
  paths_between_recip_time_nodes = paths_between_recip_time_nodes %>%
    filter(!is.na(path)) %>%
    mutate(time_vs_distance_model_stats = purrr::map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                     .x = path,
                                                     ccs=ccm@ccs,
                                                     interval_col=interval_col)) %>%
    tidyr::unnest(time_vs_distance_model_stats)



  paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% mutate(time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq)
  paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% mutate(time_dist_effect_q_val = p.adjust(time_dist_effect_pval, method="BH"))
  #paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% filter (dist_effect > 0 & dist_effect_q_val < q_val & dist_model_adj_rsq > min_dist_vs_time_r_sq) %>%
  #  group_by(to) %>%
  #  #slice_max(dist_model_score, n=3) %>%
  #  ungroup() %>% arrange(desc(dist_model_score))
  return(paths_between_recip_time_nodes)
}



get_timeseries_paths <- function(ccm,
                                 extant_cell_type_df,
                                 pathfinding_graph,
                                 start_time,
                                 stop_time,
                                 #perturbation_col,
                                 interval_col,
                                 interval_step,
                                 min_interval,
                                 max_interval,
                                 q_val = 0.01,
                                 ...)
{
  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start_time)){
    start_time = min(colData(ccm@ccs)[,interval_col])
  }
  if (is.null(stop_time)){
    stop_time = max(colData(ccm@ccs)[,interval_col])
  }

  # FIXME: shouldn't use hardcoded knockout here
  wt_timepoint_pred_df = estimate_abundances_over_interval(ccm,
                                                           start_time,
                                                           stop_time,
                                                           interval_col=interval_col,
                                                           interval_step=interval_step,
                                                           ...)
  #ko_timepoint_pred_df = estimate_abundances_over_interval(ccm, start_time, stop_time, knockout=TRUE, interval_col=interval_col, interval_step=interval_step, ...)

  timepoints = seq(start_time, stop_time, interval_step)


  paths_between_recip_time_nodes = get_paths_between_recip_time_nodes(ccm,
                                                                      pathfinding_graph,
                                                                      wt_timepoint_pred_df,
                                                                      timepoints,
                                                                      min_interval,
                                                                      max_interval,
                                                                      q_val,
                                                                      interval_col)
  return (paths_between_recip_time_nodes)
}

select_paths_from_pathfinding_graph <- function(pathfinding_graph, selected_paths, allow_cycles=FALSE)
{
  G = pathfinding_graph
  G = igraph::delete_edges(G, igraph::E(G))

  for (i in (1:nrow(selected_paths))){
    next_path = selected_paths$path[[i]] %>% select(from, to)
    next_path_graph = next_path %>% igraph::graph_from_data_frame(directed=TRUE)
    G_prime = igraph::union(G, next_path_graph)
    if (allow_cycles | length(find_cycles(G_prime)) == 0){
      G = G_prime
    }else{
      # Debug:
      #print ("skipping due to cycle:")
      #print (selected_paths[i,] %>% select(-path))
    }
  }
  igraph::E(G)$transcriptome_dist = igraph::E(pathfinding_graph)[igraph::E(G)]$weight

  #print (selected_paths)
  support_tbl= selected_paths %>% dplyr::select(path_contrast, path_score, path) %>% tidyr::unnest(path)
  edge_support_summary = support_tbl %>% group_by(from, to) %>% dplyr::select(from, to, path_score, path_contrast)  %>%
    distinct() %>% summarize(support=sum(path_score))

  edge_support_summary = edge_support_summary %>% mutate(support = ifelse(is.na(support), 0, support))
  #print (edge_support_summary %>% as.data.frame)

  #annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% left_join(edge_support)
  annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% distinct() %>% left_join(edge_support_summary)
  annotated_G = igraph::graph_from_data_frame(annotated_G, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))

  return (annotated_G)
}

break_cycles_in_state_transition_graph <- function(state_graph)
{
  print ("FAS weights")
  print (igraph::E(state_graph)$support)
  fas = igraph::feedback_arc_set(state_graph, weights=igraph::E(state_graph)$support)
  print ("feedback arc set")
  print (fas)
  state_graph = igraph::delete_edges(state_graph, fas)
  print ("checking for cycles")
  curr_cycles = find_cycles(state_graph)
  if (length(curr_cycles) > 0){
    print ("warning: cycles remain")
    print (curr_cycles)
  }
  return(state_graph)
  # updated_state_graph = state_graph
  # curr_cycles = find_cycles(state_graph)
  # while(length(curr_cycles) > 0){
  #
  #   for (cycle in curr_cycles){
  #     cycle_subgraph = igraph::subgraph(updated_state_graph, cycle)
  #     min_score_edges = which(igraph::E(cycle_subgraph)$transcriptome_dist == max(igraph::E(cycle_subgraph)$transcriptome_dist))
  #     #min_weight_edge_names = E(cycle_subgraph)[min_weight_edges]$name
  #     #print (cycle_subgraph)
  #     #print (igraph::E(cycle_subgraph)[min_score_edges]$name)
  #     #cycle_subgraph = igraph::delete_edges(cycle_subgraph, min_score_edges)
  #
  #     #cycle_subgraph = igraph::delete_edge_attr(cycle_subgraph, "transcriptome_dist")
  #     print (paste("state graph has ", length(igraph::E(updated_state_graph)), "edges"))
  #     #updated_state_graph = igraph::intersection(updated_state_graph, cycle_subgraph, byname=TRUE)
  #     updated_state_graph <<- updated_state_graph - igraph::subgraph(cycle_subgraph, min_score_edges)
  #     print (paste("\t", length(igraph::E(updated_state_graph)), "remaining"))
  #   }
  #
  #   print (paste("removed", length(curr_cycles), "cycles"))
  #   #break
  #   curr_cycles <<- find_cycles(updated_state_graph)
  #   print (paste("\t", length(curr_cycles), "remaining"))
  # }
  # return(updated_state_graph)
}

get_pcor_between_pair <- function(perturb_model_pcor_matrix, from, to){
  #edges_to_test = mapply(function(x, y) { return(c(x,y))}, from, to)
  pcor_vals = mapply(function(x, y) { return(perturb_model_pcor_matrix[x,y])}, from, to)
  #pcor_val = igraph::E(pcor_graph)[igraph::get.edge.ids(pcor_graph, edges_to_test)]$weight
  #pcor_val = as.vector(perturb_model_pcor_matrix[edges_to_test])
  return(pcor_vals)
}

score_paths_for_perturbations <- function(perturbation_ccm, paths_between_recip_time_nodes, perturbation_col="knockout", interval_col="timepoint")
{
  #net <- model$latent_network(type = type)
  #perturb_model_pcor_matrix = model(perturbation_ccm, "reduced")$latent_network(type = "partial_cor")
  # Now measure each time-flow path for a drop in perturbed cells along it
  #print ("scoring paths with model:")
  #print (perturb_model_pcor_graph)
  #print (paths_between_recip_time_nodes)

  paths_between_perturb_vs_wt_node_pairs = paths_between_recip_time_nodes %>%
    filter(!is.na(path)) %>%
    #mutate(pcor_between_node_pair =  get_pcor_between_pair(perturb_model_pcor_matrix, from, to)) %>%
    mutate(perturb_vs_distance_model_stats = purrr::map(.f = purrr::possibly(measure_perturbation_freq_along_path, NA_character_),
                                                        .x = path,
                                                        ccs=perturbation_ccm@ccs,
                                                        perturbation_col=perturbation_col,
                                                        interval_col=interval_col)) %>%
    filter(!is.na(perturb_vs_distance_model_stats)) %>%
    tidyr::unnest(perturb_vs_distance_model_stats)

  #print (head(paths_between_perturb_vs_wt_node_pairs))
  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(path_score = ifelse(perturb_dist_effect < 0, -perturb_dist_effect, 0))
  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(path_score = ifelse(perturb_dist_effect < 0, -perturb_dist_effect * perturb_dist_model_adj_rsq, 0))
  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(pertubation_path_score = ifelse(perturb_dist_effect < 0, perturb_dist_model_adj_rsq, 0))
  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(path_score = ifelse(perturb_dist_effect < 0, 1, 0))

  #print (paths_between_perturb_vs_wt_node_pairs)

  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_model_score = -perturb_dist_effect * perturb_dist_model_adj_rsq)
  #paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_effect_q_val = p.adjust(perturb_dist_effect_pval, method="BH"))
  return(paths_between_perturb_vs_wt_node_pairs)
}

compare_ko_to_wt_at_timepoint <- function(tp, perturbation_ccm, wt_pred_df, ko_pred_df, interval_col)  {
  cond_wt = wt_pred_df %>% filter(!!sym(interval_col) == tp)
  cond_ko = ko_pred_df %>% filter(!!sym(interval_col) == tp)
  return(compare_abundances(perturbation_ccm, cond_wt, cond_ko))
}

estimate_loss_timing <- function(perturbation_ccm,
                                 start_time,
                                 stop_time,
                                 interval_step,
                                 log_abund_detection_thresh=-2,
                                 delta_log_abund_loss_thresh=0,
                                 interval_col="timepoint",
                                 q_val=0.01,
                                 ...){


  wt_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
  ko_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=TRUE, interval_col=interval_col, interval_step=interval_step, ...)

  #print (wt_timepoint_pred_df)
  #print (ko_timepoint_pred_df)
  timepoints = seq(start_time, stop_time, interval_step)

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
    mutate(comp_abund = purrr::map(.f = compare_ko_to_wt_at_timepoint,
                                   .x = t1,
                                   perturbation_ccm=perturbation_ccm,
                                   interval_col=interval_col,
                                   wt_pred_df = wt_timepoint_pred_df,
                                   ko_pred_df = ko_timepoint_pred_df)) %>% tidyr::unnest(comp_abund)
    #mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
    #                              .x = comp_abund,
    #                              ccm = perturbation_ccm)) %>% tidyr::unnest(rec_edges)

  #perturb_vs_wt_nodes %>% print (n=10000)
  #print(sym(paste(interval_col, "_x", sep="")))
  #print (perturb_vs_wt_nodes)
  #print ("subsetting loss tbl")
  #print (q_val)
  #perturb_vs_wt_nodes %>% filter(cell_group == "1")  %>%
    #filter(is.na(adjusted_timepoint_x) == FALSE) %>% as.data.frame %>%
  #  select(cell_group, adjusted_timepoint_x, log_abund_x, log_abund_se_x, log_abund_y, log_abund_se_y, delta_log_abund, delta_p_value) %>%
  #  print(n=1000)


  earliest_loss_tbl =  perturb_vs_wt_nodes %>%
    filter(delta_log_abund < 0 & delta_p_value < q_val)

  #print (earliest_loss_tbl)


  peak_wt_abundance = perturb_vs_wt_nodes %>% group_by(cell_group) %>% slice_max(log_abund_x, n=1)
  loss_at_peak_wt_abundance = peak_wt_abundance %>% mutate(cell_group,
                                                           peak_wt_time = !!sym(paste(interval_col, "_x", sep="")),
                                                           lost_at_peak = delta_log_abund < 0 & delta_p_value < q_val,
                                                           peak_loss_time = ifelse(lost_at_peak, peak_wt_time, NA))
  loss_at_peak_wt_abundance = loss_at_peak_wt_abundance %>% select(cell_group, peak_wt_time, lost_at_peak, peak_loss_time)

  #loss_at_peak_wt_abundance %>% filter(cell_group == "1") %>% print()

  largest_losses = perturb_vs_wt_nodes %>%
    filter(delta_log_abund < 0 & delta_p_value < q_val) %>%
    group_by(cell_group) %>%
    slice_min(delta_log_abund, n=1) %>% select(cell_group, largest_loss_time=!!sym(paste(interval_col, "_x", sep="")), largest_loss=delta_log_abund)

  print (largest_losses)
  earliest_loss_tbl = earliest_loss_tbl %>%
    group_by(cell_group) %>%
    #mutate(largest_loss = min(delta_log_abund),
    #       largest_loss_time = !!sym(paste(interval_col, "_x", sep=""))[which(delta_log_abund == largest_loss)]) %>%
    summarize(#largest_loss_time = min(largest_loss_time),
              earliest_time = min(!!sym(paste(interval_col, "_x", sep=""))),
              latest_time = max(!!sym(paste(interval_col, "_x", sep=""))))

  loss_tbl = earliest_loss_tbl %>% left_join(largest_losses)
  loss_tbl = loss_tbl %>% left_join(loss_at_peak_wt_abundance)


  #print (earliest_loss_tbl)

  return(loss_tbl)
}

score_graph_edges_for_perturbations <- function(perturbation_ccm, state_transition_graph, start_time, stop_time, interval_step, interval_col="timepoint", q_val=0.01, ...)
{

  wt_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
  ko_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=TRUE, interval_col=interval_col, interval_step=interval_step, ...)

  compare_ko_to_wt_at_timepoint <- function(tp, wt_pred_df, ko_pred_df, interval_col)  {
    cond_wt = wt_pred_df %>% filter(!!sym(interval_col) == tp)
    cond_ko = ko_pred_df %>% filter(!!sym(interval_col) == tp)
    return(compare_abundances(perturbation_ccm, cond_wt, cond_ko))
  }

  timepoints = seq(start_time, stop_time, interval_step)

  # Find the pairs of nodes that are both lost in the perturbation at the same time
  perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
    mutate(comp_abund = purrr::map(.f = compare_ko_to_wt_at_timepoint,
                                    .x = t1,
                                    interval_col=interval_col,
                                   wt_pred_df = wt_timepoint_pred_df,
                                   ko_pred_df = ko_timepoint_pred_df)) %>%
    mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                  .x = comp_abund,
                                  ccm = perturbation_ccm)) %>% tidyr::unnest(rec_edges)

  perturb_vs_wt_edges = state_transition_graph %>% igraph::as_data_frame() %>% dplyr::select(from, to)

  #print (perturb_vs_wt_nodes)
  #print(sym(paste(interval_col, "_x", sep="")))
  earliest_loss_tbl =  perturb_vs_wt_nodes %>%
    filter(delta_log_abund < 0 & delta_p_value < q_val) %>%
    group_by(cell_group) %>% summarize(earliest_time = min(!!sym(paste(interval_col, "_x", sep=""))))

  #print (earliest_loss_tbl)

  #print (perturb_vs_wt_edges)
  #print (perturb_vs_wt_nodes)
  #corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
  perturb_vs_wt_edges = dplyr::left_join(perturb_vs_wt_edges, perturb_vs_wt_nodes %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
  #dplyr::rename(log_abund_x,
  #              to_delta_log_abund = delta_log_abund)
  perturb_vs_wt_edges = dplyr::left_join(perturb_vs_wt_edges, perturb_vs_wt_nodes %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%
  #  dplyr::rename(from_delta_log_abund = delta_log_abund)

  perturb_vs_wt_edges = dplyr::left_join(perturb_vs_wt_edges, earliest_loss_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
  perturb_vs_wt_edges = dplyr::left_join(perturb_vs_wt_edges, earliest_loss_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%

  #print (perturb_vs_wt_edges)
  # First just find the pairs of cell types that are lost concurrently,
  perturb_vs_wt_edges = perturb_vs_wt_edges %>%
    #dplyr::filter(pcor < 0) %>% # do we just want negative again?
    #dplyr::filter((from_delta_log_abund > 0 & to_delta_log_abund < 0) |
    #                (from_delta_log_abund < 0 & to_delta_log_abund > 0)) %>%
    dplyr::filter((from_delta_log_abund < 0 & to_delta_log_abund < 0)) %>%
    dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val) %>%
    dplyr::filter(from_earliest_time < to_earliest_time)

  perturb_vs_wt_edges = perturb_vs_wt_edges %>% select(from, to) %>% group_by(from, to) %>% summarize(num_intervals_supported=n())
  return(perturb_vs_wt_edges)
}

#' Identify the possible origins for each destination
#'
build_perturbation_transition_dag <- function(ccm,
                                              extant_cell_type_df,
                                              pathfinding_graph,
                                              q_val=0.01,
                                              start_time = NULL,
                                              stop_time = NULL,
                                              perturbation_col="knockout",
                                              interval_col="timepoint",
                                              interval_step = 2,
                                              min_interval = 4,
                                              max_interval = 24,
                                              min_dist_vs_time_r_sq = 0,
                                              ...) {

  paths_between_recip_time_nodes = get_timeseries_paths(ccm,
                                                        extant_cell_type_df,
                                                        pathfinding_graph,
                                                        start_time = start_time,
                                                        stop_time = stop_time,
                                                        interval_col=interval_col,
                                                        interval_step = interval_step,
                                                        min_interval = min_interval,
                                                        max_interval = max_interval,
                                                        knockout=FALSE,
                                                        ...)



  paths_between_perturb_vs_wt_node_pairs = score_paths_for_perturbations(ccm, paths_between_recip_time_nodes)
  paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(path_score = perturb_dist_model_ncells * perturb_dist_model_adj_rsq)
  paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_effect_q_val = p.adjust(perturb_dist_effect_pval, method="BH"))

  #origin_dest_node_pairs = dplyr::intersect(perturb_vs_wt_node_pairs, recip_time_node_pairs)

  origin_dest_node_pairs = dplyr::inner_join(paths_between_recip_time_nodes, paths_between_perturb_vs_wt_node_pairs)
  #print (origin_dest_node_pairs)
  selected_paths = origin_dest_node_pairs %>% filter (perturb_dist_effect < 0 & perturb_dist_effect_q_val < q_val &
                                                        time_dist_effect > 0 & time_dist_effect_q_val < q_val) %>%
    group_by(to) %>%
    #slice_max(dist_model_score, n=3) %>%
    ungroup() %>% arrange(desc(path_score)) %>%
    mutate(path_contrast="time_perturb")


  G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles=FALSE)
  return (G)
}


compute_min_path_cover <- function(ccm, G, weight_attribute="weight"){

  igraph::edge_attr(G, "weight") = igraph::edge_attr(G, weight_attribute)
  transitive.closure <- function(g,mat=FALSE,loops=TRUE){
    g <- igraph::as_adjacency_matrix(g, attr="weight")

    n <- ncol(g)

    matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
      while(i < pow) {
        z <- z %*% x
        y <- y+z
        i <- i+1
      }
      return(y)
    }

    h <- matExpIterativ(g,n)
    h <- (h>0)*1
    dimnames(h) <- dimnames(g)
    if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
    if (!mat) h = igraph::graph_from_adjacency_matrix(h, weighted=TRUE) #h <- as(h,"graphNEL")
    return(h)
  }
  G_tr = transitive.closure(G, loops=F)

  G_split = G_tr %>% igraph::as_data_frame(what="edges")
  G_split$weight = NULL

  cov_graph <- hooke:::return_igraph(model(ccm, "reduced"))
  pcor_mat = cov_graph %>% igraph::as_adjacency_matrix(attr="weight")
  pcor_mat_summ = Matrix::summary(pcor_mat)
  pcor_mat = data.frame(from      = rownames(pcor_mat)[pcor_mat_summ$i],
                        to = colnames(pcor_mat)[pcor_mat_summ$j],
                        weight      = pcor_mat_summ$x)
  G_split = left_join(G_split, pcor_mat)%>% tidyr::replace_na(list(weight = 0))
  G_split$weight = abs(G_split$weight)

  G_nodes = union(G_split$from,G_split$to)
  split_node_metadata = data.frame(id = c(stringr::str_c("left_", union(G_split$from,G_split$to)),
                                          stringr::str_c("right_", union(G_split$from,G_split$to))))
  G_split$from = stringr::str_c("left_", G_split$from)
  G_split$to = stringr::str_c("right_", G_split$to)
  G_split = igraph::graph_from_data_frame(G_split %>% dplyr::select(from, to, weight), directed=FALSE, vertices=split_node_metadata)

  igraph::V(G_split)$type <- grepl("left", igraph::V(G_split)$name)

  mbm = maxmatching::maxmatching(G_split, weighted = TRUE)
  mbm$matching
  matching <- data.frame(dest_node = mbm$matching, orig_node = names(mbm$matching))
  matching <- subset(matching, grepl("left", orig_node) & grepl("right", dest_node))
  matching = matching %>% dplyr::select(orig_node, dest_node) %>%
    mutate(orig_node = stringr::str_replace_all(orig_node, "left_", ""),
           dest_node = stringr::str_replace_all(dest_node, "right_", "")) %>% dplyr::rename(from=orig_node, to=dest_node)

  #node_dag = G_tr
  # FIXME: use real weights here:
  #igraph::E(node_dag)$weight = 1

  node_dag = G

  possible_origins = names(which(igraph::degree(node_dag, mode="in") == 0))
  possible_termini = names(which(igraph::degree(node_dag, mode="out") == 0))
  node_dag = igraph::add_vertices(node_dag, 2, attr=list("name"=c("source", "sink")))
  source_edge_df = data.frame(from="source", to=possible_origins)
  node_dag = igraph::union(node_dag, source_edge_df %>% igraph::graph_from_data_frame())
  sink_edge_df = data.frame(from=possible_termini, to="sink")
  node_dag = igraph::union(node_dag, sink_edge_df %>% igraph::graph_from_data_frame())

  # replace NA values with weight 0
  igraph::edge_attr(node_dag, "weight") = tidyr::replace_na(igraph::edge_attr(node_dag, "weight"), 0)
  #print(node_dag)
  paths_from_chains = matching %>% as_tibble() %>%
    mutate(chain_leg = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag))

  covered_graph = paths_from_chains %>% select(chain_leg) %>%
    tidyr::unnest(c(chain_leg)) %>%
    igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(G_tr)$name))

  chain_heads = names(which(igraph::degree(covered_graph, mode="in") == 0))
  chain_tails = names(which(igraph::degree(covered_graph, mode="out") == 0))

  source_edge_df = data.frame(from="source", to=chain_heads)
  sink_edge_df = data.frame(from=chain_tails, to="sink")


  paths_to_chain_heads = source_edge_df %>% as_tibble() %>%
    mutate(chain_leg = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag))

  paths_to_chain_tails = sink_edge_df %>% as_tibble() %>%
    mutate(chain_leg = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag))

  covered_graph = paths_to_chain_heads %>% select(chain_leg) %>%
    tidyr::unnest(c(chain_leg)) %>%
    igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(node_dag)$name)) %>%
    igraph::union(covered_graph)

  # covered_graph = paths_to_chain_tails %>% select(chain_leg) %>%
  #   tidyr::unnest() %>%
  #   igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(node_dag)$name)) %>%
  #   igraph::union(covered_graph)

  covered_graph = igraph::simplify(covered_graph)

  covered_graph = igraph::delete_vertices(covered_graph, c("source", "sink"))

  return (covered_graph)
}


#' returns nodes not in paga graph
#' @param  ccm
#'
not_in_paga_graph <- function(ccm) {
  cov_graph = hooke:::return_igraph(model(ccm, "reduced"))
  paga_graph = hooke:::get_paga_graph(ccm@ccs@cds)

  paga_graph = igraph::delete.vertices(igraph::simplify(paga_graph), igraph::degree(paga_graph)==0)
  igraph::delete.vertices(igraph::simplify(cov_graph), igraph::degree(cov_graph)==0)

  not_in_paga = setdiff(igraph::V(cov_graph)$name, igraph::V(paga_graph)$name)
  return(not_in_paga)
}


#' assemble a state transition graph from a timeseries
#' @export
assemble_timeseries_transitions <- function(ccm,
                                            q_val=0.01,
                                            start_time = NULL,
                                            stop_time = NULL,
                                            interval_col="timepoint",
                                            interval_step = 2,
                                            min_interval = 4,
                                            max_interval = 24,
                                            log_abund_detection_thresh=-5,
                                            percent_max_threshold=0,
                                            min_dist_vs_time_r_sq=0,
                                            weigh_by_pcor = F,
                                            ...){

  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start_time,
                                              stop_time,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              ...)

  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.

  pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)


  G = build_timeseries_transition_dag(ccm,
                                      extant_cell_type_df,
                                      pathfinding_graph,
                                      q_val,
                                      start_time,
                                      stop_time,
                                      interval_col,
                                      interval_step,
                                      min_interval,
                                      max_interval,
                                      min_dist_vs_time_r_sq,
                                      ...)

  igraph::V(G)$cell_group = ccm@ccs@info$cell_group

  return(G)
}

#' assemble a state transition graph from a timeseries
#' @export
assemble_perturbation_transitions <- function(ccm,
                                              q_val=0.01,
                                              start_time = NULL,
                                              stop_time = NULL,
                                              perturbation_col="knockout",
                                              interval_col="timepoint",
                                              interval_step = 2,
                                              min_interval = 4,
                                              max_interval = 24,
                                              log_abund_detection_thresh=-5,
                                              percent_max_threshold=0,
                                              min_dist_vs_time_r_sq=0,
                                              weigh_by_pcor = F,
                                              ...){

  # Get a table of the cell types that are in the control
  # FIXME: "knockout" is hard coded and should be a user-defined term in the model
  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start_time,
                                              stop_time,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              perturbation_col=perturbation_col,
                                              knockout=FALSE,
                                              ...)

  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.

  pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)
  pathfinding_graph = assess_graph_for_perturbation_support(ccm,
                                                                   pathfinding_graph,
                                                                   start_time,
                                                                   stop_time,
                                                                   interval_step=interval_step,
                                                                   interval_col=interval_col,
                                                                   percent_max_threshold=percent_max_threshold,
                                                                   log_abund_detection_thresh=log_abund_detection_thresh,
                                                                   q_val=q_val,
                                                                   #perturbation_col=perturbation_col,
                                                                   ...)


  return(pathfinding_graph)

  # G = build_timeseries_transition_dag(ccm,
  #                                     extant_cell_type_df,
  #                                     pathfinding_graph,
  #                                     q_val,
  #                                     start_time,
  #                                     stop_time,
  #                                     interval_col,
  #                                     interval_step,
  #                                     min_interval,
  #                                     max_interval,
  #                                     min_dist_vs_time_r_sq,
  #                                     knockout=FALSE,
  #                                     ...)
  #
  # igraph::V(G)$cell_group = ccm@ccs@info$cell_group

  return(G)

  # G = build_perturbation_transition_dag(ccm,
  #                                       extant_cell_type_df,
  #                                       pathfinding_graph,
  #                                       q_val,
  #                                       start_time,
  #                                       stop_time,
  #                                       perturbation_col,
  #                                       interval_col,
  #                                       interval_step,
  #                                       min_interval,
  #                                       max_interval,
  #                                       min_dist_vs_time_r_sq,
  #                                       ...)
  #
  # covered_G = compute_min_path_cover(ccm, G, "support")
  #
  # igraph::V(covered_G)$cell_group = ccm@ccs@info$cell_group

  return(covered_G)
}

#' #' assemble a state transition graph from a set of perturbations
#' #' @export
#' assemble_transition_graph_from_perturbations <- function(control_timeseries_ccm,
#'                                                          perturbation_ccm_tbl,
#'                                                          q_val=0.01,
#'                                                          start_time = NULL,
#'                                                          stop_time = NULL,
#'                                                          perturbation_col="knockout",
#'                                                          interval_col="timepoint",
#'                                                          interval_step = 2,
#'                                                          min_interval = 4,
#'                                                          max_interval = 24,
#'                                                          log_abund_detection_thresh=-5,
#'                                                          percent_max_threshold=0,
#'                                                          ...)
#' {
#'
#'   # Get a table of the cell types that are in the control
#'   # FIXME: "knockout" is hard coded and should be a user-defined term in the model
#'   extant_cell_type_df = get_extant_cell_types(control_timeseries_ccm,
#'                                               start_time,
#'                                               stop_time,
#'                                               interval_col=interval_col,
#'                                               percent_max_threshold=percent_max_threshold,
#'                                               log_abund_detection_thresh=log_abund_detection_thresh,
#'                                               perturbation_col=perturbation_col,
#'                                               knockout=FALSE,
#'                                               ...)
#'
#'   # Now let's set up a directed graph that links the states between which cells *could* directly
#'   # transition. If we don't know the direction of flow, add edges in both directions. The idea is
#'   # that we will find shortest paths over this graph between destination states and their plausible
#'   # origin states, then choose the best origins for each destination.
#'
#'   pathfinding_graph = init_pathfinding_graph(control_timeseries_ccm, extant_cell_type_df)
#'
#'   perturbation_ccm_tbl = perturbation_ccm_tbl %>%
#'          dplyr::mutate(perturb_pathfinding_graph = purrr::map2(.f = purrr::possibly(
#'            assess_graph_for_perturbation_support, NA_real_),
#'            .x = perturb_ccm,
#'            .y = perturb_time_window,
#'            pathfinding_graph,
#'            interval_step=interval_step,
#'            interval_col=interval_col,
#'            percent_max_threshold=percent_max_threshold,
#'            log_abund_detection_thresh=log_abund_detection_thresh,
#'            q_val=q_val,
#'            #perturbation_col=perturbation_col,
#'            ...))
#'
#'   perturbation_graph_edges = perturbation_ccm_tbl %>%
#'     select(perturb_name, perturb_pathfinding_graph) %>%
#'     #filter(!is.na(perturb_pathfinding_graph)) %>%
#'     mutate(perturb_pathfinding_graph_df = purrr::map(.f = purrr::possibly(igraph::as_data_frame, NA_real_), .x = perturb_pathfinding_graph)) %>%
#'     filter(!is.na(perturb_pathfinding_graph_df)) %>%
#'     tidyr::unnest(perturb_pathfinding_graph_df)
#'
#'   #mutate(conflicts_with_perturb = to_lost_before_from,
#'   #       supported_by_perturb = from_lost_before_to)
#'
#'   #print (perturbation_graph_edges)
#'
#'   print ("example conflicts:")
#'   perturbation_graph_edges %>% filter(conflicts_with_perturb & from %in% c("16", "22") & to %in% c("16", "22")) %>% select(-perturb_pathfinding_graph) %>% print(n=1000)
#'
#'
#'   print ("example supports:")
#'   perturbation_graph_edges %>% filter(supported_by_perturb & from %in% c("16", "22") & to %in% c("16", "22")) %>% select(-perturb_pathfinding_graph) %>% print(n=1000)
#'
#'
#'   print ("example conflicts:")
#'   perturbation_graph_edges %>% filter(conflicts_with_perturb) %>% print (n=50)
#'
#'   print ("example supports:")
#'   perturbation_graph_edges %>% filter(supported_by_perturb) %>% print (n=50)
#'
#'
#'   perturbation_graph_edge_summary = perturbation_graph_edges %>%
#'     group_by(from, to) %>%
#'     summarize(num_perturbations = length(unique(perturb_name)),
#'               num_perturbations_supporting = sum(supported_by_perturb),
#'               num_perturbations_conflicting = sum(conflicts_with_perturb),
#'               perturbations_supporting=paste0(unique(perturb_name[supported_by_perturb]), collapse = "\n"),
#'               perturbations_conflicting=paste0(unique(perturb_name[conflicts_with_perturb]), collapse = "\n"),
#'               net_perturbations_supporting = num_perturbations_supporting - num_perturbations_conflicting,
#'               perturbation_support_score = num_perturbations + num_perturbations_supporting - num_perturbations_conflicting)
#'
#'   #perturbation_graph_edge_summary = perturbation_graph_edge_summary %>% filter(net_perturbations_supporting >= 0)
#'
#'   #,
#'   #edge_name = stringr::str_c(from, to, sep="~"),
#'   #,
#'   #)
#'
#'   #print(perturbation_ccm_tbl)
#'
#'   #combined_pathfinding_graph = perturbation_ccm_tbl %>% filter(!is.na(perturb_pathfinding_graph)) %>% pull(perturb_pathfinding_graph)
#'   #combined_pathfinding_graph = do.call(igraph::union, combined_pathfinding_graph)
#'
#'   weighted_edges = hooke:::weigh_edges_by_umap_dist(control_timeseries_ccm, perturbation_graph_edge_summary)
#'   combined_pathfinding_graph = igraph::graph_from_data_frame(weighted_edges, directed=TRUE, vertices=data.frame(id=igraph::V(pathfinding_graph)$name))
#'
#'   #   #igraph::edge_attr(G, "support") = igraph::edge_attr(G, "total_perturb_path_score_supporting")
#'   # cycle_breaking_scores = igraph::edge_attr(combined_pathfinding_graph, "num_perturbations_supporting")
#'   # cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0
#'   # igraph::edge_attr(combined_pathfinding_graph, "support") = cycle_breaking_scores
#'   #
#'   # print ("breaking cycles...")
#'   # combined_pathfinding_graph = break_cycles_in_state_transition_graph(combined_pathfinding_graph)
#'   # print ("done")
#'
#'   # print (combined_pathfinding_graph)
#'   # print (igraph::E(combined_pathfinding_graph)$weight)
#'   #
#'   # G = build_timeseries_transition_dag(control_timeseries_ccm,
#'   #                                     extant_cell_type_df,
#'   #                                     combined_pathfinding_graph,
#'   #                                     q_val=q_val,
#'   #                                     start_time=start_time,
#'   #                                     stop_time=stop_time,
#'   #                                     interval_col=interval_col,
#'   #                                     interval_step=interval_step,
#'   #                                     min_interval=min_interval,
#'   #                                     max_interval=max_interval,
#'   #                                     min_dist_vs_time_r_sq=0,
#'   #                                     knockout=FALSE,
#'   #                                     ...)
#'   #
#'   # igraph::V(G)$cell_group = control_timeseries_ccm@ccs@info$cell_group
#'   # return(G)
#'
#'   return(combined_pathfinding_graph)
#' }



#' Assemble a state transition graph from a set of perturbations
#'
#' This function takes as input a control timeseries Hooke model, and a set of timeseries perturbation models.
#' Each perturbation model describes a separate experimental perturbation that may eliminate one or more
#' cell states in the experiment. The tibble must have columns "perturb_name" and "perturb_ccm", and each row
#' must have a perturbation model with a unique name.
#'
#' The function returns a state graph with edges annotated by the level of support from the perturbations.
#' @export
assemble_transition_graph_from_perturbations <- function(control_timeseries_ccm,
                                                         perturbation_ccm_tbl,
                                                         q_val=0.01,
                                                         start_time = NULL,
                                                         stop_time = NULL,
                                                         perturbation_col="knockout",
                                                         interval_col="timepoint",
                                                         interval_step = 2,
                                                         min_interval = 4,
                                                         max_interval = 24,
                                                         log_abund_detection_thresh=-5,
                                                         percent_max_threshold=0,
                                                         allow_links_between_components=FALSE,
                                                         verbose=FALSE,
                                                         ...)
{

  # Get a table of the cell types that are in the control
  # FIXME: "knockout" is hard coded and should be a user-defined term in the model
  extant_cell_type_df = get_extant_cell_types(control_timeseries_ccm,
                                              start_time,
                                              stop_time,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              perturbation_col=perturbation_col,
                                              knockout=FALSE,
                                              ...)

  if (verbose)
    message ("Setting up pathfinding graph")
  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.
  pathfinding_graph = init_pathfinding_graph(control_timeseries_ccm,
                                             extant_cell_type_df,
                                             allow_links_between_components=allow_links_between_components)


  if (verbose)
    message ("Computing paths between cell state losses")

  # Now let's find paths between nodes that are lost following each perturbation. This step
  # returns a tibble of paths, each of which links lost nodes in the pathfinding graph.
  # Paths are also consistent with the flow of time in the control, at least over the window
  # of time measured in the corresponding perturbation experiment. Note that the same path may
  # recur in multiple perturbation experiments. Such paths are returned separately in this table.
  path_tbl = perturbation_ccm_tbl %>%
    dplyr::mutate(paths_between_concordant_loss_nodes = purrr::map2(.f = purrr::possibly(
      get_perturbation_paths, NA_real_),
      .x = perturb_ccm,
      .y = perturb_time_window,
      pathfinding_graph,
      interval_col=interval_col,
      interval_step = interval_step,
      q_val=q_val,
      ...))
  path_tbl = path_tbl %>%
    filter(!is.na(paths_between_concordant_loss_nodes)) %>%
    tidyr::unnest(paths_between_concordant_loss_nodes) %>%
    filter(!is.na(path))

  if (verbose)
    message ("Assessing time flows across loss paths")

  # Measure the flow of time along each perturbation loss path
  path_tbl = path_tbl %>%
    mutate(time_vs_distance_model_stats = purrr::map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                     .x = path,
                                                     ccs=control_timeseries_ccm@ccs,
                                                     interval_col=interval_col)) %>%
    tidyr::unnest(time_vs_distance_model_stats)


  #print(path_tbl %>% head)
  #print (colnames(path_tbl))

  # Compute some scores that summarize the level of support for each path by the various perturbations.
  path_tbl = path_tbl %>% mutate(timeseries_model_path_score = ifelse(time_dist_effect > 0 & time_dist_effect_pval < q_val, time_dist_model_adj_rsq, 0))
  path_tbl = path_tbl %>% mutate(perturb_model_path_score = ifelse(perturb_dist_effect < 0 & perturb_dist_effect_pval < q_val, perturb_dist_model_adj_rsq, 0))
  path_tbl = path_tbl %>% mutate(path_score = ifelse(timeseries_model_path_score > 0 & perturb_model_path_score > 0, perturb_model_path_score * timeseries_model_path_score, 0))

  # Exclude paths that have no support from perturbations or go against the flow of time
  path_tbl =  path_tbl %>% filter (path_score > 0)

  if (verbose)
    message ("Constructing transition graph")

  # Now let's start to build up a state transition graph from the paths by taking their union across all
  # perturbations, provided they survived the above filtering steps.
  selected_paths = path_tbl %>%
    #group_by(to) %>%
    ungroup() %>% arrange(desc(path_score)) %>%
    dplyr::rename(path_contrast=perturb_name)
  G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = TRUE)

  # Now go back and score each edge in the state graph for support from perturbations, as well as support
  # the full control timeseries
  if (verbose)
    message ("Assessing perturbation support for transition graph")
  G = assess_support_for_transition_graph(control_timeseries_ccm,
                                          perturbation_ccm_tbl,
                                          path_tbl,
                                          G,
                                                    q_val,
                                                    start_time,
                                                    stop_time,
                                                    perturbation_col,
                                                    interval_col,
                                                    interval_step,
                                                    min_interval,
                                                    max_interval,
                                                    log_abund_detection_thresh,
                                                    percent_max_threshold,
                                                    ...)

  # Now that we have a graph in which every edge is supported by at least one perturbations,
  # let's add edges that may not have perturbation support but do have timeseries support
  # and also do not conflict with the perturbation-based graph. First, we will compute
  # paths between states that are reciprocally gained and lost over the timeseries
  # (i.e. one is lost and the other is gained)


  if (verbose)
    message ("Computing timeseries-only transition graph")

  paths_between_recip_time_nodes = get_timeseries_paths(control_timeseries_ccm,
                                                        extant_cell_type_df,
                                                        pathfinding_graph,
                                                        start_time = start_time,
                                                        stop_time = stop_time,
                                                        interval_col=interval_col,
                                                        interval_step = interval_step,
                                                        min_interval = min_interval,
                                                        max_interval = max_interval,
                                                        knockout=FALSE,
                                                        ...)

  # Sort these paths by the level of support from the timeseries, as quantified by how well
  # distance along a path correlates with an increase in the time cells that fall on it were
  # collected.
  selected_paths = paths_between_recip_time_nodes %>% filter (time_dist_effect > 0 & time_dist_effect_pval < q_val) %>%
    ungroup() %>% arrange(desc(time_dist_model_adj_rsq)) %>%
    mutate(path_contrast="time")

  #G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles=FALSE)

  # Add these paths, one at a time, in order of best score, provided that they do not add new
  # cycles to the graph.
  num_cycles_in_orig_G = length(find_cycles(G))
  for (i in (1:nrow(selected_paths))){
    next_path = selected_paths$path[[i]] %>% select(from, to)
    next_path_graph = next_path %>% igraph::graph_from_data_frame(directed=TRUE)
    G_prime = igraph::union(G, next_path_graph)
    if (length(find_cycles(G_prime)) == num_cycles_in_orig_G){
      G = G_prime
    }else{
      # Debug:
      #print ("skipping due to cycle:")
      #print (selected_paths[i,] %>% select(-path))
    }
  }

  # FIXME: this is gross and there is probably a cleaner way, but what we're doing here
  # is annotating that these edges from the timeseries are not supported by the perturbations
  G_num_perturbs_supporting = igraph::edge_attr(G, "num_perturbs_supporting")
  G_num_perturbs_supporting[is.na(G_num_perturbs_supporting)] = 0
  igraph::edge_attr(G, "num_perturbs_supporting") = G_num_perturbs_supporting

  G_max_timeseries_path_score_supporting = igraph::edge_attr(G, "max_timeseries_path_score_supporting")
  G_max_timeseries_path_score_supporting[is.na(G_max_timeseries_path_score_supporting)] = 0
  igraph::edge_attr(G, "max_timeseries_path_score_supporting") = G_max_timeseries_path_score_supporting

  G_total_timeseries_path_score_supporting = igraph::edge_attr(G, "total_timeseries_path_score_supporting")
  G_total_timeseries_path_score_supporting[is.na(G_total_timeseries_path_score_supporting)] = 0
  igraph::edge_attr(G, "total_timeseries_path_score_supporting") = G_total_timeseries_path_score_supporting

  G_max_perturb_path_score_supporting = igraph::edge_attr(G, "max_perturb_path_score_supporting")
  G_max_perturb_path_score_supporting[is.na(G_max_perturb_path_score_supporting)] = 0
  igraph::edge_attr(G, "max_perturb_path_score_supporting") = G_max_perturb_path_score_supporting

  G_total_perturb_path_score_supporting = igraph::edge_attr(G, "total_perturb_path_score_supporting")
  G_total_perturb_path_score_supporting[is.na(G_total_perturb_path_score_supporting)] = 0
  igraph::edge_attr(G, "total_perturb_path_score_supporting") = G_total_perturb_path_score_supporting

  G_max_path_score_supporting = igraph::edge_attr(G, "max_path_score_supporting")
  G_max_path_score_supporting[is.na(G_max_path_score_supporting)] = 0
  igraph::edge_attr(G, "max_path_score_supporting") = G_max_path_score_supporting

  G_total_path_score_supporting = igraph::edge_attr(G, "total_path_score_supporting")
  G_total_path_score_supporting[is.na(G_total_path_score_supporting)] = 0
  igraph::edge_attr(G, "total_path_score_supporting") = G_total_path_score_supporting

  G_label = igraph::edge_attr(G, "label")
  G_label[is.na(G_label)] = ""
  igraph::edge_attr(G, "label") = G_label

  return (G)
  # print ("done")
  #
  # # cycle_breaking_scores = igraph::edge_attr(G, "max_perturb_score_supporting")
  # # cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0
  # # igraph::edge_attr(G, "support") = cycle_breaking_scores
  #
  # print ("after assessment")
  # #G %>% igraph::as_data_frame() %>% filter(from == "27" & to == "24")  %>% print(n=1000)
  # # print (G %>% igraph::as_data_frame() )
  # # missing_pathfinding_edges = dplyr::setdiff(as_tibble(igraph::as_edgelist(timeseries_pathfinding_graph)), as_tibble(igraph::as_edgelist(G)))
  # # colnames(missing_pathfinding_edges) = c("from", "to")
  # # G = G + igraph::graph_from_edgelist(as.matrix(missing_pathfinding_edges))
  #
  # #igraph::edge_attr(G, "support") = igraph::edge_attr(G, "total_perturb_path_score_supporting")
  # cycle_breaking_scores = igraph::edge_attr(G, "max_perturb_score_supporting")
  # cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0
  # igraph::edge_attr(G, "support") = cycle_breaking_scores
  #
  # print (igraph::edge_attr(G, "support"))
  #
  # print ("breaking cycles...")
  # G = break_cycles_in_state_transition_graph(G)
  # print ("done")
  #
  # # print ("finding min path cover...")
  # # covered_G = compute_min_path_cover(control_timeseries_ccm, G, weight_attribute="support")
  # # print ("done. path covered graph:")
  #
  # #annotated_G = covered_G #%>% igraph::as_data_frame()
  # annotated_G = G
  # igraph::V(annotated_G)$cell_group = control_timeseries_ccm@ccs@info$cell_group
  #
  # print (annotated_G)
  # annotated_G = assess_support_for_transition_graph(control_timeseries_ccm,
  #                                         perturbation_ccm_tbl,
  #                                         path_tbl,
  #                                         annotated_G,
  #                                         q_val,
  #                                         start_time,
  #                                         stop_time,
  #                                         perturbation_col,
  #                                         interval_col,
  #                                         interval_step,
  #                                         min_interval,
  #                                         max_interval,
  #                                         log_abund_detection_thresh,
  #                                         percent_max_threshold,
  #                                         ...)

  return(annotated_G)

}

#' Assess support for a graph built via perturbations
#' @export
assess_support_for_transition_graph <- function(control_timeseries_ccm,
                                                perturbation_ccm_tbl,
                                                perturbation_path_tbl,
                                                state_transition_graph,
                                                         q_val=0.01,
                                                         start_time = NULL,
                                                         stop_time = NULL,
                                                         perturbation_col="knockout",
                                                         interval_col="timepoint",
                                                         interval_step = 2,
                                                         min_interval = 4,
                                                         max_interval = 24,
                                                         log_abund_detection_thresh=-5,
                                                         percent_max_threshold=0,
                                                         ...)
{
   path_score_tbl = perturbation_path_tbl %>% rename(origin=from, destination=to) %>% tidyr::unnest(path)

  #print (path_score_tbl)

  #edge_support = left_join(edge_support, path_score_tbl)

  edge_support_summary = path_score_tbl %>%
    ungroup () %>%
    dplyr::select(from, to, perturb_name,
                  timeseries_model_path_score,
                  perturb_model_path_score,
                  path_score)  %>%
    group_by(from, to, perturb_name) %>%
    summarize(timeseries_model_path_score = max(timeseries_model_path_score),
              perturb_model_path_score = max(perturb_model_path_score),
              path_score = max(path_score)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(from, to) %>%
    summarize(num_perturbs_supporting=length(unique(perturb_name)),
              #num_perturb_intervals_supporting=sum(num_intervals_supported),
              max_timeseries_path_score_supporting=max(timeseries_model_path_score),
              total_timeseries_path_score_supporting=sum(timeseries_model_path_score),
              max_perturb_path_score_supporting=max(perturb_model_path_score),
              total_perturb_path_score_supporting=sum(perturb_model_path_score),
              max_path_score_supporting = max(path_score),
              total_path_score_supporting=sum(path_score))

  edge_support_labels = path_score_tbl %>%
    ungroup () %>%
    dplyr::select(from, to, perturb_name)  %>%
    group_by(from, to) %>%
    distinct() %>%
    summarize(edge_name = stringr::str_c(from, to, sep="~"),
              label=paste0(perturb_name, collapse = "\n"))

  #print (edge_support_labels)
  edge_support_summary = edge_support_summary %>% left_join(edge_support_labels)

  #edge_support_summary = edge_support_summary %>% mutate(support_weight = ifelse(is.na(support_weight), 0, support_weight))
  edge_support_summary = edge_support_summary %>% dplyr::distinct()
  #print ("re-nesting")
  #edge_support = edge_support %>%
    #dplyr::select(-distance_from_root, -weight) %>%
    #tidyr::nest(support=c(perturb_name, perturb_dist_effect, perturb_dist_effect_q_val))
  #  tidyr::nest(support=c(perturb_name))

  #print (edge_support_summary)
  #print ("annotating the graph")
  annotated_state_transition_graph = igraph::as_data_frame(state_transition_graph)  %>% as_tibble %>% dplyr::select(from, to) %>% distinct #%>% left_join(edge_support)
  annotated_state_transition_graph = annotated_state_transition_graph %>% left_join(edge_support_summary)
  #annotated_state_transition_graph = annotated_state_transition_graph %>% mutate(support_weight = ifelse(is.na(support_weight), 0, support_weight))
  annotated_state_transition_graph = igraph::graph_from_data_frame(annotated_state_transition_graph, directed=TRUE, vertices=data.frame(id=igraph::V(state_transition_graph)$name))
  return(annotated_state_transition_graph)
}

#' Simplify a directed state transition graph by grouping nodes according to a specified label
#' @export
contract_state_graph <- function(ccm, state_graph, group_nodes_by){
  # Create simplified cell state graph just on cell type (not cluster):
  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs@cds) %>%
    as.data.frame %>% select(!!sym(group_nodes_by))

  cell_group_metadata$cell_group = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
    as.data.frame %>%
    dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
    dplyr::group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
  colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
  node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))

  contraction_mapping = as.factor(node_metadata$group_nodes_by)
  contraction_mapping_names = as.character(levels(contraction_mapping))
  contraction_mapping = as.numeric(contraction_mapping)
  names(contraction_mapping) = node_metadata$id
  contracted_state_graph = igraph::contract(state_graph, mapping=contraction_mapping, vertex.attr.comb="ignore")
  igraph::V(contracted_state_graph)$name = unlist(contraction_mapping_names[as.numeric(igraph::V(contracted_state_graph))])
  contracted_state_graph = igraph::simplify(contracted_state_graph)
  igraph::V(contracted_state_graph)$cell_group = group_nodes_by
  return(contracted_state_graph)
}



# augment graph functions -- still need some work

#'@param ccm
#'@param sparsity
#'@param not_in_paga_and_emerge
find_edge_to_each <- function(ccm, sparsity, not_in_paga_and_emerge) {

  ccm = select_model(ccm, criterion = "EBIC", sparsity_factor = sparsity)

  not_in_paga_cov_edges =  hooke:::return_igraph(model(ccm, "reduced")) %>%
    igraph::as_data_frame(what="edges") %>%
    dplyr::rename(pcor=weight) %>% dplyr::filter(pcor != 0.00) %>%
    filter(from %in% not_in_paga_and_emerge | to %in% not_in_paga_and_emerge)

  edge_to_each = length(setdiff(not_in_paga_and_emerge,
                                union(not_in_paga_cov_edges$from, not_in_paga_cov_edges$to))) == 0

  return(edge_to_each)

}

#' @param ccm a cell count model object
#' @param extant_cell_type_df output from get_extant_cell_types()
#' @param step_size of sparsity values to try
#' @param n the number of values to try
determine_sparsity <- function(ccm,
                               extant_cell_type_df,
                               start_time,
                               stop_time,
                               step_size = 0.01,
                               n = 10) {

  # maybe make this so that it can calculate extant cell type df on own

  # find ones not in paga that emerge after the start
  not_in_paga = not_in_paga_graph(ccm)
  initial_sparsity = ccm@sparsity

  not_in_paga_and_emerge = extant_cell_type_df %>%
    filter(cell_group %in% not_in_paga) %>%
    filter(longest_contig_start > start_time) %>%
    pull(cell_group) %>% unique()


  # first decide to go up or down
  edge_to_each = find_edge_to_each(ccm, initial_sparsity, not_in_paga_and_emerge)

  # if there isn't already an edge to each, you can add the step size instead
  if (edge_to_each) {
    step_size = -1*step_size
  }

  df = data.frame("step" = (seq(1:n)-1)*step_size) %>%
    mutate("sparsity" = initial_sparsity-step) %>%
    filter(sparsity > 0)

  df = df %>% mutate(to_each = purrr::map(.f = find_edge_to_each,
                                          .x = sparsity,
                                          ccm = ccm,
                                          not_in_paga_and_emerge = not_in_paga_and_emerge))

  selected_sparsity = df %>% filter(to_each == TRUE) %>% pull(sparsity) %>% max()

  # if don't find any, suggest trying more
  if (nrow(df %>% filter(to_each == TRUE)) == 0) {
    #print(paste0("try more n starting at ", df[n,]$sparsity))
    ccm = select_model(ccm, criterion = "EBIC", sparsity_factor = df[n,]$sparsity)
    return(ccm)
  } else {
    ccm = select_model(ccm, criterion = "EBIC", sparsity_factor = selected_sparsity)
    return(ccm)
  }
}

#' This function takes a state transition graph and adds links between disjoint groups of cells
#' @param ccm
#' @param curr_graph
#' @export
augment_pathfinding_graph = function(ccm,
                                     curr_graph,
                                     q_val = 0.01,
                                     start_time = NULL,
                                     stop_time = NULL,
                                     interval_col="timepoint",
                                     interval_step = 2,
                                     min_interval = 4,
                                     max_interval = 24,
                                     log_abund_detection_thresh=-5,
                                     percent_max_threshold=0,
                                     min_dist_vs_time_r_sq=0,
                                     ...) {

  # start with the same pathfinding graph as before
  # make sure to do this before so it doesnt change the previous sparsity??
  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start_time,
                                              stop_time,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              ...)

  # get the initial path
  timeseries_pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)

  # choose a good sparsity
  ccm = determine_sparsity(ccm, extant_cell_type_df, start_time, stop_time)

  # find ones not in paga that emerge after the start
  not_in_paga = not_in_paga_graph(ccm)
  initial_sparsity = ccm@sparsity

  not_in_paga_and_emerge = extant_cell_type_df %>%
    filter(cell_group %in% not_in_paga) %>%
    filter(longest_contig_start > start_time) %>%
    pull(cell_group) %>% unique()

  not_in_paga_cov_edges =  hooke:::return_igraph(model(ccm, "reduced")) %>%
    igraph::as_data_frame(what="edges") %>%
    dplyr::rename(pcor=weight) %>% dplyr::filter(pcor != 0.00) %>%
    filter(from %in% not_in_paga_and_emerge | to %in% not_in_paga_and_emerge)

  weighted_edges = hooke:::weigh_edges_by_umap_dist(ccm, not_in_paga_cov_edges)
  node_metadata = data.frame(id=unique(union(weighted_edges$to, weighted_edges$from)))

  pathfinding_graph = weighted_edges %>% select(from, to, weight) %>%
    igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
    igraph::as.directed()

  edges_between_concurrent_states = left_join(pathfinding_graph %>% igraph::as_data_frame(what="edges"),
                                              extant_cell_type_df, by=c("from"="cell_group")) %>%
    select(from, to, from_start=longest_contig_start, from_end=longest_contig_end)
  edges_between_concurrent_states = left_join(edges_between_concurrent_states, extant_cell_type_df, by=c("to"="cell_group")) %>%
    select(from, to, from_start, from_end, to_start=longest_contig_start, to_end=longest_contig_end)
  edges_between_concurrent_states = edges_between_concurrent_states %>%
    filter(from_start <= to_start & from_end >= to_start) %>% select(from, to) %>% distinct()

  pathfinding_graph = igraph::intersection(pathfinding_graph,
                                           edges_between_concurrent_states %>%
                                             igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata))

  # only take the top options for each
  # pathfinding_graph =
  #   igraph::as_data_frame(pathfinding_graph) %>%
  #   mutate(node = ifelse(to %in% not_in_paga_and_emerge, to, from)) %>%
  #   group_by(node) %>% top_n(1, wt = -weight) %>%
  #   add them to the current pathfinding graph
  #   rbind(igraph::as_data_frame(timeseries_pathfinding_graph)) %>%
  #   igraph::graph_from_data_frame()

  # undebug(build_transition_dag)
  G = build_transition_dag(ccm,
                           extant_cell_type_df,
                           pathfinding_graph,
                           q_val,
                           start_time,
                           stop_time,
                           interval_col,
                           interval_step,
                           min_interval,
                           max_interval,
                           min_dist_vs_time_r_sq,
                           ...)

  # take the lowest weight edges
  from_nodes = igraph::as_data_frame(G) %>%
    filter(from %in% not_in_paga_and_emerge) %>%
    group_by(from) %>% top_n(1, wt = -weight)
  to_nodes = igraph::as_data_frame(G) %>%
    filter(to %in% not_in_paga_and_emerge) %>%
    group_by(to) %>% top_n(1, wt = -weight)

  new_G = rbind(from_nodes, to_nodes) %>% igraph::graph_from_data_frame()

  covered_G = compute_min_path_cover(ccm, new_G)

  # add any new edges
  new_edges = igraph::difference(covered_G, curr_graph) %>%
    igraph::as_data_frame() %>%
    filter(from %in% not_in_paga_and_emerge | to %in% not_in_paga_and_emerge) %>%
    igraph::graph_from_data_frame()

  aug_graph = igraph::union(curr_graph, new_edges)

  return(aug_graph)

}
