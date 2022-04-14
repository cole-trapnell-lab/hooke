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
get_paga_graph <- function(cds, reduction_method = "UMAP") {

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
                                                     0.05, FALSE)


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
                                  ...){

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step=interval_step, ...)
  timepoint_pred_df = timepoint_pred_df %>%
    group_by(cell_group) %>%
    mutate(max_abundance = max(exp(log_abund)),
           percent_max_abund = exp(log_abund) / max_abundance,
           above_log_abund_thresh = log_abund - 2*log_abund_se > log_abund_detection_thresh,
           above_percent_max_thresh = percent_max_abund > percent_max_threshold,
           present_flag = ifelse(above_log_abund_thresh & above_percent_max_thresh, TRUE, NA)) %>% ungroup()

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

  undebug(longest_present_interval)
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


find_paths_to_origins <- function(ccm,
                                  neg_rec_edges_to_destinations,
                                  abund_over_time_df,
                                  traversal_graph,
                                  interval_col="timepoint",
                                  q_value_threshold = 1.0,
                                  require_presence_at_all_timepoints = TRUE) {

  interval_start = min(neg_rec_edges_to_destinations[,paste("to", interval_col, "x", sep="_")])
  interval_stop =  max(neg_rec_edges_to_destinations[,paste("to", interval_col, "y", sep="_")])
  cell_types_exant_in_contrast = abund_over_time_df %>% filter(!!sym(interval_col) >= interval_start &
                                                                 !!sym(interval_col) <= interval_stop)

  # # These are the cell groups that are missing in at least one timepoint in the interval we're looking at:
  # if (require_presence_at_all_timepoints){
  #   cell_groups_missing_in_range = cell_types_exant_in_contrast %>% filter(present_above_thresh == FALSE) %>% pull(cell_group) %>% unique
  #   cell_groups_missing_in_range = union(setdiff(row.names(ccm@ccs), cell_types_exant_in_contrast$cell_group), cell_groups_missing_in_range)
  # }else{
  #   cell_groups_present_in_range = cell_types_exant_in_contrast %>% filter(present_above_thresh == TRUE) %>% pull(cell_group) %>% unique
  #   cell_groups_missing_in_range = setdiff(row.names(ccm@ccs), cell_groups_present_in_range)
  # }
  #
  # traversal_graph = igraph::delete_vertices(traversal_graph, cell_groups_missing_in_range)

  cell_groups_present_in_range = cell_types_exant_in_contrast %>% filter(present_above_thresh == TRUE) %>% pull(cell_group) %>% unique
  cell_groups_missing_in_range = setdiff(row.names(ccm@ccs), cell_groups_present_in_range)
  traversal_graph = igraph::delete_vertices(traversal_graph, cell_groups_missing_in_range)

  edge_path = neg_rec_edges_to_destinations %>%
    dplyr::mutate(shortest_path = purrr::map2(.f =
                                                purrr::possibly(hooke:::get_shortest_path, NA_real_),
                                              .x = from, .y = to,
                                              traversal_graph)) %>%
    select(origin=from, destination=to, shortest_path) %>%
    tidyr::unnest(shortest_path) %>%
    filter(is.na(from) == FALSE & is.na(to) == FALSE) %>%
    # select(-weight) %>%
    distinct()

  return(edge_path)
}
#undebug(find_paths_to_origins)

#
# xxx_paths = find_paths_to_origins(wt_ccm_wl,
#                          compare_abundances(wt_ccm_wl,
#                                             estimate_abundances(wt_ccm_wl, tibble(expt="GAP16", timepoint=40)),
#                                             estimate_abundances(wt_ccm_wl, tibble(expt="GAP16", timepoint=72))),
#                          get_extant_cell_types(wt_ccm_wl, 18, 96, log_abund_detection_thresh = -5, percent_max_threshold=0.01),
#                          paga_edges,
#                          require_presence_at_all_timepoints = TRUE)

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

#' Function that computes a bunch of statistics on possible origin states
compute_origin_stats <- function(ccm, tp_pred_df, ext_ct_df){
  pf_gr = init_pathfinding_graph(ccm, ext_ct_df)

  #tp_pred_df = estimate_abundances_over_interval(ccm, 18, 24, interval_col="timepoint", experiment="GAP14")

  reduced_model_pcors_between_cell_groups = model(ccm, "reduced")$latent_network(type = "partial_cor") %>% as.matrix() %>% as.data.frame.table(responseName = "reduced_pcor")
  full_model_pcors_between_cell_groups = model(ccm, "full")$latent_network(type = "partial_cor") %>% as.matrix() %>% as.data.frame.table(responseName = "full_pcor")

  geodesic_distances_between_cell_groups = igraph::distances(pf_gr, mode="out")  %>% as.data.frame.table(responseName = "geodesic_dist")

  cell_group_centroids = centroids(ccm@ccs)
  euclidean_distances_between_cell_groups = as.matrix(dist(cell_group_centroids[,-1], method = "euclidean", diag = T)) %>% as.data.frame.table(responseName = "euclidean_dist")

  pcor_dist_df = left_join(geodesic_distances_between_cell_groups, euclidean_distances_between_cell_groups)
  pcor_dist_df = left_join(pcor_dist_df, reduced_model_pcors_between_cell_groups)
  pcor_dist_df = left_join(pcor_dist_df, full_model_pcors_between_cell_groups)

  pcor_dist_df = pcor_dist_df %>% as_tibble()
  pcor_dist_df = pcor_dist_df %>% dplyr::rename(from=Var1, to=Var2) %>% mutate(from = as.character(from),
                                                                               to = as.character(to))
  pcor_dist_df = pcor_dist_df %>% filter(from != to)

  pcor_dist_df = pcor_dist_df %>% filter(is.finite(geodesic_dist))
  pcor_dist_df = pcor_dist_df %>% filter(reduced_pcor != 0)
  #pcor_dist_df = pcor_dist_df %>% filter(full_pcor != 0)


  pcor_dist_df = pcor_dist_df %>% tidyr::nest(data=-from)

  pcor_dist_models = pcor_dist_df %>%
    mutate(dist_model = purrr::map(data, ~ lm(reduced_pcor ~ geodesic_dist, data=.x))) %>%
    mutate(model_summary = purrr::map(dist_model, broom::tidy)) #%>%

  origin_stats = pcor_dist_models %>%
    tidyr::unnest(model_summary) %>%
    dplyr::filter(term == "geodesic_dist") %>%
    dplyr::select(from, dist_on_pcor_effect=estimate, dist_effect_pval=p.value)

  origin_stats = left_join(origin_stats, geodesic_distances_between_cell_groups, by=c("from"="Var1")) %>% dplyr::rename(to=Var2)
  origin_stats = left_join(origin_stats, reduced_model_pcors_between_cell_groups, by=c("from"="Var1", "to"="Var2"))
  origin_stats = left_join(origin_stats, full_model_pcors_between_cell_groups, by=c("from"="Var1", "to"="Var2"))
  origin_stats = origin_stats %>% dplyr::select(from,
                                                to,
                                                dist_on_pcor_effect,
                                                dist_effect_pval,
                                                distance=geodesic_dist,
                                                reduced_pcor,
                                                full_pcor)
  return (origin_stats)

  #valid_origins = pcor_dist_models %>%
  #  tidyr::unnest(model_summary) %>%
  #  filter(term == "geodesic_dist" & estimate < 0) %>%
  #  pull(from)

}

# Helper function
find_origin_paths_from_time_contrasts <- function(ccm,
                                                  traversal_graph,
                                                  time_contrasts,
                                                  timepoint_pred_df,
                                                  extant_cell_type_df,
                                                  interval_col,
                                                  q_val = 0.01,
                                                  require_presence_at_all_timepoints=TRUE){

  select_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
    cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
    cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  time_contrasts = time_contrasts %>% mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                                                      .x = t1,
                                                                      .y = t2,
                                                                      interval_col=interval_col,
                                                                      timepoint_pred_df = timepoint_pred_df))

  time_contrasts = time_contrasts %>% mutate(neg_rec_edges = purrr::map(.f = hooke:::get_neg_dir_edges,
                                                                        .x = comp_abund,
                                                                        ccm = ccm,
                                                                        q_value_threshold = q_val))

  time_contrasts = time_contrasts%>%
    mutate(destinations_that_need_origins = purrr::map2(.f = get_destinations_that_need_origins,
                                                        .x = neg_rec_edges,
                                                        .y = t2,
                                                        interval_col=interval_col,
                                                        extant_cell_type_df)) %>%
    mutate(neg_rec_for_destinations = purrr::map2(.f = get_neg_rec_edges_for_cells_that_need_origins,
                                                  .x = neg_rec_edges,
                                                  .y = destinations_that_need_origins))
  time_contrasts = time_contrasts %>%
    mutate(path = purrr::map(.f = purrr::possibly(find_paths_to_origins, NA_real_),
                             .x = neg_rec_for_destinations,
                             ccm = ccm,
                             traversal_graph=traversal_graph,
                             interval_col=interval_col,
                             abund_over_time_df = extant_cell_type_df,
                             require_presence_at_all_timepoints=require_presence_at_all_timepoints,
                             q_value_threshold = 1.0))
  return (time_contrasts)
}


#' Select the origins for each destination and report a set of paths between them
#'
#' @export
select_timeseries_origins <- function(possible_origins,
                                         selection_policy=c("closest-origin",
                                                            "best-origin",
                                                            "max-score-origin",
                                                            "max-score-dist-ratio-origin",
                                                            "all-origins")){
  possible_origins = possible_origins %>% group_by(cell_group)

  possible_origins = possible_origins %>% mutate(origin_score = reduced_pcor,
                                                 origin_score_dist_ratio = origin_score / distance)

  if(selection_policy == "max-score-origin"){
    possible_origins = possible_origins %>% slice_min(origin_score, n=1)
  }else if(selection_policy == "max-score-dist-ratio-origin"){
    possible_origins = possible_origins %>% filter (is.finite(distance)) %>% slice_min(origin_score_dist_ratio, n=1)
  }
  else if(selection_policy == "closest-origin"){
    possible_origins = possible_origins %>% filter (is.finite(distance)) %>% slice_min(distance, n=1)
  }else if (selection_policy == "all-origins"){

  }
  possible_origins = possible_origins %>% ungroup()
  return (possible_origins)
}

#' Identify the possible origins for each destination
#'
#' @export
find_timeseries_origins <- function(ccm,
                                       extant_cell_type_df,
                                       timeseries_pathfinding_graph,
                                       q_val=0.01,
                                       start = NULL,
                                       stop = NULL,
                                       interval_col="timepoint",
                                       interval_step = 2,
                                       min_interval = 4,
                                       max_interval = 24,
                                       percent_max_threshold=0.0,
                                       log_abund_detection_thresh=-5,
                                       require_presence_at_all_timepoints=TRUE,
                                       initial_origin_policy="best-origin",
                                       ...) {

  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start)){
    start = min(colData(ccm@ccs)[,interval_col])
  }
  if (is.null(stop)){
    stop = max(colData(ccm@ccs)[,interval_col])
  }

  timepoints = seq(start, stop, interval_step)

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step, ...)


  # Now let's iterate over the cell types, finding origins for all the ones that emerge
  # during the time interval. We'll do this by first collecting all the pairwise contrasts between
  # timepoints the caller wants. Then, we'll sift through those results for each cell type that needs an
  # origin and pull out the contrasts where it's going up. Then we'll look at the cell types that are
  # going down in that contrast and have a negative pcor. Each of those is a possible origin for the
  # destination.

  destination_cell_types = extant_cell_type_df %>% select(cell_group, emerges_at=longest_contig_start) %>% distinct()

  origin_stats_df = compute_origin_stats(ccm, timepoint_pred_df, extant_cell_type_df)

  time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval)

  get_possible_origins_for_destination <- function(destination,
                                                   emerges_at,
                                                   ccm,
                                                   origin_stats_df,
                                                   time_contrasts,
                                                   timepoint_pred_df,
                                                   q_val,
                                                   origin_dist_effect_thresh=0){

    select_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
      cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
      cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
      return(compare_abundances(ccm, cond_x, cond_y))
    }

    # Let's subset the contrasts to include only those between which the destination emerges
    time_contrasts = time_contrasts %>% filter (t1 <= emerges_at & t2 >= emerges_at)

    if (nrow(time_contrasts) == 0){
      return (NULL)
    }

    # Now let's actually compare the abundances between the endpoints of all those contrasts
    # FIXME: we should be able to relax this to look at edges with positive pcor too, instead of just negative
    # provided that one of the endpoints is the destination.
    relevant_comparisons = time_contrasts %>%
      mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                      .x = t1,
                                      .y = t2,
                                      interval_col=interval_col,
                                      timepoint_pred_df = timepoint_pred_df)) %>%
      mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                        .x = comp_abund,
                                        ccm = ccm))

    # Could add additional constraints here on origins if desired?
    valid_origins = origin_stats_df %>% filter(dist_on_pcor_effect < origin_dist_effect_thresh) %>% pull(from)

    relevant_comparisons = relevant_comparisons %>%
      tidyr::unnest(rec_edges) %>%
      dplyr::filter((to %in% valid_origins & from == destination & from_delta_log_abund > 0 & to_delta_log_abund < 0) |
                    (from %in% valid_origins & to == destination & to_delta_log_abund > 0 & from_delta_log_abund < 0))
    #possible_origins = relevant_comparisons %>% dplyr::select(origin=from, origin_score = pcor) %>% distinct()
    possible_origins = relevant_comparisons %>% dplyr::select(from, to)

    backwards_edges = possible_origins %>% dplyr::filter(from == destination)
    ok_edges = possible_origins %>% dplyr::filter(to == destination)

    be_colnames = colnames(backwards_edges)
    be_colnames = unlist(lapply(be_colnames, stringr::str_replace, "^to", "from"))

    from_cols = grepl("^from", colnames(backwards_edges))
    be_colnames[from_cols] = unlist(lapply(be_colnames[from_cols], stringr::str_replace, "^from", "to"))
    colnames(backwards_edges) = be_colnames

    if (nrow(backwards_edges) != 0) {
      #backwards_edges$edge_type = "directed_from_to"
      backwards_edges = backwards_edges[colnames(ok_edges)]
    }

    possile_origins = rbind(ok_edges, backwards_edges)

    possible_origins = possible_origins %>% left_join(origin_stats_df, by=c("from", "to"))
    possible_origins = possible_origins %>% dplyr::select(-to) %>% dplyr::rename(origin=from)
    return (possible_origins)
  }
  #debug(get_possible_origins_for_destination)

  possible_origins = destination_cell_types %>%
    mutate(possible_origins = purrr::map2(.f = get_possible_origins_for_destination,
                                          .x = cell_group,
                                          .y = emerges_at,
                                          ccm,
                                          origin_stats_df,
                                          time_contrasts,
                                          timepoint_pred_df,
                                          q_val)) %>%
    tidyr::unnest(possible_origins)

  possible_origins = possible_origins %>% distinct()

  return (possible_origins)
}


#' Initialize a graph over which cells can transition
#'
#' @export
init_pathfinding_graph <- function(ccm, extant_cell_type_df){

  # There are a number of different ways we could set up this "pathfinding graph" but for now
  # Let's just use the PAGA (weighed by distance in UMAP space), subtracting edges between which
  # there is zero partial correlation in the nuisance cell count model.

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = data.frame(id=cell_groups)

  paga_graph = initial_pcor_graph(ccm@ccs) %>% igraph::graph_from_data_frame(directed = FALSE, vertices=node_metadata) %>% igraph::as.directed()
  cov_graph = hooke:::return_igraph(model(ccm, "reduced"))
  cov_graph_edges = igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight != 0.00)

  # FIXME: While we should consider this more sophisticated weighting function, let's test just using UMAP distance for now:
  #weighted_edges = hooke:::get_weighted_edges(ccm, all_edges)
  weighted_edges = hooke:::weigh_edges_by_umap_dist(ccm, cov_graph_edges)

  weighted_edges = weighted_edges %>% group_by(from, to) %>%
    mutate(adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to)) %>% ungroup()
  weighted_edges = weighted_edges %>% filter(adjacent_in_paga)

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

  #if (require_presence_at_all_timepoints){
  pathfinding_graph = igraph::intersection(pathfinding_graph,
                                           edges_between_concurrent_states %>%
                                             igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata))
  #}
  return (pathfinding_graph)
}

#debug(find_origins)

find_paths_to_origins <- function(possible_origins, timeseries_pathfinding_graph){
  possible_origins = possible_origins %>%
    mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                              .x = origin, .y = cell_group,
                              timeseries_pathfinding_graph)) %>%
    tidyr::unnest(path) %>%
    filter(is.na(from) == FALSE & is.na(to) == FALSE)
  return (possible_origins)
}

#' A function to assemble a state transition graph from a timeseries model
#' @export
assemble_timeseries_transitions <- function(ccm,
                                            q_val=0.01,
                                            start = NULL,
                                            stop = NULL,
                                            interval_col="timepoint",
                                            interval_step = 2,
                                            min_interval = 4,
                                            max_interval = 24,
                                            percent_max_threshold=0.0,
                                            log_abund_detection_thresh=-5,
                                            initial_origin_policy="closest-origin",
                                            ...){

  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start,
                                              stop,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              ...)
  timeseries_pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)


  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.


  possible_origins = find_timeseries_origins(ccm,
                                                extant_cell_type_df,
                                                timeseries_pathfinding_graph,
                                                start=start,
                                                stop=stop,
                                                interval_col=interval_col,
                                                min_interval = min_interval,
                                                percent_max_threshold=percent_max_threshold,
                                                log_abund_detection_thresh=log_abund_detection_thresh,
                                                #percent_max_threshold=0.00,

                                                ...)


  # Now we have to actually choose the right origins for each cell type. This is hard because if you
  # get them wrong you create weird backwards flows in the trajectory or undesirable shortcuts.

  # There are a bunch of policies we could use to pick the best origin(s) for a given cell state
  # The most conservative is to take the closest one or the one that has the highest ratio of pcor:geodesic
  # distance in the pathfinding graph or over the UMAP manifold. These two are very similar. Moreover, paths
  # from cell types that emerge during the time series are way easier to identify than those that are present
  # at the outset. I think we should conservatively start with the policy "max-score-dist-ratio-origin"
  # only assigning origins for emergent cell types. We could then use these paths to orient edges within the
  # pathfinding graph and then run a second round of origin finding on this graph, and assigning origins using
  # more relaxed policy (e.g. allowing multiple origins for an emergent cell type or allowing origins for cell
  # types that were present at the outset). We could also explore using the information about when cell types
  # become extinct somehow.

  selected_origins = select_timeseries_origins(possible_origins, initial_origin_policy)

  paths_to_origins = find_paths_to_origins(selected_origins, timeseries_pathfinding_graph)

  return(paths_to_origins)
}




