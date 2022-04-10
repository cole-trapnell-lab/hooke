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
select_origins <- function(ccm,
                           possible_origins,
                           selection_policy = c("acceptable-origins", "best-origin", "closest-origin", "all-origins"),
                           min_origin_pcor = 0.00,
                           fraction_origin_pcor_thresh=0.5,
                           allow_discontinuous_paths = FALSE,
                           max_origin_distance_ratio=2){

  selection_policy = match.arg(selection_policy)

  neg_rec_edges = possible_origins %>%
    select(neg_rec_edges) %>%
    tidyr::unnest(neg_rec_edges) %>%
    select(origin=from, destination=to, origin_pcor=pcor) %>% distinct() %>%
    mutate(origin_pcor = -origin_pcor)

  weighted_paths = possible_origins %>% select(path) %>%
    filter(!is.na(path)) %>%
    tidyr::unnest(path) %>%
    select(origin, destination, from, to, umap_dist=weight, edges_on_path=distance_from_root)


  #valid_origins = get_valid_origins(ccm, overlap_theshold)
  #weighted_paths = left_join(weighted_paths, valid_origins, by=c("destination"="id"))
  #weighted_paths = weighted_paths %>% filter(origin %in% unlist(possible_origins)) %>% select(-possible_origins) %>% distinct()
  weighted_paths = left_join(weighted_paths, neg_rec_edges)

  weighted_paths = weighted_paths %>% mutate(from = ifelse(is.na(from), origin, from),
                                             to = ifelse(is.na(to), destination, to))
  if (allow_discontinuous_paths == FALSE){
    weighted_paths = weighted_paths %>% filter(is.na(umap_dist) == FALSE)
  }

  weighted_paths = weighted_paths %>% filter(origin_pcor > min_origin_pcor)

  # origin_summary = weighted_paths %>% group_by(destination) %>% summarize(num_origins = length(unique(origin)),
  #                                                                         num_indirect_origins = length(unique(origin[umap_dist != -1])))

  weighted_paths = weighted_paths %>% distinct()

  weighted_paths = weighted_paths %>% group_by(origin, destination) %>% mutate(path_geodesic_umap_dist = sum(umap_dist)) %>% ungroup()

  weighted_paths = weighted_paths %>% group_by(destination) %>% mutate(max_origin_pcor = max(origin_pcor),
                                                                       fraction_max_origin_pcor = origin_pcor / max_origin_pcor)
  weighted_paths = weighted_paths %>% group_by(destination) %>% mutate(closest_origin_dist = min(path_geodesic_umap_dist),
                                                                       distance_ratio_over_closest = path_geodesic_umap_dist / closest_origin_dist)
  distance_to_best_origin_df = weighted_paths %>%
    select(origin, destination, origin_pcor, max_origin_pcor, path_geodesic_umap_dist) %>%
    slice_max(origin_pcor, with_ties=FALSE) %>%
    select(destination, distance_to_best_origin = path_geodesic_umap_dist)
  weighted_paths = left_join(weighted_paths, distance_to_best_origin_df)

  weighted_paths = weighted_paths %>% mutate (distance_ratio_over_best = path_geodesic_umap_dist / distance_to_best_origin)

  if (selection_policy == "acceptable-origins"){
    weighted_paths = weighted_paths %>% filter(distance_ratio_over_best < max_origin_distance_ratio)
    weighted_paths = weighted_paths %>% filter(fraction_max_origin_pcor > fraction_origin_pcor_thresh)
  }else if(selection_policy == "best-origin"){
    weighted_paths = weighted_paths %>% filter(fraction_max_origin_pcor == 1.0)
  }else if(selection_policy == "closest-origin"){
    weighted_paths = weighted_paths %>% filter(distance_ratio_over_closest == 1.0)
  }else if (selection_policy == "all"){

  }

  origin_df = weighted_paths %>%
    select(from, to, origin_pcor, fraction_max_origin_pcor) %>%
    group_by(from, to) %>%
    summarize(origin_pcor = sum(origin_pcor),
              fraction_max_origin_pcor = sum(fraction_max_origin_pcor),
              n=1) %>%
    group_by(to) %>%
    mutate(scaled_weight=fraction_max_origin_pcor / max(fraction_max_origin_pcor))
  #tally() %>%
  #mutate(scaled_weight = abs(n) / max(abs(n)))

  return(origin_df)
}

#' Identify the possible origins for each destination
#'
#' @export
find_origins <- function(ccm,
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

  if (is.null(start)){
    start = min(colData(ccm@ccs)[,interval_col])
  }
  if (is.null(stop)){
    stop = max(colData(ccm@ccs)[,interval_col])
  }

  timepoints = seq(start, stop, interval_step)

  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start,
                                              stop,
                                              interval_col=interval_col,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              ...)

  # FIXME: this is currently hardcoded to use "adjusted_time". Not generic and urgently needs fixing.
  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step, ...)

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = data.frame(id=cell_groups)

  paga_graph = initial_pcor_graph(ccm@ccs) %>% igraph::graph_from_data_frame(directed = FALSE, vertices=node_metadata) %>% igraph::as.directed()
  cov_graph = hooke:::return_igraph(model(ccm, "reduced"))
  cov_graph_edges = igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight != 0.00)

  # FIXME: While we should consider this more sophisticated weighting function, let's test just using UMAP distance for now:
  #weighted_edges = hooke:::get_weighted_edges(ccm, all_edges)
  weighted_edges = weigh_edges_by_umap_dist(ccm, cov_graph_edges)

  weighted_edges = weighted_edges %>% group_by(from, to) %>%
    mutate(adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to)) %>% ungroup()
  weighted_edges = weighted_edges %>% filter(adjacent_in_paga)

  traversal_graph = weighted_edges %>% select(from, to, weight) %>%
    igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
    igraph::as.directed()

  edges_between_concurrent_states = left_join(traversal_graph %>% igraph::as_data_frame(what="edges"),
                                              extant_cell_type_df, by=c("from"="cell_group")) %>%
    select(from, to, from_start=longest_contig_start, from_end=longest_contig_end)
  edges_between_concurrent_states = left_join(edges_between_concurrent_states, extant_cell_type_df, by=c("to"="cell_group")) %>%
    select(from, to, from_start, from_end, to_start=longest_contig_start, to_end=longest_contig_end)
  edges_between_concurrent_states = edges_between_concurrent_states %>%
    filter(from_start <= to_start & from_end >= to_start) %>% select(from, to) %>% distinct()

  if (require_presence_at_all_timepoints){
    traversal_graph = igraph::intersection(traversal_graph,
                                           edges_between_concurrent_states %>%
                                             igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata))
  }

  # there must be a cleaner way to specify this edge sequence, but whatev
  #traversal_graph = igraph::delete_edges(traversal_graph,
  #                                       edges_between_disjoint_states %>%
  #                                         mutate(edge_id=paste(from, to, sep="|")) %>% pull(edge_id))

  time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval)

  time_contrasts = find_origin_paths_from_time_contrasts(ccm,
                                                         traversal_graph,
                                                         time_contrasts,
                                                         timepoint_pred_df,
                                                         extant_cell_type_df,
                                                         interval_col,
                                                         q_val,
                                                         require_presence_at_all_timepoints)
  #possible_origins = select_origins(ccm, times)

  closest_origins_paths = select_origins(ccm, time_contrasts, selection_policy = initial_origin_policy)

  # May not actually need this bit below:

  origin_edge_graph = closest_origins_paths %>% select(from, to, origin_pcor) %>% igraph::graph_from_data_frame()
  A = origin_edge_graph %>% igraph::as_adjacency_matrix(attr="origin_pcor")
  B = origin_edge_graph %>% igraph::as_adjacency_matrix(attr="origin_pcor")

  A[lower.tri(A)] = 0
  B[upper.tri(B)] = 0

  # Remove the edges we want to keep from A and B below.
  # Note that this will *remove* bidirectional edges from the filtered path set if they have equally good pcor,
  # and will therefore NOT be removed later from the traversal graph (so paths will still be able to cross them)
  A[A >= Matrix::t(B)] = 0
  B[B >= Matrix::t(A)] = 0

  edges_to_drop_mat = A + B
  edges_to_drop_df = edges_to_drop_mat %>% igraph::graph_from_adjacency_matrix(weighted="origin_pcor") %>% igraph::as_data_frame()
  closest_origins_paths_filtered = inner_join(edges_to_drop_df, closest_origins_paths) %>% as_tibble()

  edges_to_remove_from_travesal_graph = closest_origins_paths_filtered %>% select(from, to)

  # there must be a cleaner way to specify this edge sequence, but whatev
  traversal_graph = igraph::delete_edges(traversal_graph,
                                         edges_to_remove_from_travesal_graph %>%
                                           mutate(edge_id=paste(from, to, sep="|")) %>% pull(edge_id))

  time_contrasts = find_origin_paths_from_time_contrasts(ccm,
                                                         traversal_graph,
                                                         time_contrasts,
                                                         timepoint_pred_df,
                                                         extant_cell_type_df,
                                                         interval_col,
                                                         q_val,
                                                         require_presence_at_all_timepoints)

  return (time_contrasts)
}
#debug(find_origins)





