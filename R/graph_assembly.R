#' Orient graph edges from a PLNnetwork using a contrast between conditions
#'
#' @param ccm A cell_count_model
#' @param cond_b_vs_a_tbl A contrast between two conditions as returned by compare_abundances()
#' @noRd
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

#' @noRd
return_igraph <- function(model, type = "partial_cor", remove.isolated=FALSE,
                          edge.color = c("#F8766D", "#00BFC4"),
                          output = "igraph"){

  net <- model$latent_network(type = type)

  if (output == "igraph") {

    G <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (remove.isolated) {
      G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
    }
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
                                  pct_dynamic_range = 0.25,
                                  pct_range_detection_thresh = pct_dynamic_range,
                                  min_cell_range = 2,
                                  ...){

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step=interval_step, ...)

  norm_mat = normalized_counts(ccm@ccs, "size_only")
  norm_mat[norm_mat == 0] = NA
  count_quantiles = sparseMatrixStats::rowQuantiles(norm_mat, probs = seq(from = 0, to = 1, by = pct_range_detection_thresh), na.rm=T)
  count_ranges = sparseMatrixStats::rowRanges(norm_mat, na.rm=T)
  row.names(count_ranges) = row.names(count_quantiles)

  #abund_range = range(timepoint_pred_df$log_abund)
  #dynamic_range = abund_range[2]-abund_range[1]

  cell_type_thresh_df = tibble(cell_group = timepoint_pred_df %>% pull(cell_group) %>% unique,
                               cell_group_pct_range_detection_thresh = count_quantiles[cell_group, 2],
                               min_count = count_ranges[cell_group, 1],
                               max_count = count_ranges[cell_group, 2])

  timepoint_pred_df = timepoint_pred_df %>% left_join(cell_type_thresh_df)

  if (is.null(log_abund_detection_thresh)) {
    log_abund_detection_thresh = abund_range[1] + pct_dynamic_range*dynamic_range
  }

  timepoint_pred_df = timepoint_pred_df %>%
    group_by(cell_group) %>%
    mutate(max_abundance = max(exp(log_abund)),
           percent_max_abund = exp(log_abund) / max_abundance,
           cell_type_prediction_range = max(log_abund)-(min(log_abund)),
           percent_cell_type_range = (log_abund - min(log_abund)) / cell_type_prediction_range,
           above_log_abund_thresh = (log_abund - 2*log_abund_se > log_abund_detection_thresh & log_abund - 2*log_abund_se > log(cell_group_pct_range_detection_thresh)) | cell_type_prediction_range < min_cell_range,
           present_flag = ifelse(above_log_abund_thresh, TRUE, NA)) %>%
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
           percent_cell_type_range,
           longest_contig_start,
           longest_contig_end,
           present_above_thresh)
  return(extant_cell_type_df)
}
#undebug(get_extant_cell_types)
#get_extant_cell_types(wt_ccm_wl, 72, 96) %>% filter(cell_group == "21")


#' @noRd
weigh_edges_by_umap_dist <- function(ccm, edges) {

  # get weighted path
  dist_df = hooke:::get_distances(ccm@ccs, matrix = F)

  weighted_edge_df = left_join(edges, dist_df, by=c("from", "to")) %>%
    mutate(weight = dist)


  return(weighted_edge_df)

}

#' @noRd
weigh_edges <- function(ccm, edges) {
  # get weighted path
  dist_df = hooke:::get_distances(ccm@ccs, matrix = F)


}

# FIXME: we need to standardize notation (use either "partition" or "component") throughout the code
#' @noRd
add_cross_component_pathfinding_links = function(ccm,
                                                 pathfinding_graph,
                                                 extant_cell_type_df,
                                                 type=c("strongest-pcor", "strong-pcor", "ctp"),
                                                 surprise_thresh=2,
                                                 components="partition"){

  pcor_graph = tibble(cell_group = row.names(counts(ccm@ccs)))
  pcor_graph = pcor_graph %>% select(cell_group) %>% tidyr::expand(cell_group, cell_group)
  colnames(pcor_graph) = c("from", "to")
  pcor_graph = pcor_graph %>% filter(from != to)

  nz_pcors = igraph::graph_from_adjacency_matrix(model(ccm, "reduced")$latent_network(type = "partial_cor"), mode = "undirected", weighted = TRUE, diag = FALSE) %>% igraph::as.directed()
  nz_pcor_graph_edges = igraph::as_data_frame(nz_pcors, what="edges") %>% dplyr::rename(pcor=weight) %>% select(from, to, pcor)#%>% dplyr::filter(pcor != 0.00)

  pcor_graph = dplyr::left_join(pcor_graph, nz_pcor_graph_edges)

  if (components == "partition"){
    component_assignments = partitions(ccm@ccs@cds)
  } else {
    component_assignments = colData(ccm@ccs@cds)[,components]
    names(component_assignments) = row.names(colData(ccm@ccs@cds))
  }

  clusters_by_partition = tibble(cell_group=ccm@ccs@metadata$cell_group_assignments$cell_group,
                                 partition=component_assignments[row.names(ccm@ccs@metadata$cell_group_assignments)])
  clusters_by_partition = clusters_by_partition %>% group_by(cell_group, partition) %>% summarize(cells_from_group_in_partition=n())
  clusters_by_partition = clusters_by_partition %>% group_by(cell_group) %>% slice_max(cells_from_group_in_partition, n=1)
  clusters_by_partition = clusters_by_partition %>% select(cell_group, partition)

  #FIXME: we should be indexing extant_cell_type_df by name, but this function doesn't know what the
  # timepoint column is called yet
  start_time = min(extant_cell_type_df[,1])
  clusters_by_partition = clusters_by_partition %>% left_join(extant_cell_type_df, by="cell_group")
  partitions_present_at_start = clusters_by_partition %>%
    filter(timepoint == start_time & present_above_thresh) %>%
    pull(partition) %>% unique
  partitions_absent_at_start = setdiff(clusters_by_partition$partition, partitions_present_at_start)
  cell_groups_in_emergent_partitions = clusters_by_partition %>%
    filter(partition %in% partitions_absent_at_start) %>% pull(cell_group) %>% unique

  # only add cross-partition links between cell groups that aren't present
  # at the start of the assembly time window
  cross_partition_map = clusters_by_partition %>% ungroup() %>% select(cell_group) %>% tidyr::expand(cell_group, cell_group)
  colnames(cross_partition_map) = c("from", "to")
  cross_partition_map = cross_partition_map %>% filter(from != to)
  # Only allow cross partition links when one cell group is in an emergent partition
  cross_partition_map = cross_partition_map %>% filter(from %in% cell_groups_in_emergent_partitions |
                                                       to %in% cell_groups_in_emergent_partitions)
  cross_partition_map = dplyr::left_join(cross_partition_map, clusters_by_partition %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
  cross_partition_map = dplyr::left_join(cross_partition_map, clusters_by_partition %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group")) #%>%
  cross_partition_map = cross_partition_map %>% filter (from_partition != to_partition)

  # only add cross-partition links between cell groups that are present
  # at the same time
  cross_partition_map = cross_partition_map %>% filter (from_timepoint == to_timepoint &
                                                          from_present_above_thresh &
                                                          to_present_above_thresh)

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
    mutate(model_fit_lower = model_fit - surprise_thresh  * model_sd,
           #model_fit_lower = ifelse(model_fit_lower > 0, model_fit_lower, min(abs_pcor)),
           model_fit_upper = model_fit + surprise_thresh * model_sd)

  cross_partition_map = cross_partition_map %>% left_join(pcor_graph, by=c("from"="from", "to"="to"))



  if (type == "strongest-pcor"){
    cross_partition_edges =  cross_partition_map %>%
      group_by(to_partition) %>%
      group_by(from_partition, to_partition) %>%
      slice_max(abs_pcor, n=1) %>% ungroup() %>%
      #filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
      select(from, to, weight=umap_dist)
  }else if(type == "strong-pcor"){
    cross_partition_edges =  cross_partition_map %>%
      filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
      select(from, to, weight=umap_dist)
  }else if (type == "ctp"){
    cross_partition_edges =  cross_partition_map %>%
      group_by(to_partition) %>%
      slice_min(to_timepoint) %>%
      slice_min(umap_dist) %>%
      ungroup() %>%
      select(from, to, weight=umap_dist)
  }

  cross_partition_graph = igraph::graph_from_data_frame(cross_partition_edges, directed=FALSE, vertices=data.frame(id=row.names(counts(ccm@ccs)))) %>% igraph::as.directed()


  #unexpectedly_strong_pcor_edges = pcor_graph %>%
  #  filter (pcor < model_fit_lower | pcor > model_fit_upper) %>%
  #  select(from, to, weight=umap_dist)
  #uspe_graph = igraph::graph_from_data_frame(unexpectedly_strong_pcor_edges, directed=TRUE, vertices=data.frame(id=row.names(counts(ccm@ccs))))
  updated_pathfinding_graph = igraph::union(pathfinding_graph, cross_partition_graph)
  updated_pathfinding_graph = igraph::simplify(updated_pathfinding_graph)
  return(updated_pathfinding_graph)
}

#' Initialize a graph over which cells can transition
#'
#' @param ccm A cell count model
#' @param extant_cell_type_df A data frame describing when cell types are present, generated by get_extant_cell_types()
#' @param allow_links_between_components Whether cells in separate partitions of the UMAP can be linked by paths
#' @export
init_pathfinding_graph <- function(ccm,
                                   extant_cell_type_df,
                                   links_between_components=c("ctp", "none", "strongest-pcor", "strong-pcor"),
                                   components="partition",
                                   weigh_by_pcor=F,
                                   edge_whitelist=NULL,
                                   edge_blacklist=NULL){

  # There are a number of different ways we could set up this "pathfinding graph" but for now
  # Let's just use the PAGA (weighed by distance in UMAP space), subtracting edges between which
  # there is zero partial correlation in the nuisance cell count model.

  links_between_components = match.arg(links_between_components)

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = data.frame(id=cell_groups)

  if (is.null(edge_whitelist)){
    message("Initializing pathfinding graph from partially correlated pairs linked in PAGA")
    paga_graph = initial_pcor_graph(ccm@ccs) %>% igraph::graph_from_data_frame(directed = FALSE, vertices=node_metadata) %>% igraph::as.directed()
    cov_graph = hooke:::return_igraph(model(ccm, "reduced"))
    cov_graph_edges = igraph::as_data_frame(cov_graph, what="edges")

    if (nrow(cov_graph_edges) > 0 & "weight" %in% colnames(cov_graph_edges)){
      cov_graph_edges = cov_graph_edges %>%
        dplyr::rename(pcor=weight) %>%
        dplyr::filter(pcor != 0.00)
    }

    cov_graph_edges$to = as.character(cov_graph_edges$to)
    cov_graph_edges$from = as.character(cov_graph_edges$from)

    weighted_edges = hooke:::weigh_edges_by_umap_dist(ccm, cov_graph_edges)

    paga_components = igraph::components(paga_graph)
    same_partition_mat = outer(paga_components$membership, paga_components$membership, FUN="==")
    weighted_edges = weighted_edges %>% group_by(from, to) %>%
      mutate(adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to),
             same_partition = same_partition_mat[from, to]) %>% ungroup()

    weighted_edges = weighted_edges %>% filter(adjacent_in_paga)
    pathfinding_graph = weighted_edges %>% select(from, to, weight) %>%
      igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
      igraph::as.directed()

  }else{
    weighted_edges = hooke:::weigh_edges_by_umap_dist(ccm, edge_whitelist)
    pathfinding_graph = weighted_edges %>% select(from, to, weight) %>%
      igraph::graph_from_data_frame(directed=TRUE, vertices=node_metadata)
  }



  if (links_between_components != "none"){
    pathfinding_graph = add_cross_component_pathfinding_links(ccm,
                                                              pathfinding_graph,
                                                              extant_cell_type_df,
                                                              type=links_between_components,
                                                              components=components)
    pathfinding_graph = igraph::as_data_frame(pathfinding_graph) %>%  select(from, to)
    pathfinding_graph = hooke:::weigh_edges_by_umap_dist(ccm, pathfinding_graph) %>%
      igraph::graph_from_data_frame(directed=FALSE, vertices=node_metadata) %>%
      igraph::as.directed() %>% igraph::simplify()
  }

  if (is.null(edge_blacklist) == FALSE){
    blacklist_graph = igraph::graph_from_data_frame(edge_blacklist, directed=TRUE, vertices=node_metadata)
    pathfinding_graph = pathfinding_graph - blacklist_graph
  }

  # if (is.null(edge_whitelist) == FALSE){
  #   # Ensuring whitelisted edges remain in pathfinding graph
  #   edge_whitelist = edge_whitelist[,c(1,2)] %>% as_tibble()
  #   colnames(edge_whitelist) = c("from", "to")
  #   pathfinding_graph = igraph::as_data_frame(pathfinding_graph) %>%  select(from, to)
  #   pathfinding_graph = pathfinding_graph %>% bind_rows(edge_whitelist) %>% distinct()
  #   pathfinding_graph = hooke:::weigh_edges_by_umap_dist(ccm, pathfinding_graph)
  # }

  return (pathfinding_graph)
}


#' Generate a blacklist of state transition relationships based on a perturbation
#' experiment.
#'
#' This function generates a blacklist of edges between cell states based on
#' the idea that if cell state X is lost, and Y is not ever lost in the experiment
#' Y can't come from X. This function is very simplistic right now. We could get
#' more sophisticated by looking at the relative timing of losses, etc. We are
#' also not accounting for power at all. We should be asking whether we have
#' power to detect a change in state Y, and if not, exclude (X,Y) from the
#' blacklist
#'
#' @export
get_discordant_loss_pairs <- function(perturbation_ccm,
                                      perturb_time_window,
                                      control_timeseries_ccm,
                                      control_time_window,
                                      interval_step,
                                      interval_col,
                                      log_abund_detection_thresh,
                                      q_val,
                                      model_for_pcors="reduced",
                                      min_pathfinding_lfc=0,
                                      ...){

  #print ("getting perturbation paths")
  #print (time_window)
  if (is.null(perturbation_ccm) || is.na(perturbation_ccm)){
    message("Error: no perturbation model")
    return(NA)
  }

  perturb_start_time = min(as.numeric(perturb_time_window$start_time))
  perturb_stop_time = min(as.numeric(perturb_time_window$stop_time))

  control_start_time = min(as.numeric(control_time_window$start_time))
  control_stop_time = min(as.numeric(control_time_window$stop_time))

  # If the perturbation experiment model covers a winder time interval than the
  # control timeseries, just limit testing to the perturbation experiment's
  # time interval
  if (perturb_start_time < control_start_time)
    perturb_start_time = control_start_time

  if (perturb_stop_time < control_stop_time)
    perturb_stop_time = control_stop_time

  wt_timepoint_pred_df = estimate_abundances_over_interval(control_timeseries_ccm, control_start_time, control_stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
  peak_wt_abundance = wt_timepoint_pred_df %>% group_by(cell_group) %>% slice_max(log_abund, n=1)
  peak_outside_perturbation_window = peak_wt_abundance %>%
    filter(!!sym(interval_col) > perturb_stop_time | !!sym(interval_col) < perturb_start_time) %>%
    pull(cell_group)

  message ("\tEstimating loss timing")
  earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
                                                   start_time=perturb_start_time,
                                                   stop_time=perturb_stop_time,
                                                   interval_step = interval_step,
                                                   interval_col=interval_col,
                                                   control_ccm=control_timeseries_ccm,
                                                   control_start_time=control_start_time,
                                                   control_stop_time=control_stop_time,
                                                   log_abund_detection_thresh=log_abund_detection_thresh,
                                                   q_val = q_val,
                                                   delta_log_abund_loss_thresh=min_pathfinding_lfc,
                                                   ...)

  lost_cell_groups = earliest_loss_tbl %>% filter (is_lost_at_peak) %>% pull(cell_group) %>% unique

  cell_groups_that_peak_within_perturb_window = earliest_loss_tbl %>%
    filter (peak_time_in_ctrl_within_perturb_time_range) %>%
    pull(cell_group) %>% unique

  unaffected_cell_groups = setdiff(cell_groups_that_peak_within_perturb_window, lost_cell_groups)

  # Exclude states that peak outside the window of this perturbation experiment,
  # as we may simply have not yet seen their loss.
  unaffected_cell_groups = setdiff(unaffected_cell_groups, peak_outside_perturbation_window)

  discordant_loss_pairs = tidyr::expand_grid(lost_cell_groups, unaffected_cell_groups)

  return (discordant_loss_pairs)
}

#' @noRd
# get_perturbation_paths <- function(perturbation_ccm,
#                                    time_window,
#                                    pathfinding_graph,
#                                    interval_step,
#                                    interval_col,
#                                    log_abund_detection_thresh,
#                                    q_val,
#                                    control_ccm=perturbation_ccm,
#                                    control_time_window=time_window,
#                                    model_for_pcors="reduced",
#                                    min_pathfinding_lfc=0,
#                                    batch_col=batch_col,
#                                    ...){
#
#   #print ("getting perturbation paths")
#   #print (time_window)
#
#   start_time = min(as.numeric(time_window$start_time))
#   stop_time = min(as.numeric(time_window$stop_time))
#
#   control_start_time = min(as.numeric(control_time_window$start_time))
#   control_stop_time = min(as.numeric(control_time_window$stop_time))
#
#   #print (start_time)
#   #print (stop_time)
#
#   message ("\tEstimating loss timing")
#   earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
#                                                    start_time=start_time,
#                                                    stop_time=stop_time,
#                                                    interval_step = interval_step,
#                                                    interval_col=interval_col,
#                                                    control_ccm=control_ccm,
#                                                    control_start_time=control_start_time,
#                                                    control_stop_time=control_stop_time,
#                                                    log_abund_detection_thresh=log_abund_detection_thresh,
#                                                    q_val = q_val,
#                                                    delta_log_abund_loss_thresh=min_pathfinding_lfc,
#                                                    ...)
#   #print (earliest_loss_tbl)
#
#   # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
#   old_omp_num_threads = single_thread_omp()
#   old_blas_num_threads = single_thread_blas()
#
#   tryCatch({
#     pln_model = model(perturbation_ccm, model_for_pcors)
#
#     change_corr_tbl = earliest_loss_tbl %>% filter (is.finite(earliest_time)) %>% select(cell_group)
#     change_corr_tbl = change_corr_tbl %>% select(cell_group) %>% tidyr::expand(cell_group, cell_group)
#     colnames(change_corr_tbl) = c("from", "to")
#     change_corr_tbl = change_corr_tbl %>% filter(from != to)
#
#     lost_cell_groups = change_corr_tbl %>% pull(from) %>% unique
#
#     #print ("change corr tbl")
#     #print (change_corr_tbl)
#
#     #corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
#     change_corr_tbl = dplyr::left_join(change_corr_tbl, earliest_loss_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
#     #dplyr::rename(log_abund_x,
#     #              to_delta_log_abund = delta_log_abund)
#     change_corr_tbl = dplyr::left_join(change_corr_tbl, earliest_loss_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%
#
#     change_corr_tbl = change_corr_tbl %>% mutate_if(is.numeric, tidyr::replace_na, replace = -Inf)
#
#     #print (change_corr_tbl)
#
#     #print ("possible loss pairs")
#     #print (change_corr_tbl)
#
#     #print ("debug pairs")
#     #change_corr_tbl %>% filter(from == "24") %>% print
#     concordant_fwd_loss_pairs = change_corr_tbl %>%
#       #dplyr::filter(pcor < 0) %>% # do we just want negative again?
#       dplyr::filter(is.finite(from_is_lost_when_present) & is.finite(to_is_lost_when_present) & from_peak_wt_time <= to_peak_wt_time)
#       #dplyr::filter(is.finite(from_peak_loss_time) & is.finite(to_peak_loss_time) & from_peak_loss_time <= to_peak_loss_time)
#
#       #dplyr::filter((is.finite(from_peak_loss_time) & is.finite(to_peak_loss_time) & from_peak_loss_time <= to_peak_loss_time) |
#       #                (is.finite(from_largest_loss_time) & is.finite(to_largest_loss_time) & from_largest_loss_time <= to_largest_loss_time) |
#       #                (is.finite(from_earliest_time) & is.finite(to_earliest_time) & from_earliest_time < to_earliest_time))
#
#     message ("\tEstimating loss timing  (again)")
#     earliest_loss_tbl_no_sig_filter = hooke:::estimate_loss_timing(perturbation_ccm,
#                                                                    start_time=start_time,
#                                                                    stop_time=stop_time,
#                                                                    interval_step = interval_step,
#                                                                    control_ccm=control_ccm,
#                                                                    control_start_time=control_start_time,
#                                                                    control_stop_time=control_stop_time,
#                                                                    interval_col=interval_col,
#                                                                    log_abund_detection_thresh=log_abund_detection_thresh,
#                                                                    q_val = 1,
#                                                                    delta_log_abund_loss_thresh=min_pathfinding_lfc,
#                                                                    ...)
#
#     lost_cell_groups
#     lost_cell_group_pathfinding_subgraph = igraph::subgraph(pathfinding_graph,
#                                                             earliest_loss_tbl %>% filter (is_lost_when_present) %>% pull(cell_group))
#
#     #lost_cell_group_pathfinding_subgraph = pathfinding_graph
#
#     # Compute the shortest paths between each of those node pairs in the pathfinding graph
#     concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>% select(from, to) %>% distinct()
#
#     message (paste("\tfinding shortest paths between", nrow(concordant_fwd_loss_pairs), "loss pairs"))
#     #print (concordant_fwd_loss_pairs)
#
#     paths_between_concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>%
#       mutate(path = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
#                                        .x = from, .y = to,
#                                        lost_cell_group_pathfinding_subgraph,
#                                        .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
#                                        .progress=TRUE))
#
#     paths_between_concordant_fwd_loss_pairs = paths_between_concordant_fwd_loss_pairs %>%
#       filter(is.na(path) == FALSE)
#
#     #print ("loss pair paths")
#     #print (paths_between_concordant_fwd_loss_pairs)
#
#     #print ("debug pair paths")
#     #paths_between_concordant_fwd_loss_pairs %>% filter(is.na(path) == FALSE) %>% rename(origin=from, destination=to) %>% tidyr::unnest(path) %>% print(n=1000)
#     #print (paths_between_concordant_fwd_loss_pairs)#  %>% tidyr::unnest(path) %>% filter(from == "24") %>% print
#
#     message ("\tscoring paths between loss pairs")
#     paths_between_concordant_fwd_loss_pairs = score_paths_for_perturbations(perturbation_ccm,
#                                                                             paths_between_concordant_fwd_loss_pairs,
#                                                                             perturbation_col="knockout",
#                                                                             interval_col=interval_col,
#                                                                             batch_col=batch_col)
#     return (paths_between_concordant_fwd_loss_pairs)
#   }, error = function(e) {
#     print (e)
#     return (NA)
#   }, finally = {
#     #RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
#     #RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
#   })
#
#
# }


get_perturbation_paths <- function(perturbation_ccm,
                                   perturb_summary_tbl,
                                   pathfinding_graph){

  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  old_omp_num_threads = single_thread_omp()
  old_blas_num_threads = single_thread_blas()

  tryCatch({

    lost_cell_groups = perturb_summary_tbl %>% filter(is_lost_when_present) %>% pull(cell_group) %>% unique
    gained_cell_groups = perturb_summary_tbl %>% filter(is_gained_when_present) %>% pull(cell_group) %>% unique

    perturbed_cell_groups = union(lost_cell_groups, gained_cell_groups)
    lost_subgraph = igraph::subgraph(pathfinding_graph, lost_cell_groups)
    nodes_without_lost_parent = igraph::V(lost_subgraph)[igraph::degree(lost_subgraph, mode="in") == 0]$name

    # We will do pathfinding to nodes that are lost, starting from nodes that
    # are either lost and have no parents that are also lost, or nodes that are
    # gained.
    if (length(lost_cell_groups) == 0)
      return (NA)
    loss_neighborhood_graph = igraph::subgraph(pathfinding_graph, lost_cell_groups)
    perturb_pathfinding_graph = loss_neighborhood_graph

    if (length(gained_cell_groups) > 0){
      gain_neighborhood_graph = do.call(igraph::union,igraph::make_ego_graph(pathfinding_graph, order = 1,
                                                                             nodes = gained_cell_groups,
                                                                             mode = "out"))
      perturb_pathfinding_graph = igraph::union(perturb_pathfinding_graph, gain_neighborhood_graph)
    }

    perturb_pathfinding_graph = igraph::subgraph(pathfinding_graph, igraph::V(perturb_pathfinding_graph)$name)

    perturb_pathfinding_graph = igraph::as_data_frame(perturb_pathfinding_graph) %>% select(from, to) %>% distinct()

    perturb_pathfinding_graph = left_join(perturb_pathfinding_graph,
                                          igraph::as_data_frame(pathfinding_graph),
                                        by=c("from", "to"))
    perturb_pathfinding_graph = igraph::graph_from_data_frame(perturb_pathfinding_graph)

    perturb_pathfinding_graph_nodes = igraph::V(perturb_pathfinding_graph)$name

    possible_origins = union(nodes_without_lost_parent, gained_cell_groups)

    #origin_dest_pairs_to_test = tibble(cell_group=loss_neighborhood_graph_nodes) %>%
    #  tidyr::expand(cell_group, cell_group)
    origin_dest_pairs_to_test =  expand.grid(possible_origins, lost_cell_groups) %>% as_tibble()
    colnames(origin_dest_pairs_to_test) = c("from", "to")
    origin_dest_pairs_to_test$from = as.character(origin_dest_pairs_to_test$from)
    origin_dest_pairs_to_test$to = as.character(origin_dest_pairs_to_test$to)
    origin_dest_pairs_to_test = origin_dest_pairs_to_test %>% filter(from != to)

    origin_dest_pairs_to_test = dplyr::left_join(origin_dest_pairs_to_test, perturb_summary_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
    origin_dest_pairs_to_test = dplyr::left_join(origin_dest_pairs_to_test, perturb_summary_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%

    concordant_fwd_loss_pairs = origin_dest_pairs_to_test %>%
      #dplyr::filter(pcor < 0) %>% # do we just want negative again?
      dplyr::filter(is.na(to_is_lost_when_present) == FALSE & to_is_lost_when_present & from_peak_wt_time <= to_peak_wt_time)

    # Compute the shortest paths between each of those node pairs in the pathfinding graph
    concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>% select(from, to) %>% distinct()
    message (paste("\tfinding shortest paths between", nrow(concordant_fwd_loss_pairs), "possible loss pairs"))

    paths_between_concordant_fwd_loss_pairs = concordant_fwd_loss_pairs %>%
      mutate(path = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                                       .x = from, .y = to,
                                       perturb_pathfinding_graph,
                                       .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                       .progress=TRUE))

    paths_between_concordant_fwd_loss_pairs = paths_between_concordant_fwd_loss_pairs %>%
      filter(is.na(path) == FALSE)

    if (nrow(paths_between_concordant_fwd_loss_pairs) == 0)
      return(NA)
    message (paste("\tscoring", nrow(paths_between_concordant_fwd_loss_pairs),"paths between loss pairs"))
    paths_between_concordant_fwd_loss_pairs = score_paths_for_perturbations(perturbation_ccm,
                                                                            paths_between_concordant_fwd_loss_pairs,
                                                                            loss_tbl=perturb_summary_tbl)
    return (paths_between_concordant_fwd_loss_pairs)
  }, error = function(e) {
    print (e)
    return (NA)
  }, finally = {
    #RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    #RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })


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
#' @noRd
measure_time_delta_along_path <- function(path_df, ccs, cells_along_path_df, interval_col="timepoint") {


  #cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df = cells_along_path_df %>% filter(cell_group %in% vertices)

  # cells_along_path_df = normalized_counts(ccs, "size_only", pseudocount = 0)[vertices,] %>%
  #   as.matrix() %>%
  #   Matrix::t() %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column() %>%
  #   tidyr::pivot_longer(!matches("rowname")) %>%
  #   rename(sample=rowname, cell_group=name, num_cells=value)

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
#' @noRd
# measure_perturbation_freq_along_path <- function(path_df, ccs, cells_along_path_df, perturbation_col="knockout", interval_col="timepoint", batch_col=NULL) {
#
#
#   #print ("measuring perturbation score")
#   #cds = ccs@cds
#   vertices = union(path_df$to, path_df$from) %>% unique()
#   path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))
#
#   cells_along_path_df = cells_along_path_df %>% filter(cell_group %in% vertices)
#
#   num_batches = 1
#   if (is.null(batch_col)){
#     cells_along_path_df = cells_along_path_df %>%
#       left_join(colData(ccs) %>% as.data.frame %>%
#                   select(sample, !!sym(perturbation_col), !!(interval_col)) %>%
#                   as_tibble(), by = c("sample" = "sample"))
#   }else{
#     cells_along_path_df = cells_along_path_df %>%
#       left_join(colData(ccs) %>% as.data.frame %>%
#                   select(sample, !!sym(perturbation_col), !!(interval_col), !!(batch_col)) %>%
#                   as_tibble(), by = c("sample" = "sample"))
#     #print(head(cells_along_path_df))
#     num_batches = cells_along_path_df %>% pull(!!sym(batch_col)) %>% unique %>% length
#   }
#
#   cells_along_path_df = cells_along_path_df %>%
#     left_join(path_df %>% select(-from), by = c("cell_group" = "to"))
#
#   #print (cells_along_path_df)
#   cells_along_path_df = cells_along_path_df %>% mutate_if(is.numeric, tidyr::replace_na, replace = 0)
#   cells_along_path_df = cells_along_path_df %>% dplyr::filter (num_cells > 0)
#   #cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
#   #cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)
#
#   cells_along_path_df$perturb = cells_along_path_df[[perturbation_col]]
#   cells_along_path_df$timepoint = cells_along_path_df[[interval_col]]
#
#   if (is.null(batch_col) || num_batches == 1){
#     path_model = glm(#perturb ~ timepoint + geodesic_dist,
#       perturb ~ geodesic_dist,
#       data=cells_along_path_df,
#       weights=cells_along_path_df$num_cells,
#       family=binomial(),
#       singular.ok = FALSE)
#   }else{
#     cells_along_path_df$batch = cells_along_path_df[[batch_col]]
#     #print (head(cells_along_path_df))
#     path_model = glm(#perturb ~ timepoint + geodesic_dist,
#       perturb ~ batch + geodesic_dist,
#       #perturb ~ geodesic_dist,
#       data=cells_along_path_df,
#       weights=cells_along_path_df$num_cells,
#       family=binomial(),
#       singular.ok = FALSE)
#     #print (summary(path_model))
#   }
#
#   path_model_tidied = broom::tidy(path_model)
#   path_model_glanced = broom::glance(path_model)
#   path_model_glanced$pseudo.r.squared = 1 - (path_model_glanced$deviance / path_model_glanced$null.deviance)
#
#   dist_param_row = nrow(path_model_tidied)
#   pm_stats = tibble(perturb_dist_effect = unlist(path_model_tidied[dist_param_row, "estimate"]),
#                     perturb_dist_effect_pval = unlist(path_model_tidied[dist_param_row, "p.value"]),
#                     perturb_dist_model_adj_rsq = unlist(path_model_glanced[1, "pseudo.r.squared"]),
#                     perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]))
#   # if (unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) { # should not happen
#   #   print (cells_along_path_df %>% as.data.frame)
#   #   print (summary(path_model))
#   # }
#
#   if (path_model$converged == FALSE | unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0)
#   {
#     return(tibble(perturb_dist_effect = 0,
#                   perturb_dist_effect_pval = 1,
#                   perturb_dist_model_adj_rsq = 0,
#                   perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"])))
#   }
#   return(pm_stats)
#   #return(coef(path_model)[["distance_from_root"]])
# }
# #debug(cells_along_path)

measure_perturbation_freq_along_path <- function(path_df, ccs, cells_along_path_df, perturbation_col="knockout", interval_col="timepoint", batch_col=NULL) {


  #print ("measuring perturbation score")
  #cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df = cells_along_path_df %>% filter(cell_group %in% vertices)

  num_batches = 1
  if (is.null(batch_col)){
    cells_along_path_df = cells_along_path_df %>%
      left_join(colData(ccs) %>% as.data.frame %>%
                  select(sample, !!sym(perturbation_col), !!(interval_col)) %>%
                  as_tibble(), by = c("sample" = "sample"))
  }else{
    cells_along_path_df = cells_along_path_df %>%
      left_join(colData(ccs) %>% as.data.frame %>%
                  select(sample, !!sym(perturbation_col), !!(interval_col), !!(batch_col)) %>%
                  as_tibble(), by = c("sample" = "sample"))
    #print(head(cells_along_path_df))
    num_batches = cells_along_path_df %>% pull(!!sym(batch_col)) %>% unique %>% length
  }

  cells_along_path_df = cells_along_path_df %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  #print (cells_along_path_df)
  cells_along_path_df = cells_along_path_df %>% mutate_if(is.numeric, tidyr::replace_na, replace = 0)
  cells_along_path_df = cells_along_path_df %>% dplyr::filter (num_cells > 0)
  #cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  #cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$perturb = cells_along_path_df[[perturbation_col]]
  cells_along_path_df$timepoint = cells_along_path_df[[interval_col]]

  if (is.null(batch_col) || num_batches == 1){
    path_model = glm(#perturb ~ timepoint + geodesic_dist,
      num_cells ~ cell_group + perturb,
      data=cells_along_path_df,
      #weights=cells_along_path_df$num_cells,
      family=quasipoisson(),
      singular.ok = FALSE)
  }else{
    cells_along_path_df$batch = cells_along_path_df[[batch_col]]
    #print (head(cells_along_path_df))
    path_model = glm(#perturb ~ timepoint + geodesic_dist,
      num_cells ~ batch + cell_group + perturb,
      #perturb ~ geodesic_dist,
      data=cells_along_path_df,
      #weights=cells_along_path_df$num_cells,
      family=quasipoisson(),
      singular.ok = FALSE)
    #print (summary(path_model))
  }

  path_model_tidied = broom::tidy(path_model)
  path_model_glanced = broom::glance(path_model)
  path_model_glanced$pseudo.r.squared = 1 - (path_model_glanced$deviance / path_model_glanced$null.deviance)

  dist_param_row = nrow(path_model_tidied)
  pm_stats = tibble(perturb_dist_effect = unlist(path_model_tidied[dist_param_row, "estimate"]),
                    perturb_dist_effect_pval = unlist(path_model_tidied[dist_param_row, "p.value"]),
                    perturb_dist_model_adj_rsq = unlist(path_model_glanced[1, "pseudo.r.squared"]),
                    perturb_dist_model_ncells = unlist(path_model_glanced[1, "nobs"]))
  # if (unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0) { # should not happen
  #   print (cells_along_path_df %>% as.data.frame)
  #   print (summary(path_model))
  # }

  if (path_model$converged == FALSE | unlist(path_model_glanced[1, "pseudo.r.squared"]) < 0)
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
#' @noRd
build_timeseries_transition_graph <- function(ccm,
                                              extant_cell_type_df,
                                              pathfinding_graph,
                                              q_val=0.01,
                                              start_time = NULL,
                                              stop_time = NULL,
                                              interval_col="timepoint",
                                              interval_step = 2,
                                              min_interval = 4,
                                              max_interval = 24,
                                              min_pathfinding_lfc=0,
                                              make_dag=FALSE,
                                              ...) {
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  #old_omp_num_threads = single_thread_omp()
  #old_blas_num_threads = single_thread_blas()

  # First, let's figure out when each cell type is present and
  # which ones emerge over the course of the caller's time interval
  if (is.null(start_time)){
    start_time = min(colData(ccm@ccs)[,interval_col])
  }
  if (is.null(stop_time)){
    stop_time = max(colData(ccm@ccs)[,interval_col])
  }

  timepoints = seq(start_time, stop_time, interval_step)

  message("Estimating abundances over time interval")
  timepoint_pred_df = estimate_abundances_over_interval(ccm, start_time, stop_time, interval_col=interval_col, interval_step=interval_step, ...)

  #' @noRd
  select_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
    cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
    cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval) %>%
    tibble::as_tibble()

  message(paste("Comparing abundances over",  nrow(time_contrasts), "timepoint contrasts"))
  relevant_comparisons = time_contrasts %>% as_tibble() %>%
    mutate(comp_abund = furrr::future_map2(.f = select_timepoints,
                                           .x = t1,
                                           .y = t2,
                                           interval_col=interval_col,
                                           timepoint_pred_df = timepoint_pred_df,
                                           .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                           .progress=TRUE))

  message(paste("Collecting relevant PLN network graph edges"))
  relevant_comparisons = relevant_comparisons %>%
    mutate(rec_edges = purrr::map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                  .x = comp_abund,
                                  ccm = ccm))

  # mutate(rec_edges = furrr::future_map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
  #                                      .x = comp_abund,
  #                                      ccm = ccm,
  #                                      .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
  #                                      .progress=TRUE))

  message(paste("Finding reciprocal pairs"))
  relevant_comparisons = relevant_comparisons %>%
    tidyr::unnest(rec_edges) %>%
    #dplyr::filter(pcor < 0) %>% # do we just want negative again?
    dplyr::filter((from_delta_log_abund > abs(min_pathfinding_lfc) & to_delta_log_abund < -abs(min_pathfinding_lfc)) |
                    (to_delta_log_abund > abs(min_pathfinding_lfc) & from_delta_log_abund < -abs(min_pathfinding_lfc))) %>%
    dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val)

  if (nrow(relevant_comparisons) == 0){
    stop ("No reciprocal pairs of nodes found")
  }

  edge_union = relevant_comparisons %>% select(from, to) %>% distinct()

  print (paste("finding shortest paths between ", nrow(edge_union), "pairs of nodes"))
  #print(head(relevant_comparisons))
  paths_for_relevant_edges = edge_union %>%
    mutate(path = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                                     .x = from, .y = to,
                                     pathfinding_graph,
                                     .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                     .progress=TRUE))

  cells_along_path_df = normalized_counts(ccm@ccs, "size_only", pseudocount = 0) %>%
    as.matrix() %>%
    Matrix::t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(!matches("rowname")) %>%
    rename(sample=rowname, cell_group=name, num_cells=value)

  paths_for_relevant_edges = paths_for_relevant_edges %>%
    filter(!is.na(path))

  if (nrow(paths_for_relevant_edges) == 0){
    stop ("No time-forward paths found")
  }

  print (paste("scoring", nrow(paths_for_relevant_edges), "paths for time flow..."))
  #print(head(paths_for_relevant_edges))

  paths_for_relevant_edges = paths_for_relevant_edges %>%
    mutate(time_vs_distance_model_stats = furrr::future_map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                            .x = path,
                                                            ccs=ccm@ccs,
                                                            cells_along_path_df=cells_along_path_df,
                                                            interval_col=interval_col,
                                                            .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                                            .progress=TRUE)) %>%
    tidyr::unnest(time_vs_distance_model_stats)
  #cells_along_path(ccs, paths_for_relevant_edges$path[[1]] %>% dplyr::select(from, to))


  print ("tabulating model stats...")
  #print(head(paths_for_relevant_edges))
  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(#time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq,
    time_dist_model_score = time_dist_model_ncells * time_dist_model_adj_rsq,
    path_score = time_dist_model_score)
  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method="BH"))

  #print (paths_for_relevant_edges)
  selected_paths = paths_for_relevant_edges %>% filter (time_dist_effect > 0 &
                                                          time_dist_effect_qval < q_val &
                                                          time_dist_model_adj_rsq > 0) %>%
    group_by(to) %>%
    #slice_max(dist_model_score, n=3) %>%
    ungroup() %>% arrange(desc(time_dist_model_score)) %>%
    dplyr::mutate(path_contrast="time")

  #selected_paths %>% select(origin=from, destination=to, path, path_score) %>% tidyr::unnest(path) %>% arrange(origin, destination) %>% filter(from %in% c("20", "38") & to %in% c("20", "38"))  %>% print(n=1000)

  #G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = FALSE)

  print ("combining paths...")
  G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles = TRUE)

  if (!is.null(G) && !is.na(G) && make_dag){

    #igraph::edge_attr(G, "support") = igraph::edge_attr(G, "total_perturb_path_score_supporting")
    cycle_breaking_scores = igraph::edge_attr(G, "total_path_score_supporting")
    cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0
    igraph::edge_attr(G, "support") = cycle_breaking_scores

    print (igraph::edge_attr(G, "support"))

    print ("breaking cycles...")
    G = hooke:::break_cycles_in_state_transition_graph(G)

    G = compute_min_path_cover(ccm, G, weight_attribute="support")

    edge_support = igraph::as_data_frame(G) %>% select(from, to)

    edge_support = left_join(edge_support,
                             selected_paths %>% select(-from, -to) %>% tidyr::unnest(path))

    #print (edge_support)
    edge_support = edge_support %>% group_by(from, to) %>% slice_max(time_dist_model_adj_rsq, n=1)

    #print (edge_support)
    G = igraph::graph_from_data_frame(edge_support, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))
  }

  print ("Finished building timeseries graph")
  return(G)
}

#' @noRd
get_paths_between_recip_time_nodes <- function(ccm,
                                               pathfinding_graph,
                                               timepoint_pred_df,
                                               timepoints,
                                               min_interval,
                                               max_interval,
                                               q_val,
                                               min_pathfinding_lfc,
                                               interval_col,
                                               verbose=FALSE){
  #timepoints = sort(unique(timepoint_pred_df[[interval_col]]))

  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  #old_omp_num_threads = single_thread_omp()
  #old_blas_num_threads = single_thread_blas()

  tryCatch({

    compare_timepoints <- function(timepoint_pred_df, t1, t2, interval_col)  {
      cond_x = timepoint_pred_df %>% filter(!!sym(interval_col) == t1)
      cond_y = timepoint_pred_df %>% filter(!!sym(interval_col) == t2)
      return(compare_abundances(ccm, cond_x, cond_y))
    }

    time_contrasts = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
      filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval)

    if (verbose)
      message(paste("\tLooking at", nrow(time_contrasts), "time contrasts"))

    # Find the nodes that undergo reciprocal fold changes between successive time points
    recip_time_node_pairs = time_contrasts %>%
      mutate(comp_abund = furrr::future_map2(.f = compare_timepoints,
                                             .x = t1,
                                             .y = t2,
                                             interval_col=interval_col,
                                             timepoint_pred_df = timepoint_pred_df,
                                             .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                             .progress=TRUE)) %>%
      mutate(rec_edges = furrr::future_map(.f = purrr::possibly(hooke:::collect_pln_graph_edges, NULL),
                                           .x = comp_abund,
                                           ccm = ccm,
                                           .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                           .progress=TRUE))
    recip_time_node_pairs = recip_time_node_pairs %>%
      tidyr::unnest(rec_edges)

    # Let's adjust all the p values to account for the fact that we've done many contrasts
    recip_time_node_pairs = recip_time_node_pairs %>%
      dplyr::mutate(from_delta_q_value = p.adjust(from_delta_p_value),
                    to_delta_p_value = p.adjust(to_delta_p_value))

    if (verbose){
      num_distinct_node_pairs = recip_time_node_pairs %>% select(from, to) %>% distinct()
      message(paste("\t Found", nrow(num_distinct_node_pairs), "reciprocal node pairs"))
    }

    recip_time_node_pairs = recip_time_node_pairs %>%
      #dplyr::filter(pcor < 0) %>% # do we just want negative again?
      dplyr::filter((from_delta_log_abund > abs(min_pathfinding_lfc) & to_delta_log_abund < -abs(min_pathfinding_lfc)) |
                      (to_delta_log_abund > abs(min_pathfinding_lfc) & from_delta_log_abund < -abs(min_pathfinding_lfc))) %>%
      dplyr::filter(from_delta_q_value < q_val & to_delta_q_value < q_val)

    if (verbose){
      num_distinct_node_pairs = recip_time_node_pairs %>% select(from, to) %>% distinct()
      message(paste("\tOf these", nrow(num_distinct_node_pairs), "survived thresholding"))
    }


    if (verbose)
      message("Computing shortest paths between reciprocal time nodes")

    # Compute the shortest paths between each of those node pairs in the pathfinding graph
    recip_time_node_pairs = recip_time_node_pairs %>% select(from, to) %>% distinct()
    paths_between_recip_time_nodes = recip_time_node_pairs %>%
      mutate(path = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                                       .x = from, .y = to,
                                       pathfinding_graph,
                                       .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                       .progress=TRUE))

    #paths_between_recip_time_nodes %>% filter (from == "15") %>% print

    cells_along_path_df = normalized_counts(ccm@ccs, "size_only", pseudocount = 0) %>%
      as.matrix() %>%
      Matrix::t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(!matches("rowname")) %>%
      rename(sample=rowname, cell_group=name, num_cells=value)

    # Score each shortest path for correlation between time and distance along it
    paths_between_recip_time_nodes = paths_between_recip_time_nodes %>%
      filter(!is.na(path))

    if (verbose)
      message(paste("Found", nrow(paths_between_recip_time_nodes), " paths through pathfinding graph. Scoring for time-flow..."))

    paths_between_recip_time_nodes = paths_between_recip_time_nodes %>%
      mutate(time_vs_distance_model_stats = furrr::future_map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                              .x = path,
                                                              ccs=ccm@ccs,
                                                              cells_along_path_df=cells_along_path_df,
                                                              interval_col=interval_col,
                                                              .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                                              .progress=TRUE)) %>%
      tidyr::unnest(time_vs_distance_model_stats)

    paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% mutate(time_dist_model_score = time_dist_effect * time_dist_model_adj_rsq)
    paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method="BH"))
    #paths_between_recip_time_nodes = paths_between_recip_time_nodes %>% filter (dist_effect > 0 & dist_effect_q_val < q_val & dist_model_adj_rsq > min_dist_vs_time_r_sq) %>%
    #  group_by(to) %>%
    #  #slice_max(dist_model_score, n=3) %>%
    #  ungroup() %>% arrange(desc(dist_model_score))
    return(paths_between_recip_time_nodes)
  }, error = function(e){
    print (e);
    return(NA)
  }, finally = {
    #RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    #RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  })
  return(NA)
}



#' @noRd
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
                                 min_pathfinding_lfc=0,
                                 verbose=FALSE,
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
                                                                      interval_col,
                                                                      min_pathfinding_lfc=min_pathfinding_lfc,
                                                                      verbose=verbose)
  return (paths_between_recip_time_nodes)
}

#' @noRd
select_paths_from_pathfinding_graph <- function(pathfinding_graph, selected_paths, allow_cycles=FALSE)
{

  support_tbl= selected_paths %>% dplyr::select(path_contrast, path_score, path) %>% tidyr::unnest(path)
  G = igraph::graph_from_data_frame(support_tbl %>% select(from, to) %>% distinct(), vertices=data.frame(id=igraph::V(pathfinding_graph)$name))
  edge_support_summary = support_tbl %>% group_by(from, to) %>% dplyr::select(from, to, path_score, path_contrast)  %>%
     distinct() %>% summarize(support=sum(path_score))
  edge_support_summary = edge_support_summary %>% mutate(support = ifelse(is.na(support), 0, support))
  annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% distinct() %>% left_join(edge_support_summary)
  annotated_G = igraph::graph_from_data_frame(annotated_G, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))
  return (annotated_G)

  # G = pathfinding_graph
  # G = igraph::delete_edges(G, igraph::E(G))
  #
  # for (i in (1:nrow(selected_paths))){
  #   next_path = selected_paths$path[[i]] %>% select(from, to)
  #   next_path_graph = next_path %>% igraph::graph_from_data_frame(directed=TRUE)
  #   G_prime = igraph::union(G, next_path_graph)
  #   if (allow_cycles | length(find_cycles(G_prime)) == 0){
  #     G = G_prime
  #   }else{
  #     # Debug:
  #     #print ("skipping due to cycle:")
  #     #print (selected_paths[i,] %>% select(-path))
  #   }
  # }
  # igraph::E(G)$transcriptome_dist = igraph::E(pathfinding_graph)[igraph::E(G)]$weight
  #
  # #print (selected_paths)
  # support_tbl= selected_paths %>% dplyr::select(path_contrast, path_score, path) %>% tidyr::unnest(path)
  # edge_support_summary = support_tbl %>% group_by(from, to) %>% dplyr::select(from, to, path_score, path_contrast)  %>%
  #   distinct() %>% summarize(support=sum(path_score))
  #
  # edge_support_summary = edge_support_summary %>% mutate(support = ifelse(is.na(support), 0, support))
  # #print (edge_support_summary %>% as.data.frame)
  #
  # #annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% left_join(edge_support)
  # annotated_G = igraph::as_data_frame(G) %>% as_tibble %>% dplyr::select(from, to) %>% distinct() %>% left_join(edge_support_summary)
  # annotated_G = igraph::graph_from_data_frame(annotated_G, directed=TRUE, vertices=data.frame(id=igraph::V(G)$name))
  #
  # return (annotated_G)
}

#' Break cycles in a state graph
#'
#' @export
#'
break_cycles_in_state_transition_graph <- function(state_graph, support_attribute)
{
  #print ("FAS weights")
  #print (igraph::E(state_graph)$support)

  cycle_breaking_scores = igraph::edge_attr(state_graph, support_attribute)
  cycle_breaking_scores[is.na(cycle_breaking_scores)] = 0

  # NOTE: squaring the scores should ensure that the cycle breaking procedure prefers to keep a few well supported edges even if that means removing lots of edges with lower support
  cycle_breaking_scores = cycle_breaking_scores^2
  #igraph::edge_attr(state_graph, "support") = cycle_breaking_scores

  #print (igraph::edge_attr(state_graph, "support"))

  fas = igraph::feedback_arc_set(state_graph, weights=cycle_breaking_scores)
  #print ("feedback arc set")
  #print (fas)
  state_graph = igraph::delete_edges(state_graph, fas)
  #print ("checking for cycles")
  #curr_cycles = find_cycles(state_graph)
  #if (length(curr_cycles) > 0){
    #print ("warning: cycles remain")
    #print (curr_cycles)
  #}
  return(state_graph)
}

#' @noRd
get_pcor_between_pair <- function(perturb_model_pcor_matrix, from, to){
  #edges_to_test = mapply(function(x, y) { return(c(x,y))}, from, to)
  pcor_vals = mapply(function(x, y) { return(perturb_model_pcor_matrix[x,y])}, from, to)
  #pcor_val = igraph::E(pcor_graph)[igraph::get.edge.ids(pcor_graph, edges_to_test)]$weight
  #pcor_val = as.vector(perturb_model_pcor_matrix[edges_to_test])
  return(pcor_vals)
}

score_path_for_loss = function(path, loss_tbl){
  cell_groups_on_path = tibble(cell_group=unique(c(path$from, path$to)))
  cell_groups_on_path = inner_join(cell_groups_on_path, loss_tbl, by="cell_group")

  # FIXME: we should find a way to incorporate effects from all the nodes along the path.
  cell_groups_on_path = cell_groups_on_path %>% filter(is_lost_when_present)
  if (nrow(cell_groups_on_path) > 0){
    cell_loss_stats_on_path = cell_groups_on_path %>%
      summarize(perturb_dist_effect = mean(loss_when_present, na.rm=TRUE),
                perturb_dist_effect_pval = mean(loss_when_present_q_val, na.rm=TRUE),
                perturb_dist_model_adj_rsq = NA,
                perturb_dist_model_ncells = NA)
  }else{
    cell_loss_stats_on_path = tibble(perturb_dist_effect = NA,
                                     perturb_dist_effect_pval = 1,
                                     perturb_dist_model_adj_rsq = NA,
                                     perturb_dist_model_ncells = NA)
  }


  return(cell_loss_stats_on_path)
}


#' @noRd
score_paths_for_perturbations <- function(perturbation_ccm,
                                          paths_between_recip_time_nodes,
                                          loss_tbl)
{
  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  #old_omp_num_threads = single_thread_omp()
  #old_blas_num_threads = single_thread_blas()

  tryCatch({
    paths_between_perturb_vs_wt_node_pairs = paths_between_recip_time_nodes %>%
      filter(!is.na(path))

    paths_between_perturb_vs_wt_node_pairs = paths_between_recip_time_nodes %>%
      #mutate(pcor_between_node_pair =  get_pcor_between_pair(perturb_model_pcor_matrix, from, to)) %>%
      mutate(perturb_vs_distance_model_stats = purrr::map(.f = purrr::possibly(score_path_for_loss, NA_character_),
                                                           .x = path,
                                                          loss_tbl#,
                                                                 #.options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                                                 #.progress=TRUE
                                                          )) %>%
      filter(!is.na(perturb_vs_distance_model_stats)) %>%
      tidyr::unnest(perturb_vs_distance_model_stats)
    paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_effect_qval = p.adjust(perturb_dist_effect_pval, method="BH"))
    return(paths_between_perturb_vs_wt_node_pairs)
  }, error = function(e) {
    print (e)
    return (NA)
  }, finally = {
    #RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    #RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  }
  )
  return (NA)
}

#' @noRd
compare_ko_to_wt_at_timepoint <- function(tp, perturbation_ccm, wt_pred_df, ko_pred_df, interval_col)  {
  cond_wt = wt_pred_df %>% filter(!!sym(interval_col) == tp)
  cond_ko = ko_pred_df %>% filter(!!sym(interval_col) == tp)
  return(compare_abundances(perturbation_ccm, cond_wt, cond_ko))
}

#' @noRd

estimate_loss_timing <- function(perturbation_ccm,
                                 start_time,
                                 stop_time,
                                 interval_step,
                                 control_ccm=perturbation_ccm,
                                 control_start_time=start_time,
                                 control_stop_time=stop_time,
                                 log_abund_detection_thresh=-5,
                                 delta_log_abund_loss_thresh=0,
                                 interval_col="timepoint",
                                 q_val=0.01,
                                 with_ties=FALSE,
                                 ...){

  fraction_of_presence_window_lost_thresh = 0.5

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

    extant_wt_tbl = get_extant_cell_types(control_ccm,
                                          control_start_time,
                                          control_stop_time,
                                          log_abund_detection_thresh=log_abund_detection_thresh,
                                          knockout=FALSE,
                                          ...)

    changes_when_present_in_wt = left_join(perturb_vs_wt_nodes %>% select(cell_group,
                                                                        wt_time_present=t1,
                                                                        delta_log_abund_when_present=delta_log_abund,
                                                                        delta_log_abund_when_present_se=delta_log_abund_se,
                                                                        log_abund_wt=log_abund_x,
                                                                        delta_q_value),
                                            extant_wt_tbl,
                                           by=c("cell_group"="cell_group", "wt_time_present"="timepoint")) %>%
      #mutate(peak_time_in_ctrl_within_perturb_time_range = tidyr::replace_na(peak_time_in_ctrl_within_perturb_time_range, FALSE)) %>%
      mutate(delta_q_value = ifelse(is.na(delta_q_value), 1, delta_q_value),
             delta_log_abund_when_present = ifelse(is.na(delta_log_abund_when_present), 0, delta_log_abund_when_present)) %>%
      mutate(is_lost_when_present = present_above_thresh & delta_log_abund_when_present < -abs(delta_log_abund_loss_thresh)) %>%
      mutate(is_gained_when_present = present_above_thresh & delta_log_abund_when_present > -abs(delta_log_abund_loss_thresh))
    #loss_when_present_in_wt = loss_when_present_in_wt %>% group_by(cell_group) %>% slice_min(peak_wt_time, n=1, with_ties=with_ties)

    # FIXME: maybe should refactor the code below into another function that summarizes contrast over intervals:
    # num samples
    n = nrow(model(perturbation_ccm)$fitted)
    # num parameters
    k = length(rownames(coef(perturbation_ccm@best_full_model)))
    df.r = n - k - 1
    #df_correction = sqrt(n / (n - k - 1))
    loss_summary_tbl = changes_when_present_in_wt %>%
      filter(is_lost_when_present) %>%
      group_by(cell_group) %>%
      summarize(loss_when_present = weighted.mean(delta_log_abund_when_present, percent_max_abund, na.rm=T),
                loss_when_present_se = weighted.mean(delta_log_abund_when_present_se, percent_max_abund, na.rm=T),
                loss_when_present_tvalue = weighted.mean(delta_log_abund_when_present/delta_log_abund_when_present_se, percent_max_abund, na.rm=T),
                loss_when_present_tvalue_df = df.r,
                loss_when_present_p_value = 2 * pt(-abs(loss_when_present_tvalue), loss_when_present_tvalue_df)) %>%
      mutate(loss_when_present_q_val = p.adjust(loss_when_present_p_value, method="bonferroni"))

    gain_summary_tbl = changes_when_present_in_wt %>%
      filter(is_gained_when_present) %>%
      group_by(cell_group) %>%
      summarize(gain_when_present = weighted.mean(delta_log_abund_when_present, percent_max_abund, na.rm=T),
                gain_when_present_se = weighted.mean(delta_log_abund_when_present_se, percent_max_abund, na.rm=T),
                gain_when_present_tvalue = weighted.mean(delta_log_abund_when_present/delta_log_abund_when_present_se, percent_max_abund, na.rm=T),
                gain_when_present_tvalue_df = df.r,
                gain_when_present_p_value = 2 * pt(-abs(gain_when_present_tvalue), gain_when_present_tvalue_df)) %>%
      mutate(gain_when_present_q_val = p.adjust(gain_when_present_p_value, method="bonferroni"))


    change_summary_tbl = changes_when_present_in_wt %>% select(cell_group) %>% distinct()
    change_summary_tbl = left_join(change_summary_tbl, loss_summary_tbl, by="cell_group")
    change_summary_tbl = left_join(change_summary_tbl, gain_summary_tbl, by="cell_group")
    change_summary_tbl = change_summary_tbl %>%
      mutate(loss_when_present_q_val = ifelse(is.na(loss_when_present_q_val), 1, loss_when_present_q_val),
             loss_when_present = ifelse(is.na(loss_when_present), NA, loss_when_present),
             is_lost_when_present = loss_when_present_q_val < q_val,
             gain_when_present_q_val = ifelse(is.na(gain_when_present_q_val), 1, gain_when_present_q_val),
             gain_when_present = ifelse(is.na(gain_when_present), NA, gain_when_present),
             is_gained_when_present = gain_when_present_q_val < q_val)

    peak_wt_abundance = estimate_abundances_over_interval(control_ccm, control_start_time, control_stop_time, interval_col=interval_col, knockout=FALSE, interval_step=interval_step, ...) %>%
          group_by(cell_group) %>% slice_max(log_abund, n=1) %>%
          select(cell_group, peak_wt_time=!!sym(interval_col))

    change_summary_tbl = left_join(change_summary_tbl,
                                   peak_wt_abundance,
                                   by=c("cell_group"))

    return(change_summary_tbl)
}

# estimate_loss_timing <- function(perturbation_ccm,
#                                  start_time,
#                                  stop_time,
#                                  interval_step,
#                                  control_ccm=perturbation_ccm,
#                                  control_start_time=start_time,
#                                  control_stop_time=stop_time,
#                                  log_abund_detection_thresh=-5,
#                                  delta_log_abund_loss_thresh=0,
#                                  interval_col="timepoint",
#                                  q_val=0.01,
#                                  with_ties=FALSE,
#                                  ...){
#
#   fraction_of_presence_window_lost_thresh = 0.5
#
#     wt_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=FALSE, interval_col=interval_col, interval_step=interval_step, ...)
#     ko_timepoint_pred_df = estimate_abundances_over_interval(perturbation_ccm, start_time, stop_time, knockout=TRUE, interval_col=interval_col, interval_step=interval_step, ...)
#
#     #print (wt_timepoint_pred_df)
#     #print (ko_timepoint_pred_df)
#     timepoints = seq(start_time, stop_time, interval_step)
#
#     # Find the pairs of nodes that are both lost in the perturbation at the same time
#     perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
#       mutate(comp_abund = purrr::map(.f = compare_ko_to_wt_at_timepoint,
#                                             .x = t1,
#                                             perturbation_ccm=perturbation_ccm,
#                                             interval_col=interval_col,
#                                             wt_pred_df = wt_timepoint_pred_df,
#                                             ko_pred_df = ko_timepoint_pred_df)) %>% tidyr::unnest(comp_abund)
#
#     earliest_loss_tbl =  perturb_vs_wt_nodes %>%
#       mutate(delta_q_value = p.adjust(delta_p_value, method="BH")) %>%
#       filter(delta_log_abund < -abs(delta_log_abund_loss_thresh) & delta_q_value <= q_val)
#
#     #print (earliest_loss_tbl)
#
#
#
#     #peak_wt_abundance = perturb_vs_wt_nodes %>% group_by(cell_group) %>% slice_max(log_abund_x, n=1, with_ties=with_ties)
#
#     extant_wt_tbl = get_extant_cell_types(control_ccm,
#                                           control_start_time,
#                                           control_stop_time,
#                                           log_abund_detection_thresh=log_abund_detection_thresh,
#                                           knockout=FALSE,
#                                           ...)
#
#     loss_when_present_in_wt = left_join(perturb_vs_wt_nodes %>% select(cell_group,
#                                                                         wt_time_present=t1,
#                                                                         delta_log_abund_when_present=delta_log_abund,
#                                                                         delta_log_abund_when_present_se=delta_log_abund_se,
#                                                                         log_abund_wt=log_abund_x,
#                                                                         delta_q_value),
#                                             extant_wt_tbl,
#                                            by=c("cell_group"="cell_group", "wt_time_present"="timepoint")) %>%
#       #mutate(peak_time_in_ctrl_within_perturb_time_range = tidyr::replace_na(peak_time_in_ctrl_within_perturb_time_range, FALSE)) %>%
#       mutate(delta_q_value = ifelse(is.na(delta_q_value), 1, delta_q_value),
#              delta_log_abund_when_present = ifelse(is.na(delta_log_abund_when_present), 1, delta_log_abund_when_present)) %>%
#       mutate(is_lost_when_present = present_above_thresh & delta_log_abund_when_present < -abs(delta_log_abund_loss_thresh) & delta_q_value <= q_val,
#              loss_when_present_time = ifelse(is_lost_when_present, wt_time_present, NA)) %>%
#       mutate(is_gained_when_present = present_above_thresh & delta_log_abund_when_present > -abs(delta_log_abund_loss_thresh) & delta_q_value <= q_val,
#              gain_when_present_time = ifelse(is_gained_when_present, wt_time_present, NA))
#     #loss_when_present_in_wt = loss_when_present_in_wt %>% group_by(cell_group) %>% slice_min(peak_wt_time, n=1, with_ties=with_ties)
#
#     loss_when_present_in_wt = loss_when_present_in_wt %>% select(cell_group,
#                                                                  present_above_thresh,
#                                                                      #peak_time_in_ctrl_within_perturb_time_range,
#                                                                  wt_time_present,
#                                                                  log_abund_wt,
#                                                                  delta_log_abund_when_present,
#                                                                  delta_log_abund_when_present_se,
#                                                                  delta_q_value,
#                                                                  is_lost_when_present,
#                                                                  loss_when_present_time,
#                                                                  is_gained_when_present,
#                                                                  gain_when_present_time) %>% distinct()
#     loss_when_present_in_wt = loss_when_present_in_wt %>% group_by(cell_group) %>%
#       filter(present_above_thresh) %>%
#       summarize(loss_when_present_q_val = min(p.adjust(delta_q_value[is_lost_when_present & log_abund_wt == max(log_abund_wt)]), na.rm=TRUE),
#                 #is_lost_when_present = sum(is_lost_when_present) / sum(present_above_thresh) > fraction_of_presence_window_lost_thresh,
#                 is_lost_when_present = sum(is_lost_when_present) > 0,
#                 largest_loss_when_present_in_wt = min(delta_log_abund_when_present, na.rm=TRUE),
#                 earliest_loss_when_present_in_wt = min(loss_when_present_time, na.rm=TRUE),
#                 latest_loss_when_present_in_wt = max(loss_when_present_time, na.rm=TRUE),
#                 gain_when_present_q_val = min(p.adjust(delta_q_value[is_gained_when_present & log_abund_wt == max(log_abund_wt)]), na.rm=TRUE),
#                 is_gained_when_present = sum(is_gained_when_present) > 0,
#                 largest_gain_when_present_in_wt = max(delta_log_abund_when_present, na.rm=TRUE),
#                 earliest_gain_when_present_in_wt = min(gain_when_present_time, na.rm=TRUE),
#                 latest_gain_when_present_in_wt = max(gain_when_present_time, na.rm=TRUE)) %>%
#       mutate(loss_when_present_q_val = ifelse(is.infinite(loss_when_present_q_val), 1, loss_when_present_q_val),
#              gain_when_present_q_val = ifelse(is.infinite(gain_when_present_q_val), 1, gain_when_present_q_val),
#              largest_loss_when_present_in_wt = ifelse(is.infinite(largest_loss_when_present_in_wt), NA, largest_loss_when_present_in_wt),
#              earliest_loss_when_present_in_wt = ifelse(is.infinite(earliest_loss_when_present_in_wt), NA, earliest_loss_when_present_in_wt),
#              latest_loss_when_present_in_wt = ifelse(is.infinite(latest_loss_when_present_in_wt), NA, latest_loss_when_present_in_wt),
#              largest_gain_when_present_in_wt = ifelse(is.infinite(largest_gain_when_present_in_wt), NA, largest_gain_when_present_in_wt),
#              earliest_gain_when_present_in_wt = ifelse(is.infinite(earliest_gain_when_present_in_wt), NA, earliest_gain_when_present_in_wt),
#              latest_gain_when_present_in_wt = ifelse(is.infinite(latest_gain_when_present_in_wt), NA, latest_gain_when_present_in_wt))
#
#     peak_wt_abundance = estimate_abundances_over_interval(control_ccm, control_start_time, control_stop_time, interval_col=interval_col, knockout=FALSE, interval_step=interval_step, ...) %>%
#       group_by(cell_group) %>% slice_max(log_abund, n=1) %>%
#       mutate(peak_time_in_ctrl_within_perturb_time_range = !!sym(interval_col) >= start_time & !!sym(interval_col) <= stop_time) %>%
#       select(cell_group, !!sym(interval_col), peak_time_in_ctrl_within_perturb_time_range)
#
#     loss_at_peak_wt_abundance = right_join(perturb_vs_wt_nodes %>% select(cell_group,
#                                                                  peak_wt_time=t1,
#                                                                  delta_log_abund_at_peak=delta_log_abund,
#                                                                  delta_q_value),
#                                   peak_wt_abundance,
#                                   by=c("cell_group", "peak_wt_time"=interval_col)) %>%
#       mutate(peak_time_in_ctrl_within_perturb_time_range = tidyr::replace_na(peak_time_in_ctrl_within_perturb_time_range, FALSE)) %>%
#       mutate(is_lost_at_peak = peak_time_in_ctrl_within_perturb_time_range & delta_log_abund_at_peak < -abs(delta_log_abund_loss_thresh) & delta_q_value <= q_val,
#              peak_loss_time = ifelse(is_lost_at_peak, peak_wt_time, NA))
#     loss_at_peak_wt_abundance = loss_at_peak_wt_abundance %>% group_by(cell_group) %>% slice_min(peak_wt_time, n=1, with_ties=with_ties)
#
#     loss_at_peak_wt_abundance = loss_at_peak_wt_abundance %>% select(cell_group,
#                                                                      peak_time_in_ctrl_within_perturb_time_range,
#                                                                      peak_wt_time,
#                                                                      delta_log_abund_at_peak,
#                                                                      is_lost_at_peak,
#                                                                      peak_loss_time)
#
#     #loss_at_peak_wt_abundance %>% filter(cell_group == "1") %>% print()
#
#     largest_losses = perturb_vs_wt_nodes %>%
#       filter(delta_log_abund < -abs(delta_log_abund_loss_thresh) & delta_q_value <= q_val) %>%
#       group_by(cell_group) %>%
#       #arrange(!!sym(paste(interval_col, "_x", sep="")), )
#       slice_min(delta_log_abund, n=1, with_ties=with_ties) %>% select(cell_group, largest_loss_time=!!sym(paste(interval_col, "_x", sep="")), largest_loss=delta_log_abund)
#
#     #print (largest_losses)
#     earliest_loss_tbl = earliest_loss_tbl %>%
#       group_by(cell_group) %>%
#       #mutate(largest_loss = min(delta_log_abund),
#       #       largest_loss_time = !!sym(paste(interval_col, "_x", sep=""))[which(delta_log_abund == largest_loss)]) %>%
#       summarize(#largest_loss_time = min(largest_loss_time),
#         earliest_time = min(!!sym(paste(interval_col, "_x", sep=""))),
#         latest_time = max(!!sym(paste(interval_col, "_x", sep=""))))
#
#     loss_tbl = loss_when_present_in_wt
#     loss_tbl = loss_tbl %>% left_join(largest_losses, by="cell_group")
#     loss_tbl = loss_tbl %>% left_join(loss_at_peak_wt_abundance, by="cell_group")
#     loss_tbl = loss_tbl %>% left_join(earliest_loss_tbl, by="cell_group")
#
#     #loss_tbl = earliest_loss_tbl %>% left_join(largest_losses, by="cell_group")
#     #loss_tbl = loss_tbl %>% left_join(loss_at_peak_wt_abundance, by="cell_group")
#     #loss_tbl = loss_tbl %>% left_join(loss_when_present_in_wt, by="cell_group")
#
#     return(loss_tbl)
# }

#' Identify the possible origins for each destination
#'
#' @noRd
# build_perturbation_transition_dag <- function(ccm,
#                                               extant_cell_type_df,
#                                               pathfinding_graph,
#                                               q_val=0.01,
#                                               start_time = NULL,
#                                               stop_time = NULL,
#                                               perturbation_col="knockout",
#                                               interval_col="timepoint",
#                                               interval_step = 2,
#                                               min_interval = 4,
#                                               max_interval = 24,
#                                               min_pathfinding_lfc=0,
#                                               ...) {
#
#   paths_between_recip_time_nodes = get_timeseries_paths(ccm,
#                                                         extant_cell_type_df,
#                                                         pathfinding_graph,
#                                                         start_time = start_time,
#                                                         stop_time = stop_time,
#                                                         interval_col=interval_col,
#                                                         interval_step = interval_step,
#                                                         min_interval = min_interval,
#                                                         max_interval = max_interval,
#                                                         min_pathfinding_lfc = min_pathfinding_lfc,
#                                                         knockout=FALSE,
#                                                         ...)
#
#
#
#   paths_between_perturb_vs_wt_node_pairs = score_paths_for_perturbations(ccm, paths_between_recip_time_nodes)
#   paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(path_score = perturb_dist_model_ncells * perturb_dist_model_adj_rsq)
#   paths_between_perturb_vs_wt_node_pairs = paths_between_perturb_vs_wt_node_pairs %>% mutate(perturb_dist_effect_qval = p.adjust(perturb_dist_effect_pval, method="BH"))
#
#   #origin_dest_node_pairs = dplyr::intersect(perturb_vs_wt_node_pairs, recip_time_node_pairs)
#
#   origin_dest_node_pairs = dplyr::inner_join(paths_between_recip_time_nodes, paths_between_perturb_vs_wt_node_pairs)
#   #print (origin_dest_node_pairs)
#   selected_paths = origin_dest_node_pairs %>% filter (perturb_dist_effect < 0 & perturb_dist_effect_q_val < q_val &
#                                                         time_dist_effect > 0 & time_dist_effect_q_val < q_val) %>%
#     group_by(to) %>%
#     #slice_max(dist_model_score, n=3) %>%
#     ungroup() %>% arrange(desc(path_score)) %>%
#     mutate(path_contrast="time_perturb")
#
#
#   G = select_paths_from_pathfinding_graph(pathfinding_graph, selected_paths, allow_cycles=FALSE)
#   return (G)
# }


#' @noRd
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
    mutate(chain_leg = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag,
                                   .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                   .progress=TRUE))

  covered_graph = paths_from_chains %>% select(chain_leg) %>%
    tidyr::unnest(c(chain_leg)) %>%
    igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(G_tr)$name))

  chain_heads = names(which(igraph::degree(covered_graph, mode="in") == 0))
  chain_tails = names(which(igraph::degree(covered_graph, mode="out") == 0))

  source_edge_df = data.frame(from="source", to=chain_heads)
  sink_edge_df = data.frame(from=chain_tails, to="sink")


  paths_to_chain_heads = source_edge_df %>% as_tibble() %>%
    mutate(chain_leg = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag,
                                   .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                   .progress=TRUE))

  paths_to_chain_tails = sink_edge_df %>% as_tibble() %>%
    mutate(chain_leg = furrr::future_map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                                   .x = from, .y = to,
                                   node_dag,
                                   .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                   .progress=TRUE))

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
#' @noRd
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
                                            min_pathfinding_lfc=0,
                                            make_dag=FALSE,
                                            links_between_components=c("ctp", "none", "strongest-pcor", "strong-pcor"),
                                            components = "partition",
                                            edge_whitelist=NULL,
                                            edge_blacklist=NULL,
                                            ...){

  message("Determining extant cell types")
  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start_time,
                                              stop_time,
                                              interval_col=interval_col,
                                              log_abund_detection_thresh=log_abund_detection_thresh,
                                              ...)

  # Now let's set up a directed graph that links the states between which cells *could* directly
  # transition. If we don't know the direction of flow, add edges in both directions. The idea is
  # that we will find shortest paths over this graph between destination states and their plausible
  # origin states, then choose the best origins for each destination.

  message("Initializing pathfinding graph")
  pathfinding_graph = init_pathfinding_graph(ccm,
                                             extant_cell_type_df,
                                             links_between_components=links_between_components,
                                             components = components,
                                             edge_whitelist=edge_whitelist,
                                             edge_blacklist=edge_blacklist)


  G = build_timeseries_transition_graph(ccm,
                                          extant_cell_type_df,
                                          pathfinding_graph,
                                          q_val,
                                          start_time,
                                          stop_time,
                                          interval_col,
                                          interval_step,
                                          min_interval,
                                          max_interval,
                                          min_pathfinding_lfc=min_pathfinding_lfc,
                                          make_dag=make_dag,
                                          ...)

    # FIXME: Consider moving this step into build_timeseries_transition_graph()?
    if (!is.null(G))
      igraph::V(G)$cell_group = ccm@ccs@info$cell_group

    return(G)
}

#' helper function for assessing the gains and losses of each cell type in a
#' perturbation
#'
#' @noRd
collect_perturb_effects = function(perturbation_ccm,
                                   time_window,
                                   control_ccm,
                                   control_time_window,
                                   interval_col,
                                   q_val=0.01,
                                   batch_col="expt",
                                   interval_step = 2,
                                   min_interval = 4,
                                   max_interval = 24,
                                   log_abund_detection_thresh=-5,
                                   min_lfc=0,
                                   ...) {
  start_time = min(as.numeric(time_window$start_time))
  stop_time = min(as.numeric(time_window$stop_time))

  control_start_time = min(as.numeric(control_time_window$start_time))
  control_stop_time = min(as.numeric(control_time_window$stop_time))

  #print (start_time)
  #print (stop_time)
  #message ("\tEstimating loss timing")
  perturb_effect_summary = hooke:::estimate_loss_timing(perturbation_ccm,
                                                        start_time=start_time,
                                                        stop_time=stop_time,
                                                        interval_step = interval_step,
                                                        interval_col=interval_col,
                                                        control_ccm=control_ccm,
                                                        control_start_time=control_start_time,
                                                        control_stop_time=control_stop_time,
                                                        log_abund_detection_thresh=log_abund_detection_thresh,
                                                        q_val = q_val,
                                                        delta_log_abund_loss_thresh=min_lfc,
                                                        ...)
  return (perturb_effect_summary)
}

#' Assess the cell type gains and losses following a perturbation.
#'
#' Returns a table of summary stats for gains and losses that's used to assemble
#' the cell types into a depencency graph
#'
#' @export
assess_perturbation_effects = function(control_timeseries_ccm,
                                       perturbation_ccm_tbl,
                                       q_val=0.01,
                                       start_time = NULL,
                                       stop_time = NULL,
                                       perturbation_col="knockout",
                                       interval_col="timepoint",
                                       batch_col="expt",
                                       interval_step = 2,
                                       min_interval = 4,
                                       max_interval = 24,
                                       log_abund_detection_thresh=-5,
                                       min_lfc=0,
                                       verbose=FALSE,
                                       ...)
{

  perturbation_ccm_tbl = perturbation_ccm_tbl %>%
    dplyr::mutate(perturb_summary_tbl = purrr::map2(.f = purrr::possibly(
      collect_perturb_effects, NA_real_),
      .x = perturb_ccm,
      .y = perturb_time_window,
      pathfinding_graph,
      control_ccm=control_timeseries_ccm,
      control_time_window=tibble(start_time=start_time, stop_time=stop_time),
      interval_col=interval_col,
      batch_col=batch_col,
      interval_step = interval_step,
      q_val=q_val,
      min_lfc=min_lfc,
      log_abund_detection_thresh=log_abund_detection_thresh,
      ...))

  # Perform a global correction for multiple testing
  perturbs = perturbation_ccm_tbl %>%
    dplyr::select(perturb_name, perturb_summary_tbl)

  # Start from the q values, as these are already corrected for the number of cell types in the model
  perturbs = perturbs %>% tidyr::unnest(cols= c(perturb_summary_tbl))%>%
    ungroup() %>%
    mutate(loss_when_present_q_val = p.adjust(loss_when_present_q_val, method="bonferroni"),
           loss_when_present_q_val = ifelse(is.na(loss_when_present_q_val), 1, loss_when_present_q_val),
           is_lost_when_present = loss_when_present_q_val < q_val,
           gain_when_present_q_val = p.adjust(gain_when_present_q_val, method="bonferroni"),
           gain_when_present_q_val = ifelse(is.na(gain_when_present_q_val), 1, gain_when_present_q_val),
           is_gained_when_present = gain_when_present_q_val < q_val)
  perturbs = perturbs %>% tidyr::nest(perturb_summary_tbl = !perturb_name)
  perturbation_ccm_tbl$perturb_summary_tbl = perturbs$perturb_summary_tbl

  return (perturbation_ccm_tbl)
}

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
                                                         batch_col="expt",
                                                         interval_step = 2,
                                                         min_interval = 4,
                                                         max_interval = 24,
                                                         log_abund_detection_thresh=-5,
                                                         min_pathfinding_lfc=0,
                                                         links_between_components=c("ctp", "none", "strongest-pcor", "strong-pcor"),
                                                         components = "partition",
                                                         verbose=FALSE,
                                                         edge_whitelist=NULL,
                                                         edge_blacklist=NULL,
                                                         ...)
{

  # Temporarily set the number of threads OpenMP & the BLAS library can use to be 1
  #old_omp_num_threads = single_thread_omp()
  #old_blas_num_threads = single_thread_blas()

  tryCatch({

    # Get a table of the cell types that are in the control
    # FIXME: "knockout" is hard coded and should be a user-defined term in the model
    extant_cell_type_df = get_extant_cell_types(control_timeseries_ccm,
                                                start_time,
                                                stop_time,
                                                interval_col=interval_col,
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
                                               links_between_components=links_between_components,
                                               components = components,
                                               edge_whitelist=edge_whitelist,
                                               edge_blacklist=edge_blacklist)

    timeseries_graph = assemble_timeseries_transitions(control_timeseries_ccm,
                                                       q_val=q_val,
                                                       start_time = start_time,
                                                       stop_time = stop_time,
                                                       perturbation_col=perturbation_col,
                                                       interval_col=interval_col,
                                                       interval_step = interval_step,
                                                       min_interval = min_interval,
                                                       max_interval = max_interval,
                                                       log_abund_detection_thresh=log_abund_detection_thresh,
                                                       min_pathfinding_lfc=min_pathfinding_lfc,
                                                       links_between_components=links_between_components,
                                                       edge_whitelist=edge_whitelist,
                                                       edge_blacklist=edge_blacklist,
                                                       ...)

    if (is.null(timeseries_graph) || is.na(timeseries_graph)){
      stop("Error: timeseries graph assembly failed. Aborting.")
    }

    if (is.null(perturbation_ccm_tbl$perturb_summary_tbl)){
      perturbation_ccm_tbl = assess_perturbation_effects(control_timeseries_ccm,
                                                         perturbation_ccm_tbl,
                                                         q_val=q_val,
                                                         start_time = start_time,
                                                         stop_time = stop_time,
                                                         perturbation_col=perturbation_col,
                                                         interval_col=interval_col,
                                                         batch_col=batch_col,
                                                         interval_step = interval_step,
                                                         min_interval = min_interval,
                                                         max_interval = max_interval,
                                                         log_abund_detection_thresh=log_abund_detection_thresh,
                                                         min_lfc=min_pathfinding_lfc,
                                                         verbose=verbose,
                                                         ...)
    }

    pathfinding_graph = igraph::intersection(pathfinding_graph, timeseries_graph)

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
        .y = perturb_summary_tbl,
        pathfinding_graph))

    # path_tbl = perturbation_ccm_tbl %>%
    #   dplyr::mutate(paths_between_concordant_loss_nodes = purrr::map2(.f = purrr::possibly(
    #     get_perturbation_paths, NA_real_),
    #     .x = perturb_ccm,
    #     .y = perturb_time_window,
    #     pathfinding_graph,
    #     control_ccm=control_timeseries_ccm,
    #     control_time_window=tibble(start_time=start_time, stop_time=stop_time),
    #     interval_col=interval_col,
    #     batch_col=batch_col,
    #     interval_step = interval_step,
    #     q_val=q_val,
    #     min_pathfinding_lfc=min_pathfinding_lfc,
    #     log_abund_detection_thresh=log_abund_detection_thresh,
    #     ...))

    path_tbl = path_tbl %>%
      filter(!is.na(paths_between_concordant_loss_nodes)) %>%
      select(-perturb_summary_tbl) %>%
      tidyr::unnest(paths_between_concordant_loss_nodes) %>%
      filter(!is.na(path))

    message (paste("Found ", nrow(path_tbl), " loss paths"))

    if (nrow(path_tbl) == 0){
      stop("No loss paths")
    }


    cells_along_path_df = normalized_counts(control_timeseries_ccm@ccs, "size_only", pseudocount = 0) %>%
      as.matrix() %>%
      Matrix::t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(!matches("rowname")) %>%
      rename(sample=rowname, cell_group=name, num_cells=value)

    if (verbose)
      message (paste("Assessing time flows across", nrow(path_tbl), " loss paths"))

    # Measure the flow of time along each perturbation loss path
    path_tbl = path_tbl %>%
      mutate(time_vs_distance_model_stats = furrr::future_map(.f = purrr::possibly(measure_time_delta_along_path, NA_character_),
                                                       .x = path,
                                                       ccs=control_timeseries_ccm@ccs,
                                                       cells_along_path_df=cells_along_path_df,
                                                       interval_col=interval_col,
                                                       .options=furrr::furrr_options(stdout=FALSE,conditions = character()),
                                                       .progress=TRUE)) %>%
      tidyr::unnest(time_vs_distance_model_stats)

    # Compute some scores that summarize the level of support for each path by the various perturbations.
    path_tbl = path_tbl %>% mutate(time_dist_effect_qval = p.adjust(time_dist_effect_pval, method="BH"))
    path_tbl = path_tbl %>% mutate(timeseries_model_path_score = ifelse(time_dist_effect > 0 & time_dist_effect_qval < q_val, time_dist_model_adj_rsq, 0))
    path_tbl = path_tbl %>% mutate(perturb_model_path_score = ifelse(perturb_dist_effect < 0 & perturb_dist_effect_qval < q_val, abs(perturb_dist_effect), 0))
    path_tbl = path_tbl %>% mutate(path_score = ifelse(timeseries_model_path_score > 0 & perturb_model_path_score > 0, perturb_model_path_score * timeseries_model_path_score, 0))

    #path_tbl = path_tbl %>% mutate(path_score = ifelse(timeseries_model_path_score > 0, timeseries_model_path_score, 0))

    # Exclude paths that have no support from perturbations or go against the flow of time
    path_tbl =  path_tbl %>% filter (path_score > 0)

    message (paste("Found ", nrow(path_tbl), "significant loss paths with forward time-flow"))

    if (nrow(path_tbl) == 0){
      stop("No significant loss paths")
    }


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
                                            ...)


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

    G_label = igraph::edge_attr(G, "support_label")
    G_label[is.na(G_label)] = ""
    igraph::edge_attr(G, "support_label") = G_label

    return (G)
  }, error = function(e){
    print (e)
    return(NA)
  }, finally = {
    #RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    #RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  }
  )

  return (NA)
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
                                                ...)
{
  # Flatten all the paths that the perturbation assembler used to link up the
  # cell states
  path_score_tbl = perturbation_path_tbl %>% rename(origin=from, destination=to) %>% tidyr::unnest(path)

  #print (path_score_tbl)

  #edge_support = left_join(edge_support, path_score_tbl)

  # Uplift some statistics collected during assembly to score edges based on
  # strength of the supporting evidence
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

  # Which perturbations support each edge?
  edge_support_perturbs = path_score_tbl %>%
    ungroup () %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(from, to, perturb_name)  %>%
    group_by(from, to) %>%
    distinct() %>% tidyr::nest(supporting_perturbs = c(perturb_name))

  # Which perturbations lead to a direct loss of each node?
  # Direct loss is defined as either:
  # - A node that is lost at its peak in the WT (and where the perturbation time series covers that peak) and NOT part of a loss path
  # - The first node in a loss path (both the first and last nodes in a loss path are lost at peak WT abundance)
  # FIXME: right now the code below does not detect the first category of direct losses. Fixing this
  # will require that we pass in the arguments needed to call estimate_loss_timing()
  node_direct_perturbs = path_score_tbl %>%
    ungroup () %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(id=origin, perturb_name) %>% distinct()

  # Which perturbations lead to an indirect loss of each node?
  # Indirect loss is all nodes on loss paths except the first.

  node_indirect_perturbs = path_score_tbl %>% group_by(origin, destination) %>%
    ungroup () %>%
    dplyr::select(origin, id=from, perturb_name)
  node_indirect_perturbs = node_indirect_perturbs %>% rbind(
    path_score_tbl %>%
      ungroup () %>%
      dplyr::select(origin, id=to, perturb_name)
  )
  node_indirect_perturbs = node_indirect_perturbs %>% filter(origin != id) %>% select(id, perturb_name)
  node_indirect_perturbs = node_indirect_perturbs %>% distinct()

  # If a node is the first in a path, but internal to another, count it as an indirect loss:
  node_direct_perturbs = dplyr::setdiff(node_direct_perturbs, node_indirect_perturbs)

  node_indirect_perturbs = node_indirect_perturbs %>%
    group_by(id) %>%
    tidyr::nest(indirect_perturb=perturb_name)

  node_direct_perturbs = node_direct_perturbs %>%
    group_by(id) %>%
    tidyr::nest(direct_perturb=perturb_name)

  node_metadata = tibble(id=igraph::V(state_transition_graph)$name)
  node_metadata = node_metadata %>% left_join(node_direct_perturbs)
  node_metadata = node_metadata %>% left_join(node_indirect_perturbs)

  # TODO: we should also annotate nodes as having "inferred" dependency on a given
  # perturbation if its all its ancestors depend on that perturbation, but the node's
  # peak abundance is outside the window of our collected measurements for that
  # perturbation

  # Construct nice labels for labeling edges according to support
  # FIXME: we should probably just move this to its own function to be
  # used in plotting functions
  edge_support_labels = path_score_tbl %>%
    ungroup () %>%
    arrange(desc(perturb_model_path_score)) %>%
    dplyr::select(from, to, perturb_name)  %>%
    group_by(from, to) %>%
    distinct() %>%
    summarize(edge_name = stringr::str_c(from, to, sep="~"),
              #supporting_perturbs = perturb_name,
              support_label=ifelse(n() > 3, paste0(c(perturb_name[1:3], paste("+", n()-3, " more", sep="")), collapse = "\n"),
                                   paste0(perturb_name, collapse = "\n")))

  #print (edge_support_labels)
  edge_support_summary = edge_support_summary %>% left_join(edge_support_labels)
  edge_support_summary = edge_support_summary %>% left_join(edge_support_perturbs)

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
  annotated_state_transition_graph = igraph::graph_from_data_frame(annotated_state_transition_graph, directed=TRUE, vertices=node_metadata)
  return(annotated_state_transition_graph)
}

#' Simplify a directed state transition graph by grouping nodes according to a specified label
#' assumes the graphs over nodes corresponding to groups in the ccm
#' @export
contract_state_graph <- function(ccs,
                                 state_graph,
                                 group_nodes_by,
                                 edge_attr_policy=list("weight"="sum",
                                                       "name"="concat",
                                                       "num_perturbs_supporting"="sum",
                                                       "max_timeseries_path_score_supporting"="sum",
                                                       "total_timeseries_path_score_supporting"="sum",
                                                       "max_perturb_path_score_supporting"="sum",
                                                       "total_perturb_path_score_supporting"="sum",
                                                       "max_path_score_supporting"="sum",
                                                       "total_path_score_supporting"="sum",
                                                       "edge_name"="ignore",
                                                       "support_label"="concat",
                                                       "supporting_perturbs"="concat",
                                                       "ignore")){
  # Create simplified cell state graph just on cell type (not cluster):
  cell_groups = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccs@cds) %>%
    as.data.frame %>% select(!!sym(group_nodes_by))

  cell_group_metadata$cell_group = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
    as.data.frame %>%
    dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
    dplyr::group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
  colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
  node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  node_metadata = node_metadata %>% mutate(id = as.numeric(id)) %>% arrange(id)
  contraction_mapping = as.factor(node_metadata$group_nodes_by)
  contraction_mapping_names = as.character(levels(contraction_mapping))
  contraction_mapping = as.numeric(contraction_mapping)
  names(contraction_mapping) = node_metadata$id
  contracted_state_graph = igraph::contract(state_graph, mapping=contraction_mapping, vertex.attr.comb="ignore")
  igraph::V(contracted_state_graph)$name = unlist(contraction_mapping_names[as.numeric(igraph::V(contracted_state_graph))])
  contracted_state_graph = igraph::simplify(contracted_state_graph,
                                            edge.attr.comb=edge_attr_policy)
  igraph::V(contracted_state_graph)$cell_group = group_nodes_by
  return(contracted_state_graph)
}

