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