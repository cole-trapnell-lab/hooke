library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)


setwd("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/mesoderm")

meso_cds = readRDS("final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")

#meso_cds = cluster_cells(meso_cds, resolution=1e-4, random_seed=42)

colData(meso_cds)$cluster = as.character(clusters(meso_cds))
meso_cds = meso_cds[, is.na(colData(meso_cds)$embryo) == FALSE]
colData(meso_cds)$cell_type = colData(meso_cds)$cell_type_sub

wt_ccs = new_cell_count_set(meso_cds,
                            sample_group = "embryo",
                            cell_group = "cluster")

#
# wt_start = 18
# wt_stop = 96
# num_time_breaks = 5
# time_breakpoints = seq(wt_start, wt_stop, length.out=num_time_breaks)
# time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
# wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")



initial_pcor_graph = function(ccs) {
  paga_edges = get_paga_graph(wt_ccs@cds) %>% igraph::as_data_frame() %>% as_tibble()
  return(paga_edges)
}

build_interval_formula <- function(ccs, num_breaks, interval_var="timepoint", interval_start=NULL, interval_stop=NULL){
  if (is.null(interval_start)){
    interval_start = as.numeric(min(colData(ccs@cds)[,interval_var]))
  }
  if (is.null(interval_stop)){
    interval_stop = as.numeric(max(colData(ccs@cds)[,interval_var]))
  }

  interval_breakpoints = seq(interval_start, interval_stop, length.out=num_breaks)
  interval_breakpoints = interval_breakpoints[2:(length(interval_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
  interval_formula_str = paste("~ splines::ns(", interval_var, ", knots=", paste("c(",paste(interval_breakpoints, collapse=","), ")", sep=""), ")")
  return(interval_formula_str)
}
#debug(build_interval_formula)

wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_start=18, interval_stop=48, num_breaks=5)

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 #main_model_formula_str = "~expt",
                                 main_model_formula_str = wt_main_model_formula_str,
                                 #nuisance_model_formula_str = "~1",
                                 nuisance_model_formula_str = "~expt",
                                 whitelist = initial_pcor_graph(wt_ccs))

wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=3)

estimate_abundances_over_interval <- function(ccm, start, stop, interval_step=2, reference_experiment="GAP16") {

  timepoint_pred_df = tibble(timepoint=seq(start, stop, interval_step))

  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      function(tp){ estimate_abundances(ccm, tibble(timepoint=tp, expt=reference_experiment))}, NA_real_),
      .x = timepoint)) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  # cell_type_assignments = colData(ccm@ccs@cds) %>%
  #   as.data.frame %>%
  #   dplyr::count(cluster, cell_type) %>%
  #   group_by(cluster) %>% slice_max(n) %>%
  #   dplyr::select(cell_group=cluster, cell_type)
  #
  # timepoint_pred_df = left_join(timepoint_pred_df, cell_type_assignments)
  # timepoint_pred_df = timepoint_pred_df %>% mutate(cell_group_label = paste(cell_type, " (", cell_group, ")", sep=""))

  return(timepoint_pred_df)
}

plot_contrast_wrapper <- function(ccm, t1, t2) {

  timepoint_pred_df = estimate_abundances_over_interval(ccm, t1, t2)

  plot_contrast(ccm, compare_abundances(ccm,
                                           timepoint_pred_df %>% filter(timepoint == t1),
                                           timepoint_pred_df %>% filter(timepoint == t2)),
                scale_shifts_by = "none",
                q_value_thresh = 0.01)

}

#t1 = 18
#t2 = 22
plot_contrast_wrapper(wt_ccm_wl, 18, 22)

# -----------------------------------------------------------------------------


# get_valid_origins <- function(ccm, start, stop, interval_step=2, percent_max_threshold = 0.1){
#   #timepoints = seq(start, stop, interval_step)
#
#   timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_step)
#   timepoint_pred_df = timepoint_pred_df %>% group_by(cell_group) %>% mutate(max_abundance = max(exp(log_abund)),
#                                                                             percent_max_abund = exp(log_abund) / max_abundance,
#                                                                             present_above_thresh = percent_max_abund > percent_max_threshold)
#   emergence_time = timepoint_pred_df %>% filter(present_above_thresh) %>% summarize(emergence_time = min(timepoint))
#
#   #cell_types_present = timepoint_pred_df %>% filter(percent_max_abund > overlap_theshold)
#   cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
#   valid_origins = tibble(id=cell_groups)
#
#   filter_origins <- function(ct, cell_types_present){
#     ct_times = cell_types_present %>% filter(cell_group == ct) %>% pull(timepoint) %>% unique()
#     possible_origins = cell_types_present %>% filter(timepoint %in% ct_times & cell_group != ct) %>% pull(cell_group) %>%unique
#     return (possible_origins)
#   }
#   valid_origins = valid_origins %>% mutate(possible_origins = purrr::map(.f = filter_origins,
#                                                                          .x = id,
#                                                                          cell_types_present=cell_types_present))
#   return(valid_origins)
# }
# debug(get_valid_origins)
# get_valid_origins(wt_ccm_wl, 18, 96)

get_extant_cell_types <- function(ccm,
                                  start,
                                  stop,
                                  interval_step = 2,
                                  log_abund_detection_thresh = -5,
                                  percent_max_threshold=NULL){

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_step)
  timepoint_pred_df = timepoint_pred_df %>% group_by(cell_group) %>% mutate(max_abundance = max(exp(log_abund)),
                                                                            percent_max_abund = exp(log_abund) / max_abundance,
                                                                            present_above_thresh = log_abund > log_abund_detection_thresh)
  if (is.null(percent_max_threshold) == FALSE){
    timepoint_pred_df = timepoint_pred_df %>%
      mutate(present_above_thresh = ifelse(present_above_thresh & percent_max_abund > percent_max_threshold, TRUE, FALSE))
  }

  extant_cell_type_df = timepoint_pred_df %>% select(timepoint, cell_group, log_abund, max_abundance, percent_max_abund, present_above_thresh)
  return(extant_cell_type_df)
}
#debug(get_extant_cell_types)
#get_extant_cell_types(wt_ccm_wl, 72, 96) %>% filter(cell_group == "21")

find_paths_to_origins <- function(ccm,
                                     cond_b_vs_a_tbl,
                                     abund_over_time_df,
                                     paga_graph_edges,
                                     q_value_threshold = 1.0,
                                     require_presence_at_all_timepoints = TRUE) {

  all_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl)
  #all_edges = paga_graph_edges
  paga_graph = igraph::graph_from_data_frame(paga_graph_edges, directed=FALSE)
  weighted_edges = hooke:::get_weighted_edges(ccm, all_edges)

  weighted_edges = weighted_edges %>% group_by(from, to) %>% mutate(adjacent_in_paga = igraph::are_adjacent(paga_graph, from, to)) %>% ungroup()
  weighted_edges = weighted_edges %>% filter(adjacent_in_paga)

  cell_types_exant_in_contrast = abund_over_time_df %>% filter(timepoint >= min(cond_b_vs_a_tbl$timepoint_x) &
                                                                timepoint <= max(cond_b_vs_a_tbl$timepoint_y))

  # These are the cell groups that are missing in at least one timepoint in the interval we're looking at:
  if (require_presence_at_all_timepoints){
    cell_groups_missing_in_range = cell_types_exant_in_contrast %>% filter(present_above_thresh == FALSE) %>% pull(cell_group) %>% unique
    cell_groups_missing_in_range = union(setdiff(row.names(ccm@ccs), cell_types_exant_in_contrast$cell_group), cell_groups_missing_in_range)
  }else{
    cell_groups_present_in_range = cell_types_exant_in_contrast %>% filter(present_above_thresh == TRUE) %>% pull(cell_group) %>% unique
    cell_groups_missing_in_range = setdiff(row.names(ccm@ccs), cell_groups_present_in_range)
  }

  weighted_edges = weighted_edges %>% filter(from %in% cell_groups_missing_in_range == FALSE &
                                               to %in% cell_groups_missing_in_range == FALSE)

  neg_rec_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(edge_type != "undirected" &
             to_delta_q_value < q_value_threshold &
             from_delta_q_value < q_value_threshold)

  assertthat::assert_that(
    tryCatch(expr = nrow(neg_rec_edges) != 0,
             error = function(e) FALSE),
    msg = "no significant negative reciprocal edges found")

  edge_path = neg_rec_edges %>%
    dplyr::mutate(shortest_path = purrr::map2(.f =
                                                purrr::possibly(hooke:::get_shortest_path, NA_real_),
                                              .x = from, .y = to,
                                              weighted_edges)) %>%
    select(origin=from, destination=to, shortest_path) %>%
    tidyr::unnest(shortest_path) %>%
    # select(-weight) %>%
    distinct()

  return(edge_path)
}
undebug(find_paths_to_origins)

#
# xxx_paths = find_paths_to_origins(wt_ccm_wl,
#                          compare_abundances(wt_ccm_wl,
#                                             estimate_abundances(wt_ccm_wl, tibble(expt="GAP16", timepoint=40)),
#                                             estimate_abundances(wt_ccm_wl, tibble(expt="GAP16", timepoint=72))),
#                          get_extant_cell_types(wt_ccm_wl, 18, 96, log_abund_detection_thresh = -5, percent_max_threshold=0.01),
#                          paga_edges,
#                          require_presence_at_all_timepoints = TRUE)

find_origins <- function(ccm,
                            paga_graph_edges,
                            q_val=0.01,
                            start = NULL,
                            stop = NULL,
                            interval_step = 2,
                            min_interval = 4,
                            max_interval = 24,
                            percent_max_threshold=NULL,
                            log_abund_detection_thresh=-5,
                            require_presence_at_all_timepoints=TRUE) {

  if (is.null(start)){
    start = min(colData(ccm@ccs)$timepoint)
  }
  if (is.null(stop)){
    stop = max(colData(ccm@ccs)$timepoint)
  }

  timepoints = seq(start, stop, interval_step)

  extant_cell_type_df = get_extant_cell_types(ccm,
                                              start,
                                              stop,
                                              percent_max_threshold=percent_max_threshold,
                                              log_abund_detection_thresh=log_abund_detection_thresh)

  timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_step)

  select_timepoints <- function(timepoint_pred_df, t1, t2)  {
    cond_x = timepoint_pred_df %>% filter(timepoint == t1)
    cond_y = timepoint_pred_df %>% filter(timepoint == t2)
    return(compare_abundances(ccm, cond_x, cond_y))
  }

  times = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2 & (t2-t1) >= min_interval & (t2-t1) <= max_interval) %>%
    mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(neg_rec_edges = purrr::map(.f = hooke:::get_neg_dir_edges,
                                      .x = comp_abund,
                                      ccm = ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(pos_rec_edges = purrr::map(.f = hooke:::get_positive_edges,
                                      .x = comp_abund,
                                      ccm = ccm,
                                      q_value_threshold = 1.0)) %>%
    mutate(path = purrr::map(.f = purrr::possibly(find_paths_to_origins, NA_real_),
                             .x = comp_abund,
                             ccm = ccm,
                             paga_graph_edges=paga_graph_edges,
                             abund_over_time_df = extant_cell_type_df,
                             require_presence_at_all_timepoints=require_presence_at_all_timepoints,
                             q_value_threshold = 1.0))
  return (times)
}
#undebug(find_origins)

wt_possible_origins = find_origins(wt_ccm_wl,
                         paga_graph_edges=initial_pcor_graph(wt_ccs),
                         start=18, stop=96,
                         min_interval = 2,
                         log_abund_detection_thresh=-5,
                         percent_max_threshold=0.01,
                         require_presence_at_all_timepoints=TRUE)

# xxx_paths = wt_tcs %>% select(t1, t2, path) %>%
#   filter(!is.na(path)) %>%
#   tidyr::unnest(path) %>%
#   select(t1, t2, origin, destination, from, to, umap_dist=weight)
# xxx_neg_rec_paths = wt_tcs %>%
#   select(neg_rec_edges) %>%
#   tidyr::unnest(neg_rec_edges) %>%
#   select(origin=from, destination=to, origin_pcor=pcor) %>% distinct() %>%
#   mutate(origin_pcor = -origin_pcor)
# xxx_paths = left_join(xxx_paths, xxx_neg_rec_paths)
#
# xxx_paths = xxx_paths %>% mutate(from = ifelse(is.na(from), origin, from),
#                                            to = ifelse(is.na(to), destination, to))
# xxx_paths = xxx_paths %>% filter(is.na(umap_dist) == FALSE)
# xxx_paths = xxx_paths %>% select(origin, destination, from, to, origin_pcor) %>% distinct()
# pcor_path_matrix = xxx_paths %>% group_by(from, to) %>% summarize(max_origin_pcor = max(origin_pcor))
# #pcor_path_graph = igraph::graph_from_data_frame(pcor_path_matrix)

select_origins <- function(ccm, possible_origins, fraction_origin_pcor_thresh=0.5, overlap_theshold=0.25, allow_discontinuous_paths = FALSE){

  neg_rec_edges = possible_origins %>%
    select(neg_rec_edges) %>%
    tidyr::unnest(neg_rec_edges) %>%
    select(origin=from, destination=to, origin_pcor=pcor) %>% distinct() %>%
    mutate(origin_pcor = -origin_pcor)

  weighted_paths = possible_origins %>% select(t1, t2, path) %>%
    filter(!is.na(path)) %>%
    tidyr::unnest(path) %>%
    select(t1, t2, origin, destination, from, to, umap_dist=weight)


  #valid_origins = get_valid_origins(ccm, overlap_theshold)
  #weighted_paths = left_join(weighted_paths, valid_origins, by=c("destination"="id"))
  #weighted_paths = weighted_paths %>% filter(origin %in% unlist(possible_origins)) %>% select(-possible_origins) %>% distinct()
  weighted_paths = left_join(weighted_paths, neg_rec_edges)

  weighted_paths = weighted_paths %>% mutate(from = ifelse(is.na(from), origin, from),
                                             to = ifelse(is.na(to), destination, to))
  if (allow_discontinuous_paths == FALSE){
    weighted_paths = weighted_paths %>% filter(is.na(umap_dist) == FALSE)
  }

  # origin_summary = weighted_paths %>% group_by(destination) %>% summarize(num_origins = length(unique(origin)),
  #                                                                         num_indirect_origins = length(unique(origin[umap_dist != -1])))

  weighted_paths = weighted_paths %>% group_by(destination) %>% mutate(max_origin_pcor = max(origin_pcor),
                                                                       fraction_max_origin_pcor = origin_pcor / max_origin_pcor)

  weighted_paths = weighted_paths %>% filter(fraction_max_origin_pcor > fraction_origin_pcor_thresh)

  origin_df = weighted_paths %>%
    select(from, to, fraction_max_origin_pcor) %>%
    group_by(from, to) %>%
    summarize(fraction_max_origin_pcor = mean(fraction_max_origin_pcor), n=1) %>%
    group_by(to) %>%
    mutate(scaled_weight=fraction_max_origin_pcor / max(fraction_max_origin_pcor))
    #tally() %>%
    #mutate(scaled_weight = abs(n) / max(abs(n)))

  return(origin_df)
}

debug(select_origins)

wt_origins = select_origins(wt_ccm_wl, wt_possible_origins, fraction_origin_pcor_thresh=0.5)

hooke:::plot_path(wt_ccm_wl, path_df = wt_origins, edge_size=0.25)

plot_state_transition_graph <- function(ccm,
                                        edges,
                                        color_nodes_by=NULL,
                                        label_nodes_by=NULL,
                                        group_nodes_by=NULL,
                                        layer_nodes_by=NULL,
                                        arrow.gap=0.03,
                                        arrow_unit = 2,
                                        bar_unit = .075,
                                        node_size = 2,
                                        num_layers=10,
                                        unlabeled_groups = c("Unknown"),
                                        hide_unlinked_nodes=TRUE,
                                        label_font_size=6,
                                        label_conn_linetype="dotted"){

  edges = hooke:::distance_to_root(edges)

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs@cds)[,c(color_nodes_by,
                                            label_nodes_by,
                                            group_nodes_by,
                                            layer_nodes_by)] %>%
    as.data.frame
  cell_group_metadata$cell_group = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  if (is.null(color_nodes_by) == FALSE){
    color_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(color_by_metadata) = c("cell_group", "color_nodes_by")
    node_metadata = left_join(node_metadata, color_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE){
    group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
    node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(layer_nodes_by) == FALSE){
    layer_by_metadata = cell_group_metadata[,c("cell_group", layer_nodes_by)] %>%
      as.data.frame
    if (is.numeric(cell_group_metadata[,c(layer_nodes_by)])){
      layer_by_metadata = layer_by_metadata %>%
        group_by(cell_group) %>%
        summarize(mean_layer_var = mean(!!sym(layer_nodes_by), na.rm=TRUE)) %>%
        mutate(layer = ntile(desc(mean_layer_var),num_layers)) %>% dplyr::select(-mean_layer_var)
    }else{
      layer_by_metadata = layer_by_metadata %>%
        count(cell_group, !!sym(layer_nodes_by)) %>%
        group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    }
    colnames(layer_by_metadata) = c("cell_group", "layer_nodes_by")
    node_metadata = left_join(node_metadata, layer_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE){
    label_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(label_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    node_metadata = left_join(node_metadata, label_by_metadata, by=c("id"="cell_group"))
  }else{
    node_metadata$label_nodes_by = node_metadata$id
  }
  node_metadata = node_metadata %>% distinct() %>% as.data.frame
  row.names(node_metadata) = node_metadata$id
  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  G = edges %>% select(from, to, scaled_weight) %>% distinct()  %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  # run sugiyama layout
  layers = NULL
  if (is.null(layer_nodes_by) == FALSE) {
    layers=igraph::V(G)$layer_nodes_by
  }
  lay1 <- igraph::layout_with_sugiyama(G, layers=layers, maxiter=1000)

  g = ggnetwork::ggnetwork(G, layout = lay1$layout, arrow.gap = arrow.gap)

  # add level information
  #g = g %>% left_join(level_df %>% rownames_to_column("id"), by = c("vertex.names"="id"))
  #g = g %>% left_join(regulator_score_df, by = c("vertex.names" = "gene_id") )


  p <- ggplot(mapping = aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
    # draw activator edges
    ggnetwork::geom_edges(data = g,
                          arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))
  if (is.null(group_nodes_by) == FALSE){
    p = p + ggforce::geom_mark_rect(aes(fill = group_nodes_by,
                                        label=group_nodes_by,
                                        filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                    size=0,
                                    label.fontsize=label_font_size,
                                    con.linetype=label_conn_linetype,
                                    data=g)
  }

  if (is.null(color_nodes_by) == FALSE) {

    # if numerical
    if (is.numeric(g[[color_nodes_by]])) {
      p = p + ggnetwork::geom_nodelabel(data = g,
                                        aes(fill = as.factor(color_nodes_by),
                                            label = label_nodes_by),
                                        size = node_size) +
        labs(fill = color_nodes_by) +
        scale_fill_gradient2(low = "darkblue", mid = "white", high="red4")
    }
    else {
      # if categorical
      p = p + ggnetwork::geom_nodelabel(data = g,
                                        aes(fill = color_nodes_by,
                                            label = label_nodes_by),
                                        size = node_size) +
        labs(fill = color_nodes_by)

    }

  } else {
    p = p + ggnetwork::geom_nodelabel(data = g,
                                      aes(label = label_nodes_by),
                                      size = node_size)
  }

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position="none")
  return(p)
}
undebug(plot_state_transition_graph)
plot_state_transition_graph(wt_ccm_wl, wt_origins, color_nodes_by = "timepoint", group_nodes_by="cell_type_broad")


plot_state_transition_graph(wt_ccm_wl, pos_edge_paths, color_nodes_by = "timepoint", group_nodes_by="cell_type")

