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

  # FIXME: this strategy doesn't really work well. The problem is that it creates lots of connections between
  # PAGA partitions, which causes problems later on. We really need a way to link partitions with only a very
  # few good edges.
  if (allow_links_between_components){
    # Allow links between cell groups in different partitions, but only if they have negative partial correlation
    # weighted_edges = weighted_edges %>% filter(adjacent_in_paga | (!adjacent_in_paga & weight < 0))
    weighted_edges = weighted_edges %>% filter(adjacent_in_paga | (same_partition == FALSE & pcor < 0))
  }else{
    weighted_edges = weighted_edges %>% filter(adjacent_in_paga)
  }


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
cells_along_path <- function(path_df, ccs, interval_col="timepoint") {


  #cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% arrange(distance_from_root) %>% mutate(geodesic_dist = cumsum(weight))

  cells_along_path_df = counts(ccs)[vertices,] %>%
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

  cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$y = cells_along_path_df[[interval_col]]

  path_model = lm(y ~ geodesic_dist, data=cells_along_path_df, weights=cells_along_path_df$num_cells)

  path_model_tidied = broom::tidy(path_model)
  path_model_glanced = broom::glance(path_model)
  pm_stats = tibble(dist_effect = unlist(path_model_tidied[2, "estimate"]),
                    dist_effect_pval = unlist(path_model_tidied[2, "p.value"]),
                    dist_model_adj_rsq = unlist(path_model_glanced[1, "adj.r.squared"]),
                    dist_model_ncells = unlist(path_model_glanced[1, "nobs"]))
  return(pm_stats)
  #return(coef(path_model)[["distance_from_root"]])
}
#debug(cells_along_path)



#' Identify the possible origins for each destination
#'
build_transition_dag <- function(ccm,
                                 extant_cell_type_df,
                                 timeseries_pathfinding_graph,
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


  paths_for_relevant_edges = edge_union %>%
    mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NA_character_),
                              .x = from, .y = to,
                              timeseries_pathfinding_graph))

  paths_for_relevant_edges = paths_for_relevant_edges %>%
    filter(!is.na(path)) %>%
    mutate(time_vs_distance_model_stats = purrr::map(.f = purrr::possibly(cells_along_path, NA_character_),
                                                     .x = path,
                                                     ccs=ccm@ccs,
                                                     interval_col=interval_col)) %>%
    tidyr::unnest(time_vs_distance_model_stats)
  #cells_along_path(ccs, paths_for_relevant_edges$path[[1]] %>% dplyr::select(from, to))


  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(dist_model_score = dist_effect * dist_model_adj_rsq)
  paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(dist_effect_q_val = p.adjust(dist_effect_pval, method="BH"))

  selected_paths = paths_for_relevant_edges %>% filter (dist_effect > 0 & dist_effect_q_val < q_val & dist_model_adj_rsq > min_dist_vs_time_r_sq) %>%
    group_by(to) %>%
    #slice_max(dist_model_score, n=3) %>%
    ungroup() %>% arrange(desc(dist_model_score))

  # paths_for_relevant_edges = paths_for_relevant_edges %>%
  #   #filter (dist_effect > 0 & dist_effect_pval < 0.01 & to %in% emergent_cell_types) %>%
  #   filter (dist_effect > 0 & dist_effect_pval < 0.01) %>%
  #   arrange(desc(dist_model_adj_rsq))


  #G = igraph::make_empty_graph()
  #G = igraph::add_vertices(G, length(igraph::V(timeseries_pathfinding_graph)), attr=list("name"=igraph::V(timeseries_pathfinding_graph)$name))

  G = timeseries_pathfinding_graph
  G = igraph::delete_edges(G, igraph::E(G))

  for (i in (1:nrow(selected_paths))){
    next_path = selected_paths$path[[i]] %>% select(from, to)
    next_path_graph = next_path %>% igraph::graph_from_data_frame(directed=TRUE)
    G_prime = igraph::union(G, next_path_graph)
    if (length(find_cycles(G_prime)) == 0){
      G = G_prime
    }else{
      # Debug:
      #print ("skipping:")
      #print (selected_paths[i,] %>% select(-path))
    }
  }
  igraph::E(G)$weight = igraph::E(timeseries_pathfinding_graph)[igraph::E(G)]$weight

  return (G)
}

compute_min_path_cover <- function(ccm, G){

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
  igraph::E(node_dag)$weight = tidyr::replace_na(igraph::E(node_dag)$weight, 0)

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

  timeseries_pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)


  G = build_transition_dag(ccm,
                           extant_cell_type_df,
                           timeseries_pathfinding_graph,
                           q_val,
                           start_time,
                           stop_time,
                           interval_col,
                           interval_step,
                           min_interval,
                           max_interval,
                           min_dist_vs_time_r_sq,
                           ...)

  covered_G = compute_min_path_cover(ccm, G)

  igraph::V(covered_G)$cell_group = ccm@ccs@info$cell_group

  return(covered_G)
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
    print(paste0("try more n starting at ", df[n,]$sparsity))
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
