library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)

#setwd("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/")
#devtools::load_all("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/")


setwd("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/mesoderm")

meso_cds = readRDS("final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")


# let's test this on mesoderm first
#meso_cds = readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")

colData(meso_cds)$cluster = as.character(clusters(meso_cds))
meso_cds = meso_cds[, is.na(colData(meso_cds)$embryo) == FALSE]
colData(meso_cds)$cell_type = colData(meso_cds)$cell_type_sub

staging_df = colData(meso_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))

staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
staging_df$predicted_timepoint = predict(staging_model)

colData(meso_cds)$adjusted_timepoint = staging_df$predicted_timepoint

meso_wt_ccs = new_cell_count_set(meso_cds,
                                 sample_group = "embryo",
                                 cell_group = "cluster")

meso_wt_main_model_formula_str = "~ splines::ns( adjusted_timepoint , knots= c(24,48,60), Boundary.knots=c(20,92) )"

meso_wt_ccm_wl = new_cell_count_model(meso_wt_ccs,
                                      #main_model_formula_str = "~expt",
                                      main_model_formula_str = meso_wt_main_model_formula_str,
                                      #nuisance_model_formula_str = "~1",
                                      #nuisance_model_formula_str = "~expt",
                                      whitelist = initial_pcor_graph(meso_wt_ccs))
meso_wt_ccm_wl = select_model(meso_wt_ccm_wl, criterion = "EBIC", sparsity_factor= 0.01)


ccs = meso_wt_ccs
ccm = meso_wt_ccm_wl
interval_col = "adjusted_timepoint"
interval_step = 2
min_interval = 4
max_interval = 24
start_time = NULL
stop_time = NULL
log_abund_detection_thresh = -2
percent_max_threshold = 0
q_value_thresh = 0.01

# 1) Build a directed graph on the nodes by first computing the shortest paths between negative pcor
# reciprocal pairs through the existing path finding graph and then taking the union of all of these paths


# First, let's figure out when each cell type is present and
# which ones emerge over the course of the caller's time interval
if (is.null(start_time)){
  start_time = min(colData(ccm@ccs)[,interval_col])
}
if (is.null(stop_time)){
  stop_time = max(colData(ccm@ccs)[,interval_col])
}

timepoints = seq(start_time, stop_time, interval_step)

timepoint_pred_df = estimate_abundances_over_interval(ccm, start_time, stop_time, interval_col=interval_col, interval_step)

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
  dplyr::filter(from_delta_q_value < q_value_thresh & to_delta_q_value < q_value_thresh)

edge_union = relevant_comparisons %>% select(from, to) %>% distinct()

# currently just take a subset for testing
# edge_union = edge_union %>% head(25)

hooke:::plot_path(meso_wt_ccm_wl, edge_union, edge_size = 0.25)


# does the edge go backwards in time?

# find cells that fall along a path

extant_cell_type_df = hooke:::get_extant_cell_types(ccm,
                                            start_time,
                                            stop_time,
                                            interval_col=interval_col,
                                            percent_max_threshold=percent_max_threshold,
                                            log_abund_detection_thresh=log_abund_detection_thresh)

emergent_cell_types = extant_cell_type_df %>%
  select(cell_group, emerges_at=longest_contig_start) %>%
  distinct() %>%
  filter(emerges_at > min(emerges_at)) %>% pull(cell_group)

timeseries_pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)

#paths_to_origins = find_paths_to_origins(selected_origins, timeseries_pathfinding_graph)

#interval_start = min(neg_rec_edges_to_destinations[,paste("to", interval_col, "x", sep="_")])
#interval_stop =  max(neg_rec_edges_to_destinations[,paste("to", interval_col, "y", sep="_")])



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


test_cycle = data.frame("from" = c(12, 15),
                        "to" = c(15, 12))

g = igraph::graph_from_data_frame(test_cycle)
#g = igraph::delete_edges(g, "15|12")
plot(g)
length(find_cycles(g))


test_cycle = data.frame("from" = c(12, 31, 15),
                        "to" = c(31, 15, 12))

g = igraph::graph_from_data_frame(test_cycle)
#g = igraph::delete_edges(g, "15|12")
plot(g)
length(find_cycles(g))


#' score a path based on fitting a linear model of time ~ geodeisic distance
#' @param ccs
#' @param path_df
cells_along_path <- function(path_df, ccs) {
  cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()
  path_df = path_df %>% mutate(geodesic_dist = cumsum(weight))
  #path_df = hooke:::distance_to_root(path_df)
  cds_along_path = cds[, colData(cds)[[ccs@info$cell_group]] %in% vertices]
  colData(cds_along_path)$cell_group = colData(cds_along_path)[[ccs@info$cell_group]]

  cells_along_path_df = colData(cds_along_path) %>% as.data.frame %>% select(cell, cell_group, adjusted_timepoint) %>% as_tibble() %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  cells_along_path_df$geodesic_dist = tidyr::replace_na(cells_along_path_df$geodesic_dist, 0)

  cells_along_path_df$y = cells_along_path_df[[interval_col]]

  path_model = lm(y ~ geodesic_dist, data=cells_along_path_df)

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

paths_for_relevant_edges = edge_union %>%
  #head (5) %>%
  mutate(path = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                            .x = from, .y = to,
                            timeseries_pathfinding_graph))


paths_for_relevant_edges = paths_for_relevant_edges %>%
  mutate(time_vs_distance_model_stats = purrr::map(.f = purrr::possibly(cells_along_path, NULL),
                            .x = path,
                            ccs)) %>%
  tidyr::unnest(time_vs_distance_model_stats)
#cells_along_path(ccs, paths_for_relevant_edges$path[[1]] %>% dplyr::select(from, to))


paths_for_relevant_edges = paths_for_relevant_edges %>% mutate(dist_model_score = dist_effect * dist_model_adj_rsq)

selected_paths = paths_for_relevant_edges %>% filter (dist_effect > 0 & dist_effect_pval < 0.01 & dist_model_adj_rsq > 0.0) %>%
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
    G <<- G_prime
  }else{
    # Debug:
    print ("skipping:")
    print (selected_paths[i,] %>% select(-path))
  }
}
igraph::E(G)$weight = igraph::E(timeseries_pathfinding_graph)[igraph::E(G)]$weight

hooke:::plot_path(meso_wt_ccm_wl, G %>% igraph::as_data_frame(), edge_size = 0.25)


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
pcor_mat_summ = summary(pcor_mat)
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

paths_from_chains = matching %>% as_tibble() %>%
  mutate(chain_leg = purrr::map2(.f = purrr::possibly(hooke:::get_shortest_path, NULL),
                            .x = from, .y = to,
                            node_dag))

covered_graph = paths_from_chains %>% select(chain_leg) %>%
  tidyr::unnest() %>%
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
  tidyr::unnest() %>%
  igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(node_dag)$name)) %>%
  igraph::union(covered_graph)

# covered_graph = paths_to_chain_tails %>% select(chain_leg) %>%
#   tidyr::unnest() %>%
#   igraph::graph_from_data_frame(vertices=data.frame(id=igraph::V(node_dag)$name)) %>%
#   igraph::union(covered_graph)

covered_graph = igraph::simplify(covered_graph)

covered_graph = igraph::delete_vertices(covered_graph, c("source", "sink"))

hooke:::plot_path(meso_wt_ccm_wl, covered_graph %>% igraph::as_data_frame(), edge_size = 0.25)

plot_state_transition_graph(meso_wt_ccm_wl, covered_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_sub", group_nodes_by="cell_type_broad")


#covered_graph = igraph::graph_from_data_frame(matching)


hooke:::plot_path(meso_wt_ccm_wl, paths_from_chains %>% select(path) %>% tidyr::unnest() %>% select(from, to), edge_size = 0.25)

# 3) Compute a min path cover on the resulting DAG


G1 <- igraph::graph(c(1, 2, 1, 3, 1, 4, 3, 4, 3, 5,
                      5, 6, 6, 7, 7, 8, 8, 9, 3, 8, 5, 8), directed = FALSE)

mbm = maxmatching::maxmatching(G1, weighted = F)
mbm$matching

G2 <- igraph::graph(c(1, 5, 1, 6, 1, 7, 2, 5, 2, 8,
                      3, 6, 3, 7, 3, 8, 4, 6, 4, 7, 4, 8), directed = FALSE)
plot(G2)
mbm = maxmatching::maxmatching(G2, weighted = FALSE)
mbm$matching


G3 <- igraph::graph(c(1, 5, 1, 6, 1, 7, 2, 5, 2, 8,
                      3, 6, 3, 7, 3, 8, 4, 6, 4, 7, 4, 8), directed = FALSE)
igraph::E(G3)$weight <- runif(length(igraph::E(G3)), 0, 1)
plot(G3)
mbm = maxmatching::maxmatching(G3, weighted = TRUE)
data.frame(mbm$matching)

# from https://github.com/cole-trapnell-lab/cicero/issues/15
find_best_ccan_match <- function(ccans1, ccans2) {
  names(ccans1) <- c("Peak", "CCAN1")
  names(ccans2) <- c("Peak", "CCAN2")
  comps_mem <- merge(ccans1, ccans2, all=TRUE)
  comps_mem$value <- 1

  comps_mem$CCAN1 <- paste0("C1_", comps_mem$CCAN1)
  comps_mem$CCAN2 <- paste0("C2_", comps_mem$CCAN2)

  comps_mem_temp <- subset(comps_mem, CCAN1 != "C1_NA" | CCAN2 != "C2_NA")
  comps_mem_dc <- reshape2::dcast(CCAN1 ~ CCAN2, value.var = "value", data = comps_mem_temp)

  comps_mem_dc <- subset(comps_mem_dc, CCAN1 != "C1_NA")
  comps_mem_dc <- comps_mem_dc[,1:(ncol(comps_mem_dc) - 1)]
  row.names(comps_mem_dc) <- comps_mem_dc$CCAN1
  comps_mem_dc$CCAN1 <- NULL

  comps_mem_dc_melt <- comps_mem_dc
  comps_mem_dc_melt$CCAN1 <- row.names(comps_mem_dc_melt)
  comps_mem_dc_melt <- melt(comps_mem_dc_melt, id.vars = "CCAN1")
  names(comps_mem_dc_melt) <- c("from", "to", "weight")
  comps_mem_dc_melt <- subset(comps_mem_dc_melt, from != to)
  grph <- igraph::graph_from_data_frame(subset(comps_mem_dc_melt, weight > 0), directed = FALSE)
  igraph::V(grph)$type <- grepl("orig", igraph::V(grph))
  mbm <- maxmatching::maxmatching(grph)

  matching <- data.frame(CCAN1 = mbm$matching, CCAN2 = names(mbm$matching))
  matching <- subset(matching, grepl("C1", CCAN1) & grepl("C2", CCAN2))
  matching$CCAN1 <- gsub("C1_", "", matching$CCAN1)
  matching$CCAN2 <- gsub("C2_", "", matching$CCAN2)
  return(matching)
}


