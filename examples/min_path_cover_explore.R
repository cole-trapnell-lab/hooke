library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)

setwd("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/")
devtools::load_all("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/")

# let's test this on mesoderm first
meso_cds = readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")

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

# 1) Build a directed graph on the nodes by first computing the shortest paths between negative pcor
# reciprocal pairs through the existing path finding graph and then taking the union of all of these paths


# First, let's figure out when each cell type is present and
# which ones emerge over the course of the caller's time interval
if (is.null(start)){
  start = min(colData(ccm@ccs)[,interval_col])
}
if (is.null(stop)){
  stop = max(colData(ccm@ccs)[,interval_col])
}

timepoints = seq(start, stop, interval_step)

timepoint_pred_df = estimate_abundances_over_interval(ccm, start, stop, interval_col=interval_col, interval_step)

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
  dplyr::filter(pcor < 0) %>% # do we just want negative again?
  dplyr::filter((from_delta_log_abund > 0 & to_delta_log_abund < 0) |
                  (to_delta_log_abund > 0 & from_delta_log_abund < 0))

edge_union = relevant_comparisons %>% select(from, to) %>% distinct()

# currently just take a subset for testing
# edge_union = edge_union %>% head(25)

plot_path(meso_wt_ccm_wl, edge_union, edge_size = 1)


# does the edge go backwards in time?

# find cells that fall along a path

paths_to_origins = find_paths_to_origins(selected_origins, timeseries_pathfinding_graph)

interval_start = min(neg_rec_edges_to_destinations[,paste("to", interval_col, "x", sep="_")])
interval_stop =  max(neg_rec_edges_to_destinations[,paste("to", interval_col, "y", sep="_")])



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
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}


test_cycle = data.frame("from" = c(12, 31, 15),
                        "to" = c(31, 15, 12))

g = igraph::graph_from_data_frame(test_cycle)
plot(g)
find_cycles(g)


#' score a path based on fitting a linear model of time ~ geodeisic distance
#' @param ccs
#' @param path_df
cells_along_path <- function(ccs, path_df) {
  cds = ccs@cds
  vertices = union(path_df$to, path_df$from) %>% unique()

  path_df = distance_to_root(path_df)
  cds_along_path = cds[, colData(cds)[[ccs@info$cell_group]] %in% vertices]
  colData(cds_along_path)$cell_group = colData(cds_along_path)[[ccs@info$cell_group]]

  cells_along_path_df = colData(cds_along_path) %>% as.data.frame %>% select(cell, cell_group, adjusted_timepoint) %>% as_tibble() %>%
    left_join(path_df %>% select(-from), by = c("cell_group" = "to"))

  cells_along_path_df$distance_from_root = tidyr::replace_na(cells_along_path_df$distance_from_root, 0)
  cells_along_path_df$y = cells_along_path_df[[interval_col]]

  path_model = lm(y ~ distance_from_root, data=cells_along_path_df)

  return(coef(path_model)[["distance_from_root"]])

}



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


