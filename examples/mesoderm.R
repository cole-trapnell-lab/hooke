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


ggplot(aes(x=expt, y=mean_nn_time),
       data=colData(meso_cds) %>% as.data.frame) +
  geom_boxplot() + facet_wrap(~timepoint, scale="free") +
  monocle3:::monocle_theme_opts()

ggplot(aes(x=as.numeric(timepoint), y=mean_nn_time, color=temp),
       data=colData(meso_cds) %>% as.data.frame) + geom_jitter() +
  facet_wrap(~expt) + geom_abline(color="grey") + geom_smooth(method="lm", color="black") +
  ylim(c(16,100)) + monocle3:::monocle_theme_opts()

staging_df = colData(meso_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))

staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
staging_df$predicted_timepoint = predict(staging_model)

colData(meso_cds)$adjusted_timepoint = staging_df$predicted_timepoint

# JUST FOR DEBUGGING:
colData(meso_cds)$timepoint = NULL

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


# Use custom breaks because the sampling over the time range is so uneven
#wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_var="adjusted_timepoint", interval_start=18, interval_stop=96, num_breaks=6)
wt_main_model_formula_str = "~ splines::ns( adjusted_timepoint , knots= c(30,45,60), Boundary.knots=c(20,92) )"

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 #main_model_formula_str = "~expt",
                                 main_model_formula_str = wt_main_model_formula_str,
                                 #nuisance_model_formula_str = "~1",
                                 #nuisance_model_formula_str = "~expt",
                                 whitelist = initial_pcor_graph(wt_ccs))

wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=5)

xxx_interval_abundances = estimate_abundances_over_interval(wt_ccm_wl, 18, 96, interval_col="adjusted_timepoint")
qplot(adjusted_timepoint, log_abund, color = log_abund - 3*log_abund_se > 0, geom="point", data=xxx_interval_abundances) +
  #geom_ribbon(aes(ymin=log_abund - 2*log_abund_se, ymax=log_abund + 2*log_abund_se), alpha=0.25) +
  facet_wrap(~cell_group)

plot_contrast_wrapper <- function(ccm, t1, t2, q_val=0.01) {

  timepoint_pred_df = estimate_abundances_over_interval(ccm, t1, t2, interval_col="adjusted_timepoint", expt="GAP16")

  plot_contrast(ccm, compare_abundances(ccm,
                                           timepoint_pred_df %>% filter(adjusted_timepoint == t1),
                                           timepoint_pred_df %>% filter(adjusted_timepoint == t2)),
                scale_shifts_by = "none",
                q_value_thresh = q_val)

}

#t1 = 18
#t2 = 22
#debug(plot_contrast_wrapper)
plot_contrast_wrapper(wt_ccm_wl, 18, 22)

# -----------------------------------------------------------------------------


xxx_extant_cell_type_df = get_extant_cell_types(wt_ccm_wl,
                                                18,
                                                96,
                                                interval_col="adjusted_timepoint",
                                                percent_max_threshold=0.00,
                                                log_abund_detection_thresh=-2)

qplot(adjusted_timepoint, log_abund, color = present_above_thresh, geom="point", data=xxx_extant_cell_type_df) +
  #geom_ribbon(aes(ymin=log_abund - 2*log_abund_se, ymax=log_abund + 2*log_abund_se), alpha=0.25) +
  facet_wrap(~cell_group)


# This function chooses the sparsity parameter by finding the largest value that provides at least one
# possible origin for cell types that emerge in the time interval
get_emergent_cell_types <- function(ccm, start, stop, interval_col="timepoint", ...){
  extant_cell_type_df = get_extant_cell_types(ccm, start, stop, interval_col, ...)
  emergent_cell_types = extant_cell_type_df %>% ungroup() %>%
    group_by(cell_group) %>%
    filter (present_above_thresh) %>%
    arrange(cell_group) %>%
    #filter(cell_group == "25") %>%
    slice_min(!!sym(interval_col)) %>% filter(!!sym(interval_col) > start) %>%
    pull(cell_group) %>% unique()
  return (emergent_cell_types)
}
debug(get_emergent_cell_types)
emergent_cell_types = get_emergent_cell_types(wt_ccm_wl, start=18, stop=96, interval_col="adjusted_timepoint", percent_max_threshold=0.25)



wt_possible_origins = find_origins(wt_ccm_wl,
                         start=18, stop=96,
                         interval_col="adjusted_timepoint",
                         min_interval = 2,
                         log_abund_detection_thresh=-2,
                         percent_max_threshold=0.01,
                         require_presence_at_all_timepoints=TRUE,
                         initial_origin_policy="all-origins")

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


wt_origins = select_origins(wt_ccm_wl, wt_possible_origins, selection_policy = "acceptable-origins")
hooke:::plot_path(wt_ccm_wl, path_df = wt_origins, edge_size=0.25)

origin_edge_graph = wt_origins %>% select(from, to, origin_pcor) %>% igraph::graph_from_data_frame()


A = origin_edge_graph %>% igraph::as_adjacency_matrix(attr="origin_pcor")
B = origin_edge_graph %>% igraph::as_adjacency_matrix(attr="origin_pcor")

A[lower.tri(A)] = 0
B[upper.tri(B)] = 0

A[A < t(B)] = 0
B[B < t(A)] = 0

final_pcor_adj_mat = A + B
final_pcor_graph_df = final_pcor_adj_mat %>% igraph::graph_from_adjacency_matrix(weighted="origin_pcor") %>% igraph::as_data_frame()
wt_origins_filtered = inner_join(final_pcor_graph_df, wt_origins) %>% as_tibble()


hooke:::plot_path(wt_ccm_wl, path_df = wt_origins_filtered, edge_size=0.25)

undebug(plot_state_transition_graph)
plot_state_transition_graph(wt_ccm_wl, wt_origins, color_nodes_by = "cell_type_sub", group_nodes_by="cell_type_broad")


plot_state_transition_graph(wt_ccm_wl, pos_edge_paths, color_nodes_by = "timepoint", group_nodes_by="cell_type")

