library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)


# setwd("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/mesoderm")
setwd("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/examples/")
devtools::load_all("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model")

# full meso cds
meso_cds = readRDS("~/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/mesoderm_all_cells_projected_940k_anno_cds.RDS")


colData(meso_cds)$sample = NULL
# save this clustering ?
# meso_cds = cluster_cells(meso_cds, resolution=1e-4, random_seed=42)

plot_cells(meso_cds)

# meso_cds = readRDS("final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")
# mesocds = readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")

meso_wt_cds = meso_cds[, grepl("ctrl", colData(meso_cds)$gene_target)]

meso_wt_cds = cluster_cells(meso_wt_cds, resolution=1e-4, random_seed=42)
colData(meso_wt_cds)$cluster = as.character(clusters(meso_wt_cds))
meso_wt_cds = meso_wt_cds[, is.na(colData(meso_wt_cds)$embryo) == FALSE]
colData(meso_wt_cds)$cell_type = colData(meso_wt_cds)$cell_type_sub


ggplot(aes(x=expt, y=mean_nn_time),
       data=colData(meso_cds) %>% as.data.frame) +
  geom_boxplot() + facet_wrap(~timepoint, scale="free") +
  monocle3:::monocle_theme_opts()

ggplot(aes(x=as.numeric(timepoint), y=mean_nn_time, color=temp),
       data=colData(meso_cds) %>% as.data.frame) + geom_jitter() +
  facet_wrap(~expt) + geom_abline(color="grey") + geom_smooth(method="lm", color="black") +
  ylim(c(16,100)) + monocle3:::monocle_theme_opts()

staging_df = colData(meso_wt_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))

staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
staging_df$predicted_timepoint = predict(staging_model)

colData(meso_wt_cds)$adjusted_timepoint = staging_df$predicted_timepoint

# JUST FOR DEBUGGING:
colData(meso_cds)$timepoint = NULL

wt_ccs = new_cell_count_set(meso_wt_cds,
                            sample_group = "embryo",
                            cell_group = "cluster")


# wt_start = 18
# wt_stop = 96
# num_time_breaks = 5
# time_breakpoints = seq(wt_start, wt_stop, length.out=num_time_breaks)
# time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
# wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")


# Use custom breaks because the sampling over the time range is so uneven
#wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_var="adjusted_timepoint", interval_start=18, interval_stop=96, num_breaks=6)
wt_main_model_formula_str = "~ splines::ns( adjusted_timepoint , knots= c(24,48,60), Boundary.knots=c(20,92) )"

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 #main_model_formula_str = "~expt",
                                 main_model_formula_str = wt_main_model_formula_str,
                                 #nuisance_model_formula_str = "~1",
                                 #nuisance_model_formula_str = "~expt",
                                 whitelist = initial_pcor_graph(wt_ccs))

wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=1)

# xxx_interval_abundances = estimate_abundances_over_interval(wt_ccm_wl, 18, 96, interval_col="adjusted_timepoint")
# qplot(adjusted_timepoint, log_abund, color = log_abund - 3*log_abund_se > 0, geom="point", data=xxx_interval_abundances) +
#   #geom_ribbon(aes(ymin=log_abund - 2*log_abund_se, ymax=log_abund + 2*log_abund_se), alpha=0.25) +
#   facet_wrap(~cell_group)

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


# Reboot
# -----------------------------------------------------------------------------


state_transition_graph = assemble_timeseries_transitions(wt_ccm_wl,
                                                         start_time=18, stop_time=96,
                                                         interval_col="adjusted_timepoint",
                                                         min_interval = 2,
                                                         log_abund_detection_thresh=-2,
                                                         experiment="GAP14")


plot_state_transition_graph(wt_ccm_wl, state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


# benchmark using contract graph -----------------------------------------------


meso_stg_ctb = contract_state_graph(wt_ccm_wl, state_transition_graph, group_nodes_by = "cell_type_broad")

# how to plot these ?

# score against lineage graph
lit_tree <- load_lineage_tree()

# subset lit_tree to just mesoderm

meso_ctb_groups = unique(meso_cds@colData$cell_type_broad)

lit_tree_meso = lit_tree %>%
  igraph::as_data_frame() %>%
  filter(from %in%  meso_ctb_groups | to %in% meso_ctb_groups) %>%
  igraph::graph_from_data_frame()

# to do: support nodes that don't have a place in the original graph
#' @param truth_graph
#' @param pred_graph
#' @param legend_position
plot_graph_comparison <- function(truth_graph, pred_graph, legend_position = "right") {

  missing_edges = igraph::difference(truth_graph, pred_graph) # blue
  new_edges = igraph::difference(pred_graph, truth_graph) # red
  correct_edges = igraph::intersection(truth_graph, pred_graph) # green


  igraph::E(missing_edges)$type = "missing"
  igraph::E(new_edges)$type = "new"
  igraph::E(correct_edges)$type = "correct"

  # new edges that have nodes in the original

  new_edges = new_edges %>% igraph::as_data_frame() %>%
    filter(from %in% igraph::V(G)$name, to %in% igraph::V(G)$name) %>%
    igraph::graph_from_data_frame()

  union_graph = igraph::union(missing_edges, new_edges, correct_edges) %>%
    igraph::as_data_frame() %>% mutate(type = case_when(
      !is.na(type_1) ~ type_1,
      !is.na(type_2) ~ type_2,
      !is.na(type_3) ~ type_3,
    )) %>% select(-c(type_1, type_2, type_3)) %>%  igraph::graph_from_data_frame()


  # use the correct one as the base for the plot
  G = truth_graph

  # extra nodes just add

  # new_nodes = igraph::V(meso_stg_ctb)$name[!igraph::V(meso_stg_ctb)$name %in% igraph::V(G)$name]
  # shared_nodes = igraph::V(meso_stg_ctb)$name[igraph::V(meso_stg_ctb)$name %in% igraph::V(G)$name]

  # for (n in new_nodes) {
  #   G <- igraph::add_vertices(G, 1, name = n)
  # }

  # G = G %>%
  #   igraph::as_data_frame() %>%
  #   left_join(rbind(igraph::as_data_frame(missing_edges),
  #                   igraph::as_data_frame(correct_edges)), by = c("from", "to")) %>%
  #   igraph::graph_from_data_frame()

  G = union_graph

  G_nel = graph::graphAM(igraph::get.adjacency(G) %>% as.matrix(),
                         edgemode = 'directed') %>%
    as("graphNEL")

  gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot") #, subGList=subgraph_df$subgraph)
  gvizl_coords = cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)

  beziers = lapply(gvizl@renderInfo@edges$splines, function(bc) {
    bc_segments = lapply(bc, Rgraphviz::bezierPoints)
    bezier_cp_df = do.call(rbind, bc_segments) %>% as.data.frame
    colnames(bezier_cp_df) = c("x", "y")
    bezier_cp_df
  })
  bezier_df = do.call(rbind, beziers)
  bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]

  missing_edge_df = igraph::as_data_frame(missing_edges, what = "edges") %>% mutate(edge_name = paste0(from, "~", to))
  correct_edge_df = igraph::as_data_frame(correct_edges, what = "edges") %>% mutate(edge_name = paste0(from, "~", to))
  new_edge_df = igraph::as_data_frame(new_edges, what = "edges") %>% mutate(edge_name = paste0(from, "~", to))

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  edge_colors = c("green4", "blue4", "red4")
  names = c("correct", "missing", "new")

  bezier_df = left_join(bezier_df,
            rbind(missing_edge_df,correct_edge_df, new_edge_df) %>% select(c(type, edge_name)),
            by = c("edge_name"))

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, color=type),
                       data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed")) +
    scale_color_manual(values = edge_colors)

  p = p + ggnetwork::geom_nodelabel(data = g,
                            aes(x, y,
                                # fill = name,
                                label = name),
                            size = node_size)

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)

  return(p)
}

plot_graph_comparison(lit_tree_meso, meso_stg_ctb)


# benchmark with genetics -----------------------------------------------------

# meso_mt_cds = meso_cds[, !grepl("ctrl", colData(meso_cds)$gene_target) ]

colData(meso_cds)$sample = NULL
meso_ccs = new_cell_count_set(meso_cds,
                              sample_group = "embryo",
                              cell_group = "cell_type_sub")

tbx16_ccm = fit_genotype_ccm(genotype = "tbx16",
                            ccs = meso_ccs,
                            ctrl_ids = c("ctrl-uninj", "ctrl-inj"),
                            # multiply = F,
                            # whitelist = wt_state_transition_graph %>% igraph::as_data_frame(),
                            sparsity_factor = 0.1)

tbx16_eff = collect_genotype_effects(ccm = tbx16_ccm,
                                    timepoint = 36,
                                    expt = "GAP16")

plot_contrast(tbx16_ccm, tbx16_eff, plot_edges = "none",
              switch_label = "cell_type_broad", repel_labels = T)


# plot on the graph

parent_child_foldchanges(lit_tree_meso, tbx16_eff)


# plot_FC_graph <- function() {
#
#   cond_b_v_a_tbl = cond_b_v_a_tbl %>%
#     mutate(delta_log_abund = ifelse(delta_q_value < qval_threshold, delta_log_abund, 0))
#
#   G = lit_tree_meso
#   G_nel = graph::graphAM(igraph::get.adjacency(G) %>% as.matrix(),
#                          edgemode = 'directed') %>%
#     as("graphNEL")
#
#   gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot")
#   gvizl_coords = cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)
#
#   beziers = lapply(gvizl@renderInfo@edges$splines, function(bc) {
#     bc_segments = lapply(bc, Rgraphviz::bezierPoints)
#     bezier_cp_df = do.call(rbind, bc_segments) %>% as.data.frame
#     colnames(bezier_cp_df) = c("x", "y")
#     bezier_cp_df
#   })
#   bezier_df = do.call(rbind, beziers)
#   bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]
#
#   g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)
#
#   g = left_join(g, cond_b_v_a_tbl, by = c("name" = "cell_group"))
#
#   p <- ggplot() +
#     ggplot2::geom_path(aes(x, y, group=edge_name),
#                        data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))
#
#   p = p + ggnetwork::geom_nodelabel(data = g,
#                                     aes(x, y,
#                                         fill = delta_log_abund,
#                                         label = name),
#                                     size = node_size) +
#     scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
#
#   p = p + scale_size_identity() +
#     monocle3:::monocle_theme_opts() +
#     ggnetwork::theme_blank() +
#     theme(legend.position=legend_position)
#
# }


