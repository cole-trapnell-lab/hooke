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

# ----------------------------------------------------------------------------

wt_cds = readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")

wt_cds = cluster_cells(wt_cds, resolution = 1e-4) #257 clusters

plot_cells(wt_cds, color_cells_by = "cluster")

colData(wt_cds)$cluster = as.character(clusters(wt_cds))
wt_cds = wt_cds[, is.na(colData(wt_cds)$embryo) == FALSE]
colData(wt_cds)$cell_type = colData(wt_cds)$cell_type_sub

staging_df = colData(wt_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))

staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
staging_df$predicted_timepoint = predict(staging_model)

colData(wt_cds)$adjusted_timepoint = staging_df$predicted_timepoint

wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "embryo",
                            cell_group = "cluster")
wt_main_model_formula_str = "~ splines::ns( adjusted_timepoint, knots= c(24,48,60), Boundary.knots=c(20,92) )"

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 main_model_formula_str = wt_main_model_formula_str,
                                 whitelist = initial_pcor_graph(wt_ccs))

wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=1)


state_transition_graph = assemble_timeseries_transitions(wt_ccm_wl,
                                                         start_time=18, stop_time=96,
                                                         interval_col="adjusted_timepoint",
                                                         min_interval = 2,
                                                         log_abund_detection_thresh=-2,
                                                         experiment="GAP14")


plot_state_transition_graph(wt_ccm_wl, state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


stg_ctb = contract_state_graph(wt_ccm_wl,
                               state_transition_graph,
                               group_nodes_by = "cell_type_broad")

wt_ccm_wl_contract = contract_ccm(wt_ccm_wl, group_nodes_by = "cell_type_broad")


plot_state_transition_graph(wt_ccm_wl_contract,
                            stg_ctb %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


# score against lineage graph
lit_tree <- load_lineage_tree()

plot_state_transition_graph(wt_ccm_wl_contract,
                            lit_tree %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")



# now done by sub umaps --------------------------------------------------------

subset_gap <- function(cds, major_group) {

  sub_cds = cds[, colData(cds)$major_group == major_group]

  reducedDims(sub_cds)[["UMAP"]] = sub_cds@colData %>% as.data.frame %>%
    select(subumap3d_1, subumap3d_2, subumap3d_3) %>%
    as.matrix

  return(sub_cds)
}

make_state_transition_graph <- function(cds)  {

  cds = cluster_cells(cds, resolution = 1e-4)
  colData(cds)$cluster = as.character(clusters(cds))
  ccs = new_cell_count_set(cds,
                           sample_group = "embryo",
                           cell_group = "cluster")

  main_model_formula_str = "~ splines::ns( adjusted_timepoint, knots= c(24,48,60), Boundary.knots=c(20,92) )"

  ccm_wl = new_cell_count_model(ccs,
                                main_model_formula_str = main_model_formula_str,
                                whitelist = initial_pcor_graph(ccs))

  ccm_wl = select_model(ccm_wl, criterion = "EBIC", sparsity_factor=1)


  state_transition_graph = assemble_timeseries_transitions(ccm_wl,
                                                           start_time=18, stop_time=96,
                                                           interval_col="adjusted_timepoint",
                                                           min_interval = 2,
                                                           log_abund_detection_thresh=-2,
                                                           experiment="GAP14")

  return(list(ccs = ccs,
              ccm_wl = ccm_wl,
              state_transition_graph = state_transition_graph))

}


# mesoderm --------------------------------------------------------------------


meso_cds = subset_gap(wt_cds, major_group = "mesoderm")

meso_results = make_state_transition_graph(meso_cds)

plot_state_transition_graph(meso_results$ccm_wl,
                            meso_results$state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

meso_contract = igraph::union(contract_state_graph(meso_results$ccm_wl,
                                                  meso_results$state_transition_graph,
                                                  group_nodes_by = "cell_type_broad"))

meso_ccm_contract = contract_ccm(meso_results$ccm_wl, group_nodes_by = "cell_type_broad")

plot_state_transition_graph(meso_ccm_contract,
                            meso_contract %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

# CNS -------------------------------------------------------------------------

cns_cds = subset_gap(wt_cds, major_group = "CNS")

cns_results = make_state_transition_graph(cns_cds)


plot_state_transition_graph(cns_results$ccm_wl,
                            cns_results$state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

cns_contract = igraph::union(contract_state_graph(cns_results$ccm_wl,
                                                  cns_results$state_transition_graph,
                                                  group_nodes_by = "cell_type_broad"))

cns_ccm_contract = contract_ccm(cns_results$ccm_wl, group_nodes_by = "cell_type_broad")

plot_state_transition_graph(cns_ccm_contract,
                            cns_contract %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

# periderm / other ------------------------------------------------------------

peri_cds = subset_gap(wt_cds, major_group = "periderm-other")

peri_results = make_state_transition_graph(peri_cds)

plot_state_transition_graph(peri_results$ccm_wl,
                            peri_results$state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

peri_contract = igraph::union(contract_state_graph(peri_results$ccm_wl,
                                                   peri_results$state_transition_graph,
                                                   group_nodes_by = "cell_type_broad"))

peri_ccm_contract = contract_ccm(peri_results$ccm_wl, group_nodes_by = "cell_type_broad")

plot_state_transition_graph(peri_ccm_contract,
                            peri_contract %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

# mesenchyme-fin --------------------------------------------------------------

mesen_cds = subset_gap(wt_cds, major_group = "mesenchyme-fin")

mesen_results = make_state_transition_graph(mesen_cds)

plot_state_transition_graph(mesen_results$ccm_wl,
                            mesen_results$state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

mesen_contract = igraph::union(contract_state_graph(mesen_results$ccm_wl,
                                                    mesen_results$state_transition_graph,
                                                    group_nodes_by = "cell_type_broad"))

mesen_ccm_contract = contract_ccm(mesen_results$ccm_wl, group_nodes_by = "cell_type_broad")

plot_state_transition_graph(mesen_ccm_contract,
                            mesen_contract %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")

# combine all sub umap graphs

combined_state_transition_graph = igraph::union(meso_results$state_transition_graph,
                                                 peri_results$state_transition_graph,
                                                 cns_results$state_transition_graph,
                                                 mesen_results$state_transition_graph)

plot_state_transition_graph(wt_ccm_wl,
                            combined_state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


combined_cell_transition_graph = igraph::union(meso_contract,
                                               mesen_contract,
                                               cns_contract,
                                               peri_contract)


# to fix -- incompatible nodes
plot_state_transition_graph(wt_ccm_wl_contract,
                            combined_cell_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


