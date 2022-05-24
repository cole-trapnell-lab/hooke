library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)

setwd("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model/examples/")
devtools::load_all("/Users/maddyduran/OneDrive/UW/Trapnell/hooke_split_model")

# kidney_cds = readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/R_objects/kidney.cds.cole.RDS")
# kidney_cds = detect_genes(kidney_cds)
#
# # Drop clusters that are likely to be multiplets that got stuck together during the initial labeling.
# kidney_cds = kidney_cds[,clusters(kidney_cds) %in% c(10) == FALSE]
# colData(kidney_cds)$orig_cluster = colData(kidney_cds)$cluster
# set.seed(42)
#
# kidney_cds = cluster_cells(kidney_cds, resolution=1e-3, random_seed=42)
# colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
# colData(kidney_cds)$new_cluster = colData(kidney_cds)$cluster
#
# colData(kidney_cds)$experiment = colData(kidney_cds)$expt
# colData(kidney_cds)$sample = NULL

# colData(kidney_cds)$cell_type = colData(kidney_cds)$kidney.celltype
# colData(kidney_cds)$genotype = colData(kidney_cds)$gene_target1
# colData(kidney_cds)$genotype[colData(kidney_cds)$genotype == "ctrl"] = "wt"
# colData(kidney_cds)$Size_Factor = size_factors(kidney_cds)

# pronephros_classifier <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pronephros_classifier.RDS")
# colData(kidney_cds)$garnett_cluster = clusters(kidney_cds)
# kidney_cds = classify_cells(kidney_cds, pronephros_classifier, db="none", cluster_extend=TRUE)

# colData(kidney_cds)$segment = case_when(
#   colData(kidney_cds)$cell_type %in% c("Renal progenitors") ~ -1,
#   colData(kidney_cds)$cell_type %in% c("Podocyte") ~ 0,
#   colData(kidney_cds)$cell_type %in% c("Neck") ~ 1,
#   colData(kidney_cds)$cell_type %in% c("Proximal Convoluted Tubule") ~ 2,
#   colData(kidney_cds)$cell_type %in% c("Proximal Straight Tubule") ~ 3,
#   colData(kidney_cds)$cell_type %in% c("Distal Early") ~ 4,
#   colData(kidney_cds)$cell_type %in% c("Corpuscles of Stannius") ~ 5,
#   colData(kidney_cds)$cell_type %in% c("Distal Late") ~ 6,
#   colData(kidney_cds)$cell_type %in% c("Multiciliated cells")  ~ 7,
#   colData(kidney_cds)$cell_type %in% c("Cloaca")  ~ 8
# )

# kidney_cds = kidney_cds[,is.na(colData(kidney_cds)$Oligo) == FALSE &
#                           is.na(colData(kidney_cds)$timepoint.1) == FALSE &
#                           colData(kidney_cds)$timepoint <= 48]
#
# kid_cds = readRDS("~/OneDrive/UW/Trapnell/fish-crew/updated_R_objects/all_Kidney_cells.RDS")
# kidney_df = colData(kid_cds) %>%
#   as.data.frame() %>%
#   mutate(cell_type_sub = case_when(
#     cluster == 1 ~ "PCT 1",
#     cluster == 2 ~ "(late) Distal late",
#     cluster == 3 ~ "Distal early",
#     cluster == 4 ~ "(early) Distal late",
#     cluster == 5 ~ "PST",
#     cluster == 6 ~ "Early podocyte",
#     cluster == 7 ~ "Early duct",
#     cluster == 8 ~ "Cloaca",
#     cluster == 9 ~ "Late neck",
#     cluster == 10 ~ "PCT 2",
#     cluster == 11 ~ "Podocyte",
#     cluster == 12 ~ "CS",
#     cluster == 13 ~ "Early neck",
#     cluster == 14 ~ "MCCs (duct)",
#     cluster == 15 ~ "PCT 3",
#     cluster == 16 ~ "Distal early",
#     TRUE ~ NA_character_
#   )) %>%
#   select(cell, cell_type_sub) %>%
#   mutate(cell_type_broad = "pronephros",
#          tissue = "Kidney",
#          germ_layer = "mesoderm",
#          major_group = "periderm-other")
#
# kid_coldata = colData(kidney_cds) %>%
#   as.data.frame %>%
#   select(-c(cell_type_sub, cell_type_broad, germ_layer, tissue)) %>%
#   left_join(kidney_df, by = "cell")
#
# colData(kidney_cds)$cell_type_sub = kid_coldata$cell_type_sub
# colData(kidney_cds)$cell_type_broad = kid_coldata$cell_type_broad
# colData(kidney_cds)$germ_layer = kid_coldata$germ_layer
# colData(kidney_cds)$tissue = kid_coldata$tissue

# plot_cells(kidney_cds, color_cells_by = "cell_type")
# plot_cells(kidney_cds, color_cells_by = "cell_type_sub")

# saveRDS(kidney_cds, "~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pronephros.RDS")

kidney_cds <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pronephros.RDS")

# ----------------------------------------------------------------------------

wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
                      colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP16", "GAP18", "HF4")  ]


colData(wt_cds)$cluster = clusters(wt_cds)

wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "Oligo",
                            cell_group = "cluster")

wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_var="timepoint", interval_start=18, interval_stop=48, num_breaks=4)


undebug(new_cell_count_model)

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 main_model_formula_str = wt_main_model_formula_str,
                                 nuisance_model_formula_str = "~experiment",
                                 whitelist = initial_pcor_graph(wt_ccs))


wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.1)

# do we need to fix select model?


kidney_cell_type_abundances = get_extant_cell_types(wt_ccm_wl, start = 18, stop = 48,
                                                    log_abund_detection_thresh=-2,
                                                    percent_max_threshold=0.01, experiment="GAP14")


ggplot(aes(timepoint, log_abund, color=present_above_thresh), data=kidney_cell_type_abundances) +
  geom_point() + facet_wrap(~cell_group, scale="free_y")


wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.01)

wt_state_transition_graph = assemble_timeseries_transitions(wt_ccm_wl,
                                                            start=18, stop=48,
                                                            interval_col="timepoint",
                                                            min_interval = 2,
                                                            log_abund_detection_thresh=-2,
                                                            min_dist_vs_time_r_sq=0.0,
                                                            experiment="GAP14")

hooke:::plot_path(wt_ccm_wl, path_df = wt_state_transition_graph %>% igraph::as_data_frame(), edge_size=1)

plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")


# wt_state_transition_graph_weighted = assemble_timeseries_transitions(wt_ccm_wl,
#                                                             start=18, stop=48,
#                                                             interval_col="timepoint",
#                                                             min_interval = 2,
#                                                             log_abund_detection_thresh=-2,
#                                                             min_dist_vs_time_r_sq=0.0,
#                                                             experiment="GAP14", weigh_by_pcor=T)
#
# hooke:::plot_path(wt_ccm_wl, path_df = wt_state_transition_graph %>% igraph::as_data_frame(), edge_size=1)
#
# plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph_weighted %>% igraph::as_data_frame(),
#                             color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")

# round 2 ----------------------------------------------------------------------
# augment the graph to include things that are not in the paga, but emerge after

not_in_paga = not_in_paga_graph(wt_ccm_wl)
# "20" "24" "39"

colData(kidney_cds) %>%
  as.data.frame() %>%
  filter(cluster %in% not_in_paga) %>%
  select(cluster, cell_type) %>%
  distinct()
# 20 = CS
# 24 - MCCs ?
# 39 - mix of stuff

colData(kidney_cds) %>%
  as.data.frame() %>% group_by(cluster, cell_type) %>% tally() %>%
  filter(cluster == 4) %>%
  # select(cluster, cell_type) %>%
  distinct()

# timepoint 36 to 48 hpf
kidney_cell_type_abundances %>%
  filter(cell_group %in% not_in_paga) %>%
  select(cell_group, longest_contig_start, longest_contig_end) %>%
  distinct()

# cell_group longest_contig_start longest_contig_end
# <chr>                     <dbl>              <dbl>
#   1 20                           34                 46
#   2 24                           36                 48
#   3 39                           18                 48

# going forward we want to link 20 and 24 back


# -----------------------------------------------------------------------------

aug_pathfinding_graph <- function(ccm,
                                  curr_graph,
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

  timeseries_pathfinding_graph = init_pathfinding_graph(ccm, extant_cell_type_df)

  # find ones not in paga that emerge after the start
  not_in_paga = not_in_paga_graph(ccm)

  not_in_paga_and_emerge = extant_cell_type_df %>%
    filter(cell_group %in% not_in_paga) %>%
    filter(longest_contig_start > start_time) %>%
    pull(cell_group) %>% unique()

  edge_to_each = FALSE
  sparsity = 0.11

  while(!edge_to_each & sparsity >=0.01) {
    sparsity = sparsity-0.01
    ccm = select_model(ccm, criterion = "EBIC", sparsity_factor= sparsity)

    not_in_paga_cov_edges =  hooke:::return_igraph(model(ccm, "reduced")) %>%
      igraph::as_data_frame(what="edges") %>%
      dplyr::rename(pcor=weight) %>% dplyr::filter(pcor != 0.00) %>%
      filter(from %in% not_in_paga_and_emerge | to %in% not_in_paga_and_emerge)
    # print(nrow(not_in_paga_cov_edges))

    # node_metadata = data.frame(id=unique(union(not_in_paga_cov_edges$to, not_in_paga_cov_edges$from)))

    # check if at least one edge to each exists
    edge_to_each = length(setdiff(not_in_paga_and_emerge,
                                  union(not_in_paga_cov_edges$from, not_in_paga_cov_edges$to))) == 0

  }

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
  pathfinding_graph = igraph::as_data_frame(pathfinding_graph) %>%
    mutate(node = ifelse(to %in% not_in_paga_and_emerge, to, from)) %>%
    group_by(node) %>% top_n(1, wt = -weight) %>%
    # add them to the current pathfinding graph
    rbind(igraph::as_data_frame(timeseries_pathfinding_graph)) %>%
    igraph::graph_from_data_frame()


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
                           # experiment = "GAP14")
                           ...)

  covered_G = compute_min_path_cover(ccm, G)

  # add any new edges
  new_edges = igraph::difference(covered_G, curr_graph) %>%
    igraph::as_data_frame() %>%
    filter(from %in% not_in_paga_and_emerge | to %in% not_in_paga_and_emerge) %>%
    igraph::graph_from_data_frame()


  aug_graph = igraph::union(curr_graph, new_edges)



  return(aug_graph)
}

# ------------------------------------------------------------------------------

# run the model with the previous state graph as a whitelist
# wt_ccm_wl_2 = new_cell_count_model(wt_ccs,
#                                    main_model_formula_str = wt_main_model_formula_str,
#                                    nuisance_model_formula_str = "~experiment",
#                                    whitelist = wt_state_transition_graph %>% igraph::as_data_frame())

augmented_graph = aug_pathfinding_graph(wt_ccm_wl,
                                        wt_state_transition_graph,
                                        start_time = 18,
                                        stop_time = 48,
                                        experiment = "GAP14")

plot_path(wt_ccm_wl, path_df = wt_state_transition_graph %>% igraph::as_data_frame(), edge_size=1, x=2, y=3)

plot_path(wt_ccm_wl, path_df = augmented_graph %>% igraph::as_data_frame(), edge_size=1, x=2, y=3)

plot_path(wt_ccm_wl, path_df = augmented_graph %>% igraph::as_data_frame(), edge_size=1)

plot_state_transition_graph(wt_ccm_wl, augmented_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")




# genetics --------------------------------------------------------------------


kid_ccs = new_cell_count_set(kidney_cds[,colData(kidney_cds)$expt %in% c("GAP14", "GAP18", "GAP16")],
                             sample_group = "Oligo",
                             cell_group = "cluster")

smo_ccm = fit_genotype_ccm(genotype = "smo",
                           ccs = kid_ccs,
                           ctrl_ids = c("wt", "ctrl-inj"),
                           multiply = F,
                           whitelist = wt_state_transition_graph %>% igraph::as_data_frame(),
                           sparsity_factor = 0.1)

smo_eff = collect_genotype_effects(ccm = smo_ccm,
                                   timepoint = 36,
                                   expt = "GAP16")

plot_contrast(wt_ccm_wl, smo_eff, plot_edges = "none",
              switch_label = "cell_type_sub", repel_labels = T)


# plot different state transition graphs --------------------------------------

# original by cell type
plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type",
                            group_nodes_by="cell_type",
                            layer_nodes_by="timepoint")

# by log fc
plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame(),
                            cond_b_v_a_tbl = cond_b_v_a_tbl,
                            fc_limits = c(-1,1),
                            color_nodes_by = "cell_type",
                            group_nodes_by="cell_type",
                            layer_nodes_by="timepoint")

# by gene expression
plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame(),
                            cond_b_v_a_tbl = cond_b_v_a_tbl,
                            genes = c("smo"),
                            fc_limits = c(-1,1),
                            fract_expr = 0.01,
                            color_nodes_by = "cell_type",
                            group_nodes_by="cell_type",
                            layer_nodes_by="timepoint")

# test labeling edges

labels = c('test', "test2", "test3")
names(labels) = c("12~3", "17~11", "22~2")

plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type",
                            labels = labels,
                            group_nodes_by="cell_type",
                            layer_nodes_by="timepoint")




