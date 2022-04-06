library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)

#kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/kidney/kidney.cds.cole.RDS")

kidney_cds = detect_genes(kidney_cds)

# assign best celltype column and reduce dims
colData(kidney_cds)$cell_type = colData(kidney_cds)$kidney.celltype
colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
colData(kidney_cds)$genotype = colData(kidney_cds)$gene_target1
colData(kidney_cds)$genotype[colData(kidney_cds)$genotype == "ctrl"] = "wt"
colData(kidney_cds)$Size_Factor = size_factors(kidney_cds)
colData(kidney_cds)$cell_type_dk = case_when(
  colData(kidney_cds)$cluster == 1 ~ "Proximal Convoluted Tubule (1)",
  colData(kidney_cds)$cluster == 2 ~ "Mature Distal late",
  colData(kidney_cds)$cluster == 3 ~ "Immature Proximal Tubule" , #Proximal Convoluted Tubule (1)",
  colData(kidney_cds)$cluster == 4 ~ "Mature Distal Early",
  colData(kidney_cds)$cluster == 5 ~ "Immature Distal late", #
  colData(kidney_cds)$cluster == 6 ~ "Proximal Straight Tubule", #"Early duct",
  colData(kidney_cds)$cluster == 7 ~ "Immature podocyte",
  colData(kidney_cds)$cluster == 8 ~ "Cloaca",
  colData(kidney_cds)$cluster == 9 ~ "Mature neck",
  colData(kidney_cds)$cluster == 10 ~ "Proximal Convoluted Tubule (10)",
  colData(kidney_cds)$cluster == 11 ~ "Podocyte",
  colData(kidney_cds)$cluster == 12 ~ "Immature Neck", #
  colData(kidney_cds)$cluster == 13 ~ "Corpuscles of Stannius",
  colData(kidney_cds)$cluster == 14 ~ "Multiciliated cells",
  colData(kidney_cds)$cluster == 15 ~ "Proximal Convoluted Tubule (15)",
  colData(kidney_cds)$cluster == 16 ~ "Unknown",
  TRUE ~ "Unknown"
)

colData(kidney_cds)$orig_cluster = colData(kidney_cds)$cluster

plot_cells(kidney_cds, color_cells_by="orig_cluster", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_orig_clusters.png", width=4, height=4)



set.seed(42)

kidney_cds = cluster_cells(kidney_cds, resolution=1e-3, random_seed=42)
colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
colData(kidney_cds)$new_cluster = colData(kidney_cds)$cluster
plot_cells(kidney_cds)

colData(kidney_cds)$cell_type_ct = case_when(
  colData(kidney_cds)$cluster %in% c(12,1,3,5,26,21) ~ "Proximal Convoluted Tubule",
  colData(kidney_cds)$cluster %in% c(7,16,13) ~ "Distal Early",
  colData(kidney_cds)$cluster %in% c(6,9,4) ~ "Distal Late", #
  colData(kidney_cds)$cluster %in% c(2,25) ~ "Proximal Straight Tubule", #"Early duct",
  colData(kidney_cds)$cluster %in% c(29,14) ~ "Cloaca",
  colData(kidney_cds)$cluster %in% c(19,31,23,11) ~ "Podocyte",
  colData(kidney_cds)$cluster %in% c(36,22,39,8) ~ "Neck", #
  colData(kidney_cds)$cluster %in% c(20) ~ "Corpuscles of Stannius",
  colData(kidney_cds)$cluster %in% c(28) ~ "Multiciliated cells",
  colData(kidney_cds)$cluster %in% c(10,15,27) ~ "Renal progenitors",
  #colData(kidney_cds)$cluster == 16 ~ "Unknown",
  TRUE ~ "Unknown"
)
plot_cells(kidney_cds, color_cells_by="cell_type_ct")

pronephros_classifier <- train_cell_classifier(cds = kidney_cds,
                                         marker_file = "./examples/pronephros_cell_types.txt",
                                         db="none",#org.Dr.eg.db::org.Dr.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
saveRDS(pronephros_classifier, "pronephros_classifier.RDS")

kidney_cds = classify_cells(kidney_cds, pronephros_classifier, db="none")

colData(kidney_cds)$segment = case_when(
  colData(kidney_cds)$cell_type %in% c("Renal progenitors") ~ -1,
  colData(kidney_cds)$cell_type %in% c("Podocyte") ~ 0,
  colData(kidney_cds)$cell_type %in% c("Neck") ~ 1,
  colData(kidney_cds)$cell_type %in% c("Proximal Convoluted Tubule") ~ 2,
  colData(kidney_cds)$cell_type %in% c("Proximal Straight Tubule") ~ 3,
  colData(kidney_cds)$cell_type %in% c("Distal Early") ~ 4,
  colData(kidney_cds)$cell_type %in% c("Corpuscles of Stannius") ~ 5,
  colData(kidney_cds)$cell_type %in% c("Distal Late") ~ 6,
  colData(kidney_cds)$cell_type %in% c("Multiciliated cells")  ~ 7,
  colData(kidney_cds)$cell_type %in% c("Cloaca")  ~ 8
  )

kidney_markers = c("slc20a1a", # PCT
                   "slc12a3", "pppr1b", # Late Distal late
                   "slc12a1", # Distal early
                   "clcnk", "mecom", "slc12a3",  # Early Distal late
                   "trpm7", "slc13a1", # PST
                   "wt1a", # Early podocyte
                   "rfx2", # Early duct
                   "aqp3a", "dnase1l4.1", "evx1", # Cloaca
                   "pax2a", "rfx2", # Late neck
                   "slc20a1a", "pdzk1", # PCT*
                   "wt1b", # Podocyte
                   "stc1l", # CS Early neck
                   "rfx2", "pax2a", "odf3b", #Early neck
                   "odf3b", "rfx2b", # MCCs (duct)
                   "slc4a4a", "slc4a2a", "slc26a2" # PCT*
) %>% unique
plot_cells(kidney_cds, genes=kidney_markers) +
  ggsave("kidney_cell_marker_genes.png", width=12, height=12)


early_kidney_markers = c("hoxb4a", "hoxb7a", "lhx1a", "lhx1b", "osr1",  # AIM
                         "hoxa11a", "hoxc11a", "hoxc11b", "hoxd11a", "hoxd11b", "osr1", "eya1",  # PIM
                         "pax2a", "pax2b", "pax8", # MM
                         "hoxb7a", "pax2b", "pax8", "lhx1a", "lhx1b", # WD
                         "emx2", "gata3", "ret", "wnt11", "wnt9b", "lhx1a", "lhx1b", "bmp7a", "bmp7b", "gfra1a", #UB
                         "cyp26a1", "zulu", "aldh1a2"
) %>% unique
plot_cells(kidney_cds, genes=early_kidney_markers) +
  ggsave("kidney_early_marker_genes.png", width=12, height=12)


colData(kidney_cds)$experiment = colData(kidney_cds)$expt
colData(kidney_cds)$sample = NULL

plot_cells(kidney_cds, color_cells_by="cluster", show_trajectory_graph=FALSE)



plot_cells(kidney_cds, color_cells_by="segment", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_segment.png", width=7, height=6)


plot_cells(kidney_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_cell_types.png", width=4, height=4)

plot_cells(kidney_cds, color_cells_by="cell_type_ct", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_cell_types_ct.png", width=4, height=4)


#plot_cells(kidney_cds, color_cells_by="gene_target1", show_trajectory_graph=FALSE) + facet_wrap(~gene_target1)

#plot_cells(cds)

kidney_cds = kidney_cds[,is.na(colData(kidney_cds)$Oligo) == FALSE & is.na(colData(kidney_cds)$timepoint.1) == FALSE & colData(kidney_cds)$timepoint <= 48]


plot_cells(kidney_cds, color_cells_by="timepoint.1", show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #theme(legend.position="none") +
  ggsave("kidney_time.png", width=7, height=6)

###### Stuff from Maddy:


wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
                      colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP18", "HF4")  ]


colData(wt_cds)$cluster = clusters(wt_cds)


wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "Oligo",
                            cell_group = "cluster")


get_valid_origins <- function(wt_ccm){
  timepoints = seq(18, 48, 2)

  timepoint_pred_df = get_timepoint_pred(wt_ccm)
  timepoint_pred_df = timepoint_pred_df %>% group_by(cell_group) %>% mutate(max_abundance = max(exp(log_abund)),
                                                                            percent_max_abund = exp(log_abund) / max_abundance)
  cell_types_present = timepoint_pred_df %>% filter(percent_max_abund > 0.01)
  valid_origins = tibble(cell_group = unique(colData(wt_ccm@ccs)$cell_group))

  filter_origins <- function(ct, cell_types_present){
    ct_times = cell_types_present %>% filter(cell_group == ct) %>% pull(timepoint) %>% unique()
    possible_origins = cell_types_present %>% filter(timepoint %in% ct_times & cell_group != ct) %>% pull(cell_group) %>%unique
    return (possible_origins)
  }
  valid_origins = valid_origins %>% mutate(possible_origins = purrr::map(.f = filter_origins,
                                                                         .x = cell_group,
                                                                         cell_types_present=cell_types_present))
  return(valid_origins)
}


wt_time_start = 18
wt_time_stop = 48
num_time_breaks = 3
time_breakpoints = seq(wt_time_start, wt_time_stop, length.out=num_time_breaks)
time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")

wt_ccm  = new_cell_count_model(wt_ccs,
                               main_model_formula_str = wt_main_model_formula_str,
                               nuisance_model_formula_str = "~experiment",
                               whitelist = NULL )

valid_origins = get_valid_origins(wt_ccm)
paga_edges = get_paga_graph(wt_ccs@cds) %>% igraph::as_data_frame() %>% as_tibble()

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 main_model_formula_str = wt_main_model_formula_str,
                                 nuisance_model_formula_str = "~experiment",
                                 whitelist = wl )

wt_ccm_01 = select_model(wt_ccm, criterion = "EBIC", sparsity_factor=0.1)
wt_ccm_02 = select_model(wt_ccm, criterion = "EBIC", sparsity_factor=0.2)

wt_ccm_wl_001 = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.01)
wt_ccm_wl_01 = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.1)
wt_ccm_wl_02 = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.2)

wt_ccm_wl_bl = new_cell_count_model(wt_ccs,
                                    main_model_formula_str = wt_main_model_formula_str,
                                    nuisance_model_formula_str = "~experiment",
                                    whitelist = wl,
                                    blacklist = bl)



get_timepoint_pred <- function(wt_ccm) {

  timepoint_pred_df = tibble(timepoint=seq(18,48))

  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      function(tp){ estimate_abundances(wt_ccm, tibble(timepoint=tp, experiment="GAP14"))}, NA_real_),
      .x = timepoint)) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  cell_type_assignments = colData(wt_ccm@ccs@cds) %>%
    as.data.frame %>%
    dplyr::count(cluster, cell_type) %>%
    group_by(cluster) %>% slice_max(n) %>%
    dplyr::select(cell_group=cluster, cell_type)

  timepoint_pred_df = left_join(timepoint_pred_df, cell_type_assignments)
  timepoint_pred_df = timepoint_pred_df %>% mutate(cell_group_label = paste(cell_type, " (", cell_group, ")", sep=""))

  return(timepoint_pred_df)

}


timepoint_pred = get_timepoint_pred(wt_ccm)
timepoint_pred_01 = get_timepoint_pred(wt_ccm_01)
timepoint_pred_02 = get_timepoint_pred(wt_ccm_02)

timepoint_pred_wl = get_timepoint_pred(wt_ccm_wl)
timepoint_pred_wl_01 = get_timepoint_pred(wt_ccm_wl_01)
timepoint_pred_wl_02 = get_timepoint_pred(wt_ccm_wl_02)


# -----------------------------------------------------------------------------


plot_contrast_wrapper <- function(wt_ccm, t1, t2) {

  timepoint_pred_df = get_timepoint_pred(wt_ccm)

  plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                           timepoint_pred_df %>% filter(timepoint == t1),
                                           timepoint_pred_df %>% filter(timepoint == t2)),
                scale_shifts_by = "sender",
                q_value_thresh = 0.01)

}

t1 = 18
t2 = 22
plot_contrast_wrapper(wt_ccm, t1, t2)
plot_contrast_wrapper(wt_ccm_01, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_02, t1, t2)

plot_contrast_wrapper(wt_ccm_wl, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_001, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_01, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_02, t1, t2)

# -----------------------------------------------------------------------------


find_union_path <- function(wt_ccm, q_val=0.01, time_diff = 4) {

  timepoints = seq(18, 48, 2)

  timepoint_pred_df = get_timepoint_pred(wt_ccm)

  select_timepoints <- function(timepoint_pred_df, t1, t2)  {
    cond_x = timepoint_pred_df %>% filter(timepoint == t1)
    cond_y = timepoint_pred_df %>% filter(timepoint == t2)
    return(compare_abundances(wt_ccm, cond_x, cond_y))
  }

  times = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2, (t2-t1) >= time_diff) %>%
    mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(neg_rec_edges = purrr::map(.f = hooke:::get_neg_dir_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(pos_rec_edges = purrr::map(.f = hooke:::get_positive_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(path = purrr::map(.f = purrr::possibly(hooke:::get_path, NA_real_),
                             .x = comp_abund,
                             ccm = wt_ccm,
                             q_value_threshold = q_val,
                             origin_policy="best"))

  distinct_edges = times %>%
    select(neg_rec_edges) %>%
    tidyr::unnest(neg_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  pos_edges = times %>%
    select(pos_rec_edges) %>%
    tidyr::unnest(pos_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  paths = times %>% select(path) %>%
    filter(!is.na(path)) %>%
    tidyr::unnest(path) %>%
    group_by(from, to) %>%
    tally() %>%
    mutate(scaled_weight = abs(n) / max(abs(n)))

  return(list(neg_rec_edges = distinct_edges,
              pos_edges = pos_edges,
              paths = paths))

}
#debug(find_union_path)

find_parsimonius_origins <- function(wt_ccm, q_val=0.01, time_diff = 4) {

  timepoints = seq(18, 48, 2)

  timepoint_pred_df = get_timepoint_pred(wt_ccm)

  select_timepoints <- function(timepoint_pred_df, t1, t2)  {
    cond_x = timepoint_pred_df %>% filter(timepoint == t1)
    cond_y = timepoint_pred_df %>% filter(timepoint == t2)
    return(compare_abundances(wt_ccm, cond_x, cond_y))
  }

  times = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2, (t2-t1) >= time_diff) %>%
    mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(neg_rec_edges = purrr::map(.f = hooke:::get_neg_dir_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(pos_rec_edges = purrr::map(.f = hooke:::get_positive_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(path = purrr::map(.f = purrr::possibly(hooke:::get_path, NA_real_),
                             .x = comp_abund,
                             ccm = wt_ccm,
                             q_value_threshold = q_val,
                             origin_policy="best"))

  valid_origins = get_valid_origins(wt_ccm)

  distinct_edges = times %>%
    select(neg_rec_edges) %>%
    tidyr::unnest(neg_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  pos_edges = times %>%
    select(pos_rec_edges) %>%
    tidyr::unnest(pos_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  paths = times %>% select(path) %>%
    filter(!is.na(path)) %>%
    tidyr::unnest(path)

  indirect_paths = paths
  indirect_paths$path_type = "indirect"
  indirect_paths = indirect_paths %>% select(origin, destination, path_type, from, to)

  direct_paths = distinct_edges # these are all negative reciprocal
  direct_paths = direct_paths %>% mutate (origin=from,
                         destination=to,
                         path_type = "direct")
  direct_paths = direct_paths %>% select(origin, destination, path_type, from, to)


  paths = rbind(indirect_paths, direct_paths)

  paths = left_join(paths, valid_origins, by=c("destination"="cell_group"))
  paths = paths %>% filter(origin %in% unlist(possible_origins)) %>% select(-possible_origins) %>% distinct()

  destinations_with_indirect_origins = paths %>%
    group_by(destination) %>%
    summarize(num_origin = n()) %>%
    ungroup %>%
    filter(num_origin_types > 1) %>%
    pull(destination) %>% unique()

  paths = paths %>% filter((destination %in% destinations_with_indirect_origins & path_type == "indirect") |
                            (destination %in% destinations_with_indirect_origins == FALSE & path_type == "direct"))
  paths = paths %>%
    group_by(from, to) %>%
    tally() %>%
    mutate(scaled_weight = abs(n) / max(abs(n)))

  return(list(neg_rec_edges = distinct_edges,
              pos_edges = pos_edges,
              paths = paths))

}
debug(find_parsimonius_origins)

results_00 = find_union_path(wt_ccm)
results_01 = find_union_path(wt_ccm_01)
results_02 = find_union_path(wt_ccm_02)

results_wl = find_union_path(wt_ccm_wl)
results_wl_01 = find_union_path(wt_ccm_wl_01)
results_wl_02 = find_union_path(wt_ccm_wl_02)


hooke:::plot_path(wt_ccm, path_df = results_00$paths)
hooke:::plot_path(wt_ccm, path_df = results_01$paths)
hooke:::plot_path(wt_ccm, path_df = results_02$paths)

hooke:::plot_path(wt_ccm_01, path_df = results_wl$paths)
hooke:::plot_path(wt_ccm_01, path_df = results_wl_01$paths)
hooke:::plot_path(wt_ccm_01, path_df = results_wl_02$paths)




# how many clusters are connected
union(results_wl$paths$to, results_wl$paths$from) %>% unique %>% length()

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
                                        hide_unlinked_nodes=TRUE){

  # # remove any edge duplicates
  # ade = edges %>%
  #   group_by(from,to) %>%
  #   arrange(edge_type) %>% # needs to deal w this later, only ggnetwork can't handle it
  #   slice(1) %>% # currently chooses the activator if multiple
  #   ungroup() %>%
  #   mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))

  edges = hooke:::distance_to_root(edges)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs)[,c("cell_group",
                                            color_nodes_by,
                                            label_nodes_by,
                                            group_nodes_by,
                                            layer_nodes_by)] %>%
    as.data.frame
  node_metadata = tibble(id=unique(cell_group_metadata$cell_group))
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

  #state_order = path %>% select(to, distance_from_root) %>%
  #  rbind(data.frame("to"="4", distance_from_root=0)) %>%
  #  mutate("cell_group" = paste0("cluster_", to))

  # level_df = data.frame("name" = V(G)$name) %>%
  #   #left_join(gene_id_level, by = c("name" = "gene_id")) %>%
  #   left_join(state_order, by = "cell_group") %>%
  #   group_by(distance_from_root) %>%
  #   mutate(rn = row_number()) %>%
  #   mutate(group = cut(rn, num_levels, labels=F)) %>%
  #   mutate(group_label = as.numeric(distance_from_root) + (group-1)*(0.75/num_levels)) %>%
  #   ungroup() %>%
  #   tibble::column_to_rownames("name")

  # run sugiyama layout
  layers = NULL
  if (is.null(layer_nodes_by) == FALSE) {
    layers=igraph::V(G)$layer_nodes_by
  }
  lay1 <- igraph::layout_with_sugiyama(G, layers=layers)

  g = ggnetwork::ggnetwork(G, layout = lay1$layout, arrow.gap = arrow.gap)

  # add level information
  #g = g %>% left_join(level_df %>% rownames_to_column("id"), by = c("vertex.names"="id"))
  #g = g %>% left_join(regulator_score_df, by = c("vertex.names" = "gene_id") )


  p <- ggplot(mapping = aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
    # draw activator edges
    ggnetwork::geom_edges(data = g,
               arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))
  if (is.null(group_nodes_by) == FALSE){
    p = p + ggforce::geom_mark_rect(aes(fill = group_nodes_by, label=group_nodes_by, filter = group_nodes_by %in% unlabeled_groups == FALSE), size=0, data=g)
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
    ggnetwork::theme_blank()
  return(p)
}
undebug(plot_state_transition_graph)
plot_state_transition_graph(wt_ccm_01, results_wl$paths, color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="segment")
plot_state_transition_graph(wt_ccm_wl_01, results_wl_01$paths, color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="segment")
plot_state_transition_graph(wt_ccm_wl_02, results_wl_02$paths, color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="segment")

# ----------------------------------------------------------------------------



#gp = my_plot_cells(kidney_ccs, color_cells_by = "timepoint")

# edges = distinct_edges %>%
#   add_umap_coords(centroids(wt_ccs))
#
# edges = pos_edges %>%
#   add_umap_coords(centroids(wt_ccs))

my_plot_path <- function(wt_ccs, edges, x = 1, y = 2) {

  gp = hooke:::my_plot_cells(wt_ccs, color_cells_by = "timepoint", x=x, y=y)

  edges = edges %>%
    hooke:::add_umap_coords(centroids(wt_ccs))

  gp +
    geom_segment(data = edges,
                 aes(x = get(paste0("umap_to_", x)),
                     y = get(paste0("umap_to_", y)),
                     xend =  get(paste0("umap_from_", x)),
                     yend =  get(paste0("umap_from_", y)),
                     size = scaled_from_weight,
                 ),
                 color="black") +
    geom_segment(data = edges,
                 aes(x = get(paste0("umap_from_", x)),
                     y = get(paste0("umap_from_", y)),
                     xend = ( get(paste0("umap_to_", x))+ get(paste0("umap_from_", x)) )/2,
                     yend = ( get(paste0("umap_to_", y))+ get(paste0("umap_from_", y)) )/2,
                     size = scaled_from_weight,
                 ),
                 color="black",
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
    scale_size_identity()

  # return(p)

}



# gp = my_plot_cells(wt_ccs, color_cells_by = "timepoint")
my_plot_path(wt_ccs, results_wl$paths) %>%
  ggsave(filename = "kidney_wt_path.png")



#########
# Plot changes over time in wild type

###################

# New anaysis

# # WT kinetics first:
#
# wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
#                       colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP18", "HF4")  ]
#
#
# wt_ccs = new_cell_count_set(wt_cds,
#                          sample_group = "Oligo",
#                          cell_group = "cell_type")
#
#
# wt_ccm  = new_cell_count_model(wt_ccs,
#                             main_model_formula_str = "~splines::ns(timepoint,df=4)",
#                             nuisance_model_formula_str = "~experiment")
#
# #ccm = select_model(ccm, criterion="StARS", sparsity_factor=0.1)
#
# timepoint_pred_df = tibble(timepoint=seq(18,48))
#
# timepoint_pred_df = timepoint_pred_df %>%
#   dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
#     function(tp){ estimate_abundances(wt_ccm, tibble(timepoint=tp, experiment="GAP14"))}, NA_real_), .x = timepoint)) %>% tidyr::unnest()
# qplot(timepoint, log_abund, data=timepoint_pred_df, geom="line", color=cell_group) +
#   geom_ribbon(aes(timepoint, ymin=log_abund-2*log_abund_se, ymax=log_abund+2*log_abund_se), alpha=I(0.1)) +
#   facet_wrap(~cell_group, scales="free_y")
#
# plot_cells(wt_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE) +
#   theme(axis.line=element_blank(),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank()) +
#   facet_wrap(~timepoint)+
#   ggsave("kidney_cell_types_time_faceted.png", width=10, height=10)
#
#
# kidney_cell_type_markers = top_markers(wt_cds, group_cells_by = "cell_type")
#
# #plot_genes_by_group(wt_cds, )
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                  timepoint_pred_df %>% filter(timepoint == 18),
#                                  timepoint_pred_df %>% filter(timepoint == 36)),
#               q_value_thresh = 0.01,
#               sender_cell_groups="Mixed")
#


### Microclustering of the WT:


wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
                      colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP18", "HF4")  ]
colData(wt_cds)$cluster = clusters(wt_cds)

kidney_cell_type_top_markers = top_markers(wt_cds, group_cells_by = "cell_type")


plot_genes_by_group(wt_cds,
                    markers=kidney_cell_type_top_markers %>%
                      group_by(cell_group) %>%
                      arrange(marker_score) %>%
                      slice_tail(n=3) %>%
                      pull(gene_short_name),
                    group_cells_by="cell_type")

# Plot the classic markers
plot_genes_by_group(wt_cds, markers=kidney_markers, group_cells_by="cell_type")

wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "Oligo",
                            cell_group = "cluster")

wt_time_start = 18
wt_time_stop = 48
num_time_breaks = 3
time_breakpoints = seq(wt_time_start, wt_time_stop, length.out=num_time_breaks)
time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")


wt_ccm  = new_cell_count_model(wt_ccs,
                               main_model_formula_str = wt_main_model_formula_str,
                               nuisance_model_formula_str = "~experiment")
wt_ccm = select_model(wt_ccm, criterion = "EBIC", sparsity_factor=0.2)

timepoint_pred_df = tibble(timepoint=seq(18,48))

timepoint_pred_df = timepoint_pred_df %>%
  dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
    function(tp){ estimate_abundances(wt_ccm, tibble(timepoint=tp, experiment="GAP14"))}, NA_real_), .x = timepoint)) %>% tidyr::unnest()

cell_type_assignments = colData(wt_ccs@cds) %>%
  as.data.frame %>%
  dplyr::count(cluster, cell_type) %>%
  group_by(cluster) %>% slice_max(n) %>%
  dplyr::select(cell_group=cluster, cell_type)
timepoint_pred_df = left_join(timepoint_pred_df,cell_type_assignments)
timepoint_pred_df = timepoint_pred_df %>% mutate(cell_group_label = paste(cell_type, " (", cell_group, ")", sep=""))

qplot(timepoint, log_abund, data=timepoint_pred_df, geom="line", color=cell_group) +
  #geom_ribbon(aes(timepoint, ymin=log_abund-2*log_abund_se, ymax=log_abund+2*log_abund_se), alpha=I(0.1)) +
  facet_wrap(~cell_group_label, scales="free_y")


timepoint_pred_df %>% filter(grepl("Podocyte", cell_group_label)) %>%
  mutate(cell_group_label=forcats::fct_relevel(cell_group_label, c("Podocyte (11)", "Podocyte (23)", "Podocyte (31)", "Podocyte (19)"))) %>%
ggplot(aes(timepoint, log_abund, color=cell_group_label)) +
  geom_line() +
  #geom_ribbon(aes(timepoint, ymin=log_abund-2*log_abund_se, ymax=log_abund+2*log_abund_se), alpha=I(0.1)) +
  facet_wrap(~cell_group_label, ncol=1, scales="free_y") + monocle3:::monocle_theme_opts() + theme(legend.position="none") +
  ggsave("podocyte_kinetics.png", width=3, height=4)


timepoint_pred_df %>% filter(grepl("Renal progenitors \\(15\\)|Convoluted", cell_group_label)) %>%
  mutate(cell_group_label=forcats::fct_relevel(cell_group_label, c("Renal progenitors (15)",
                                                                   "Proximal Convoluted Tubule (3)",
                                                                   "Proximal Convoluted Tubule (21)",
                                                                   "Proximal Convoluted Tubule (5)",
                                                                   "Proximal Convoluted Tubule (1)",
                                                                   "Proximal Convoluted Tubule (33)",
                                                                   "Proximal Convoluted Tubule (12)"))) %>%
  ggplot(aes(timepoint, log_abund, color=cell_group_label)) +
  geom_line() +
  #geom_ribbon(aes(timepoint, ymin=log_abund-2*log_abund_se, ymax=log_abund+2*log_abund_se), alpha=I(0.1)) +
  facet_wrap(~cell_group_label, ncol=1, scales="free_y") + monocle3:::monocle_theme_opts() + theme(legend.position="none") +
  ggsave("PCT_kinetics.png", width=3, height=7)


plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                         timepoint_pred_df %>% filter(timepoint == 18),
                                         timepoint_pred_df %>% filter(timepoint == 36)),
              q_value_thresh = 0.01) +
  theme(legend.position="none") +
  ggtitle("WT, 18 vs. 36 hpf") +
  ggsave("wt_pronephros_lfc_18_v_36.png", width=4, height=4)


plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                         timepoint_pred_df %>% filter(timepoint == 18),
                                         timepoint_pred_df %>% filter(timepoint == 24)),
              q_value_thresh = 0.01) +
  theme(legend.position="none") +
  ggtitle("WT, 18 vs. 24 hpf") +
  ggsave("wt_pronephros_lfc_18_v_24.png", width=4, height=4)


plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                         timepoint_pred_df %>% filter(timepoint == 24),
                                         timepoint_pred_df %>% filter(timepoint == 36)),
              q_value_thresh = 0.01) +
  theme(legend.position="none") +
  ggtitle("WT, 24 vs. 36 hpf") +
  ggsave("wt_pronephros_lfc_24_v_36.png", width=4, height=4)


plot_pronephros_transitions <- function(ccm,
                          #umap_centers,
                          cond_b_vs_a_tbl,
                          log_abundance_thresh = -5,
                          scale_shifts_by=c("receiver", "sender", "none"),
                          #cell_group="cluster",
                          edge_size=2,
                          cell_size=1,
                          q_value_thresh = 1.0,
                          group_label_size=2,
                          plot_labels = c("significant", "all", "none"),
                          fc_limits=c(-3,3),
                          sender_cell_groups=NULL,
                          receiver_cell_groups=NULL,
                          plot_edges = TRUE,
                          model_for_pcors="reduced"){

  space_time_summary = colData(ccm@ccs@cds) %>% as.data.frame %>% group_by(cluster) %>% summarise(mean_hpf = mean(timepoint),
                                                                                             mean_segment = mean(segment, na.rm=TRUE))
  cell_type_assignments = colData(ccm@ccs@cds) %>% as.data.frame %>% dplyr::count(cluster, cell_type) %>% group_by(cluster) %>% slice_max(n)
  space_time_summary = left_join(space_time_summary, cell_type_assignments)
  space_time_summary$cell_group = space_time_summary$cluster
  space_time_summary$cluster = NULL
  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))
  space_time_summary = dplyr::left_join(space_time_summary, cond_b_vs_a_tbl, by=c("cell_group"="cell_group"))
  space_time_summary = space_time_summary %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  space_time_summary = space_time_summary %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))

  corr_edge_coords_umap_delta_abund = hooke:::collect_pln_graph_edges(ccm,
                                                              space_time_summary,
                                                              log_abundance_thresh,
                                                              model_for_pcors=model_for_pcors)
  directed_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("directed_to_from", "directed_from_to"))
  undirected_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("undirected"))

  if (is.null(sender_cell_groups) == FALSE){
    directed_edge_df = directed_edge_df %>% dplyr::filter(from %in% sender_cell_groups)
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(from %in% sender_cell_groups | to %in% sender_cell_groups)
  }

  if (is.null(receiver_cell_groups) == FALSE){
    directed_edge_df = directed_edge_df %>% dplyr::filter(to %in% receiver_cell_groups)
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(from %in% receiver_cell_groups | to %in% receiver_cell_groups)
  }

  if (scale_shifts_by == "sender"){
    directed_edge_df = directed_edge_df %>%
      dplyr::group_by(to) %>%
      dplyr::mutate(flow_factor = -pmin(0, pcor),
                    total_weight = sum(flow_factor),
                    scaled_weight  = flow_factor / total_weight)
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }else if (scale_shifts_by == "receiver"){
    directed_edge_df = directed_edge_df %>%
      dplyr::group_by(from) %>%
      dplyr::mutate(flow_factor = -pmin(0, pcor),
                    total_weight = sum(flow_factor),
                    scaled_weight  = flow_factor / total_weight)
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }else{
    directed_edge_df = directed_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }

  gp = ggplot(aes(x=mean_hpf, y=mean_segment), data=space_time_summary)

  if (plot_edges) {
    gp = gp  +
      geom_segment(data = undirected_edge_df,
                   aes(x = to_mean_hpf,
                       y = to_mean_segment,
                       xend=from_mean_hpf,
                       yend = from_mean_segment,
                       size=edge_size * scaled_weight),
                   #size=edge_size / 4,
                   color="lightgray") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = to_mean_hpf,
                       y = to_mean_segment,
                       xend=from_mean_hpf,
                       yend = from_mean_segment,
                       size=edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = to_mean_hpf,
                       y = to_mean_segment,
                       xend=(to_mean_hpf+from_mean_hpf)/2,
                       yend = (to_mean_segment+from_mean_segment)/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = from_mean_hpf,
                       y = from_mean_segment,
                       xend=to_mean_hpf,
                       yend = to_mean_segment,
                       size=edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = from_mean_hpf,
                       y = from_mean_segment,
                       xend=(from_mean_hpf+to_mean_hpf)/2,
                       yend = (from_mean_segment+to_mean_segment)/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_size_identity()
  }

  if (plot_labels != "none") {
    label_df = space_time_summary
    if (plot_labels == "significant")
      label_df = label_df %>% filter(delta_log_abund != 0)
    gp <- gp + ggrepel::geom_label_repel(data = label_df,
                                         mapping = aes(mean_hpf, mean_segment, label=paste(cell_group," - ", cell_type)),
                                         size=I(group_label_size),
                                         fill = "white")
  }
  return(gp)
}
#debug(plot_pronephros_transitions)
plot_pronephros_transitions(wt_ccm,
                            compare_abundances(wt_ccm,
                                               timepoint_pred_df %>% filter(timepoint == 18),
                                               timepoint_pred_df %>% filter(timepoint == 24)),
                            q_value_thresh = 0.1,
                            #q_value_thresh = 1,
                            plot_labels = "all")

plot_pronephros_transitions(wt_ccm,
                            compare_abundances(wt_ccm,
                                               timepoint_pred_df %>% filter(timepoint == 24),
                                               timepoint_pred_df %>% filter(timepoint == 30)),
                            q_value_thresh = 0.1,
                            #scale_shifts_by = "none",
                            #q_value_thresh = 1,
                            plot_labels = "all")

plot_pronephros_transitions(wt_ccm,
                            compare_abundances(wt_ccm,
                                               timepoint_pred_df %>% filter(timepoint == 30),
                                               timepoint_pred_df %>% filter(timepoint == 36)),
                            q_value_thresh = 0.1,
                            #q_value_thresh = 1,
                            plot_labels = "all")

plot_pronephros_transitions(wt_ccm,
                            compare_abundances(wt_ccm,
                                               timepoint_pred_df %>% filter(timepoint == 36),
                                               timepoint_pred_df %>% filter(timepoint == 48)),
                            q_value_thresh = 0.1,
                            #q_value_thresh = 1,
                            plot_labels = "all")

PCT_cds = wt_cds[,colData(wt_cds)$cell_type %in% c("Renal progenitors", "Proximal Convoluted Tubule")]
podocyte_cds = wt_cds[,colData(wt_cds)$cell_type %in% c("Renal progenitors", "Podocyte")]
neck_cds = wt_cds[,colData(wt_cds)$cell_type %in% c("Renal progenitors", "Neck")]

renal_progenitor_markers = c("gfra1a", "lhx1a", "ret", "pax2a", "pax2b", "pax8")
PCT_markers = c("slc20a1a", "slc4a4a", "slc4a2a", "slc26a2")
Podocyte_markers = c("wt1b", "wt1a")
Neck_markers = c("rfx2", "pax2a", "odf3b")

monocle3::plot_percent_cells_positive(PCT_cds[rowData(wt_cds)$gene_short_name %in% union(renal_progenitor_markers, PCT_markers),],
                                      group_cells_by="timepoint")

monocle3::plot_percent_cells_positive(podocyte_cds[rowData(wt_cds)$gene_short_name %in% union(renal_progenitor_markers, Podocyte_markers),],
                                      group_cells_by="timepoint")

monocle3::plot_percent_cells_positive(neck_cds[rowData(wt_cds)$gene_short_name %in% union(renal_progenitor_markers, Neck_markers),],
                                      group_cells_by="timepoint")



agg_expr_mat = monocle3::aggregate_gene_expression(wt_ccs@cds,
                                                   cell_group_df = tibble::rownames_to_column(wt_ccs@metadata[["cell_group_assignments"]]),
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0,
                                                   cell_agg_fun="mean")
agg_coldata = wt_ccs@metadata[["cell_group_assignments"]] %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]

pseudobulk_metadata_summary = colData(wt_ccs) %>% as.data.frame %>%
  group_by(cell_group) %>% dplyr::summarize(timepoint = mean(timepoint),
                                        segment = mean(segment))
agg_coldata = left_join(agg_coldata,
                        pseudobulk_metadata_summary,
                        by=c("cell_group"="cell_group"))
row.names(agg_coldata) = agg_coldata$group_id
agg_coldata = agg_coldata[colnames(agg_expr_mat),]

# gene_set = msigdbr(species = "Danio rerio", subcategory = "GO:MF")
# transcription_regulators = gene_set %>%
#   dplyr::select(gs_id, gene_symbol, gs_name) %>%
#   dplyr::filter(grepl("Transcription", gs_name, ignore.case=TRUE)) %>%
#   pull(gene_symbol) %>% unique %>% sort

pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(wt_ccs@cds) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
#pseudobulk_cds = preprocess_cds(pseudobulk_cds)
#pseudobulk_cds = reduce_dimension(pseudobulk_cds)
#rowData(pseudobulk_cds)$TF = rowData(pseudobulk_cds)$gene_short_name %in% transcription_regulators

#plot_cells(pseudobulk_cds, color_cells_by="cell_group")
pseudobulk_cds = pseudobulk_cds[,colData(pseudobulk_cds)$segment >= 0 &
                                 colData(pseudobulk_cds)$num_cells_in_group > 1 &
                                  is.na(colData(pseudobulk_cds)$segment) == FALSE ]
pseudo_fit = fit_models(pseudobulk_cds,
                        "~splines::ns(segment, df=3) + splines::ns(timepoint, df=3)",
                        #"~1",
                        weights=colData(pseudobulk_cds)$num_cells_in_group,
                        cores=4)
pseudo_coefs = coefficient_table(pseudo_fit)
position_varying_degs = pseudo_coefs %>% filter(grepl("segment", term) & q_value < 0.01) %>% pull(gene_short_name) %>% unique

plot_genes_by_group(wt_cds[,colData(wt_cds)$timepoint == 24 &
                            colData(wt_cds)$segment >= 0 & is.na(colData(wt_cds)$segment) == FALSE],
                    markers = position_varying_degs, group_cells_by="segment",
                    ordering_type="none")


#late_cell_cds = wt_cds[,colData(wt_cds)$timepoint == 36]
monocle3::plot_percent_cells_positive(wt_cds[grepl("hoxb", rowData(wt_cds)$gene_short_name)],
                                      group_cells_by="segment")


### Genetic analysis


ccs = new_cell_count_set(kidney_cds[,colData(kidney_cds)$experiment %in% c("GAP14", "GAP18", "GAP16")],
                         sample_group = "Oligo",
                         cell_group = "cluster")


# ccm  = new_cell_count_model(ccs,
#                             main_model_formula_str = "~splines::ns(timepoint,df=3) + genotype",
#                             nuisance_model_formula_str = "1")
#
# ccm = select_model(ccm, criterion="StARS", sparsity_factor=0.1)


fit_genotype_ccm = function(genotype,
                            ccs,
                            ctrl_ids=c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                            num_time_breaks=3){
  subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]

  colData(subset_ccs)$knockout = colData(subset_ccs)$gene_target == genotype
  knockout_time_start = min(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  knockout_time_stop = max(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  subset_ccs = subset_ccs[,colData(subset_ccs)$timepoint >= knockout_time_start & colData(subset_ccs)$timepoint <= knockout_time_stop]
  time_breakpoints = c()
  if (num_time_breaks > 2 & knockout_time_stop > knockout_time_start){
    time_breakpoints = seq(knockout_time_start, knockout_time_stop, length.out=num_time_breaks)
    time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
    main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ") + knockout")
  }else{
    main_model_formula_str = "~ knockout"
  }

  if (length(unique(colData(subset_ccs)$experiment)) > 1)
    nuisance_model_formula_str = "~experiment"
  else
    nuisance_model_formula_str = "~1"

  # FIXME: Should prob use splines
  #  compound_ccm = new_cell_count_model(subset_ccs,
  #                                      main_model_formula_str = "~ splines::ns(compound_dose, knots=c(1,3))")
  genotype_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
                                                       main_model_formula_str = main_model_formula_str,
                                                       #main_model_formula_str = "~ splines::ns(timepoint, knots=c(24, 30, 36)) + knockout",
                                                       #main_model_formula_str = "~ as.factor(timepoint) + knockout",

                                                       nuisance_model_formula_str = nuisance_model_formula_str
                                                       ))
  genotype_ccm = select_model(genotype_ccm, sparsity_factor = 0.2)
  return(genotype_ccm)
}
undebug(fit_genotype_ccm)


genotype_df = colData(ccs@cds) %>% as_tibble() %>% dplyr::select(gene_target) %>% distinct()

crispant_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("smo", "noto", "tbxta", "cdx4", "egr2b", "mafba", "epha4a"))
# Fit a model for each genotype:
#

crispant_models_tbl = crispant_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, ctrl_ids=c("wt", "ctrl-inj"))) #%>%
#tidyr::unnest(states)

noto_mut_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("noto-mut"))
# Fit a model for each genotype:

noto_models_tbl = noto_mut_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, ctrl_ids=c("wt", "ctrl-noto"))) #%>%
#tidyr::unnest(states)


genotype_models_tbl = rbind(noto_models_tbl, crispant_models_tbl)


collect_genotype_effects = function(ccm, timepoint=24, experiment="GAP16"){
  control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=timepoint, experiment=experiment))
  knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=timepoint, experiment=experiment))
  genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
}
#debug(collect_genotype_effects)

genotype_models_tbl = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=18)) #%>%
#tidyr::unnest(states)

plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[5]], genotype_models_tbl$genotype_eff[[5]], q_value_thresh = 0.05)


# noto
plot_contrast(genotype_models_tbl$genotype_ccm[[8]], genotype_models_tbl$genotype_eff[[8]], q_value_thresh = 0.05)+
  theme(legend.position="none") +
  ggtitle("noto vs. ctrl 18hpf") +
  ggsave("noto_pronephros_lfc_18hpf.png", width=4, height=4)


# smo
plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05)+
  theme(legend.position="none") +
  ggtitle("smo vs. ctrl 18hpf") +
  ggsave("smo_pronephros_lfc_18hpf.png", width=4, height=4)


# egr2b
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05)+
  theme(legend.position="none") +
  ggtitle("egr2b vs. ctrl 18hpf") +
  ggsave("egr2b_pronephros_lfc_18hpf.png", width=4, height=4)



#state_transitions = hooke:::collect_pln_graph_edges(ccm, cond_ra_vs_wt_tbl) %>% as_tibble %>%
#  filter(edge_type == "directed_from_to" & to_delta_p_value < 0.05 & from_delta_p_value < 0.05)

#
# collect_state_transitions = function(ccm, dose_contrast){
#   transitions = hooke:::collect_pln_graph_edges(ccm, dose_contrast)
#   transitions = transitions %>% dplyr::filter(edge_type != "undirected")
#   return(transitions)
# }
#
# compound_models_tbl = compound_models_tbl %>%
#   dplyr::mutate(state_transitions = purrr::map2(.f = purrr::possibly(
#     collect_state_transitions, NA_real_), .x = drug_ccm,
#     .y = dose_eff)) #%>%
# #tidyr::unnest(states)
#
# state_transition_effects = compound_models_tbl %>%
#   dplyr::select(catalog_number, pathway_level_1, pathway_level_2, product_name, target, state_transitions) %>%
#   tidyr::unnest(state_transitions)
#
# sig_state_transition_effects = state_transition_effects %>% dplyr::filter(from_delta_q_value < 0.01 & to_delta_q_value < 0.01)

## Let's look at the kinetics of some of the crispants:


collect_genotype_kinetics = function(ccm, experiment="GAP16", start_time=18, stop_time=36){
  timepoint_pred_df = tibble(timepoint=seq(start_time,stop_time))
  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      function(tp){
        control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=tp, experiment=experiment))
        knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=tp, experiment=experiment))
        genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
        return(genotype_comparison_tbl)
      }, NA_real_), .x = timepoint)) %>% tidyr::unnest()
}

genotype_models_kinetics = genotype_models_tbl %>%
  dplyr::mutate(genotype_kinetics = purrr::map(.f = purrr::possibly(
    collect_genotype_kinetics, NA_real_), .x = genotype_ccm)) %>%
  dplyr::select(gene_target, genotype_kinetics) %>% tidyr::unnest()
#tidyr::unnest(states)

cell_type_assignments = colData(ccs@cds) %>%
  as.data.frame %>%
  dplyr::count(cluster, cell_type) %>%
  group_by(cluster) %>% slice_max(n) %>%
  dplyr::select(cell_group=cluster, cell_type)
genotype_models_kinetics = left_join(genotype_models_kinetics,cell_type_assignments)
genotype_models_kinetics = genotype_models_kinetics %>% mutate(cell_group_label = paste(cell_type, " (", cell_group, ")", sep=""))
qplot(timepoint, log_abund_y, data=genotype_models_kinetics, geom="line", color=gene_target) + facet_wrap(~cell_group_label, scales="free_y")


