library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(monocle3)
library(hooke)
library(garnett)
library(msigdbr)

#kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

setwd("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/kidney")
kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/kidney/kidney.cds.cole.RDS")

kidney_cds = detect_genes(kidney_cds)

# Drop clusters that are likely to be multiplets that got stuck together during the initial labeling.
kidney_cds = kidney_cds[,clusters(kidney_cds) %in% c(10) == FALSE]


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
  colData(kidney_cds)$cluster %in% c(4, 11, 3, 19, 12) ~ "Proximal Convoluted Tubule",
  colData(kidney_cds)$cluster %in% c(23, 1, 16) ~ "Distal Early",
  colData(kidney_cds)$cluster %in% c(6, 17, 2, 8, 31) ~ "Distal Late", #
  colData(kidney_cds)$cluster %in% c(15, 5) ~ "Proximal Straight Tubule", #"Early duct",
  colData(kidney_cds)$cluster %in% c(14, 26) ~ "Cloaca",
  colData(kidney_cds)$cluster %in% c(18, 27, 22, 9) ~ "Podocyte",
  colData(kidney_cds)$cluster %in% c(34, 21, 7, 38) ~ "Neck", #
  colData(kidney_cds)$cluster %in% c(20) ~ "Corpuscles of Stannius",
  colData(kidney_cds)$cluster %in% c(25) ~ "Multiciliated cells",
  colData(kidney_cds)$cluster %in% c(28, 13, 10, 29) ~ "Renal progenitors",
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

# Commenting this out so I don't actually overwrite the classifier:
#saveRDS(pronephros_classifier, "pronephros_classifier.RDS")

colData(kidney_cds)$garnett_cluster = clusters(kidney_cds)
kidney_cds = classify_cells(kidney_cds, pronephros_classifier, db="none", cluster_extend=TRUE)

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
  ggsave("kidney_cell_marker_genes_1_2.png", width=12, height=12)

plot_cells(kidney_cds, 1, 3, genes=kidney_markers) +
  ggsave("kidney_cell_marker_genes_1_3.png", width=12, height=12)

plot_cells(kidney_cds, 2, 3, genes=kidney_markers) +
  ggsave("kidney_cell_marker_genes_2_3.png", width=12, height=12)


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

plot_cells(kidney_cds, color_cells_by="cluster_ext_type", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_cell_types_cluster_ext_type.png", width=4, height=4)


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

# Use the Garnett labels
colData(kidney_cds)$cell_type = colData(kidney_cds)$cluster_ext_type

plot_cells(kidney_cds, color_cells_by="timepoint.1", show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
  #theme(legend.position="none") +
ggsave("kidney_time.png", width=7, height=6)

###### Stuff from Maddy:


wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
                      colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP16", "GAP18", "HF4")  ]


colData(wt_cds)$cluster = clusters(wt_cds)


wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "Oligo",
                            cell_group = "cluster")

#wt_time_start = 18
#wt_time_stop = 48
#num_time_breaks = 3
#time_breakpoints = seq(wt_time_start, wt_time_stop, length.out=num_time_breaks)
#time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
#wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")

wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_var="timepoint", interval_start=18, interval_stop=48, num_breaks=4)


wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 main_model_formula_str = wt_main_model_formula_str,
                                 nuisance_model_formula_str = "~experiment",
                                 whitelist = initial_pcor_graph(wt_ccs) )

kidney_cell_type_abundances = get_extant_cell_types(wt_ccm_wl, start = 18, stop = 48,
                                                    log_abund_detection_thresh=-2,
                                                    percent_max_threshold=0.01, experiment="GAP14")
ggplot(aes(timepoint, log_abund, color=present_above_thresh), data=kidney_cell_type_abundances) + geom_point() + facet_wrap(~cell_group, scale="free_y")


wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.01)

# -----------------------------------------------------------------------------


plot_contrast_wrapper <- function(ccm, t1, t2, q_val=0.01, model_for_pcors="reduced") {

  timepoint_pred_df = estimate_abundances_over_interval(ccm, t1, t2, interval_col="timepoint", experiment="GAP14")

  plot_contrast(ccm, compare_abundances(ccm,
                                        timepoint_pred_df %>% filter(timepoint == t1),
                                        timepoint_pred_df %>% filter(timepoint == t2)),
                scale_shifts_by = "sender",
                q_value_thresh = q_val,
                model_for_pcors=model_for_pcors)

}

t1 = 18
t2 = 22
plot_contrast_wrapper(wt_ccm_wl, 24, 36, model_for_pcors="reduced")

wt_state_transition_graph = assemble_timeseries_transitions(wt_ccm_wl,
                                                   start=18, stop=48,
                                                   interval_col="timepoint",
                                                   min_interval = 2,
                                                   log_abund_detection_thresh=-2,
                                                   min_dist_vs_time_r_sq=0.0,
                                                   experiment="GAP14")

#plot_origins(wt_ccm_wl, paths_to_origins, edge_size=0.25) + facet_wrap(~destination)
#plot_origins(wt_ccm_wl, paths_to_origins %>% filter(emerges_at > 18), edge_size=0.25) + facet_wrap(~destination)

#wt_origins = select_origins(wt_ccm_wl, wt_possible_origins, selection_policy = "acceptable-origins")
hooke:::plot_path(wt_ccm_wl, path_df = wt_state_transition_graph %>% igraph::as_data_frame() , edge_size=0.25)

plot_state_transition_graph(wt_ccm_wl, wt_state_transition_graph %>% igraph::as_data_frame() , color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")

# ----------------------------------------------------------------------------
# Differential analysis across space and time in WT
# (right now this is a total mess and may not be worth keeping)

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
plot_pronephros_transitions(wt_ccm_wl,
                            compare_abundances(wt_ccm_wl,
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

# ---------------------------------------------------------------------------

# gene cell_type self parents siblings children interpretation
# XXXX XXXX

marker_ids = rowData(wt_ccm_wl@ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% kidney_markers) %>% pull(id)


cell_type_graph = contract_state_graph(wt_ccm_wl, wt_state_transition_graph, "cell_type")

#debug(classify_genes_over_graph)
debug_marker_ids = rowData(wt_ccm_wl@ccs@cds) %>% as.data.frame() %>% filter(gene_short_name %in% c("prdm1a", "foxa3")) %>% pull(id)
gene_patterns_over_state_graph = classify_genes_over_graph(wt_ccm_wl, wt_state_transition_graph, debug_marker_ids, abs_expr_thresh=1e-3)

#gene_patterns_over_state_graph = classify_genes_over_graph(wt_ccm_wl, wt_state_transition_graph, marker_ids, abs_expr_thresh=1e-3)
gene_patterns_over_state_graph = gene_patterns_over_state_graph %>% tidyr::unnest(gene_classes) %>% tidyr::unnest(interpretation)
gene_patterns_over_state_graph = left_join(gene_patterns_over_state_graph,
                                           rowData(wt_ccm_wl@ccs@cds) %>%
                                             as_tibble %>%
                                             select(id, gene_short_name), by=c("gene_id"="id"))

cell_type_graph = contract_state_graph(wt_ccm_wl, wt_state_transition_graph, "cell_type")

# Edit the graph to fix up the errors that currently exist:
cell_type_graph = igraph::add_edges(cell_type_graph,
                                   c("Renal progenitors", "Neck",
                                     "Renal progenitors", "Corpuscles of Stannius",
                                     "Renal progenitors", "Multiciliated cells",
                                     "Renal progenitors", "Cloaca",
                                     "Renal progenitors", "Proximal Straight Tubule"))
#cell_type_graph = igraph::delete_edges(cell_type_graph,
#                                       c("Proximal Straight Tubule|Proximal Convoluted Tubule",
#                                         "Proximal Convoluted Tubule|Proximal Straight Tubule",
#                                         "Podocyte|Proximal Convoluted Tubule"))
cell_type_graph = igraph::delete_edges(cell_type_graph,
                                       c("Proximal Straight Tubule|Proximal Convoluted Tubule"))

cell_type_graph = igraph::delete_vertices(cell_type_graph, "Unknown")

gene_id_sample = rowData(wt_ccm_wl@ccs@cds) %>% as.data.frame %>% dplyr::sample_n(250) %>% pull(id)
gene_patterns_over_cell_type_graph = classify_genes_over_graph(wt_ccm_wl,
                                                               cell_type_graph,
                                                               #gene_id_sample,
                                                               #marker_ids,
                                                               group_nodes_by="cell_type",
                                                               log_fc_thresh = 1,
                                                               abs_expr_thresh=1e-3,
                                                               cores=4)
gene_patterns_over_cell_type_graph = gene_patterns_over_cell_type_graph %>% tidyr::unnest(gene_classes)  %>% tidyr::unnest(interpretation)
gene_patterns_over_cell_type_graph = left_join(gene_patterns_over_cell_type_graph,
                                           rowData(wt_ccm_wl@ccs@cds) %>%
                                             as_tibble %>%
                                             select(id, gene_short_name), by=c("gene_id"="id"))

# Gut check a few patterns:
sel_activated = gene_patterns_over_cell_type_graph %>%
  filter(interpretation == "Selectively activated") %>%
  dplyr::sample_n(6) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=sel_activated)

spec_activated = gene_patterns_over_cell_type_graph %>%
  filter(interpretation == "Specifically activated") %>%
  dplyr::sample_n(6) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=spec_activated)

mlp = gene_patterns_over_cell_type_graph %>%
  filter(interpretation == "Selectively maintained") %>%
  dplyr::sample_n(6) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=mlp)

sel_activated = gene_patterns_over_cell_type_graph %>%
  filter(interpretation == "Selectively deactivated") %>%
  dplyr::sample_n(6) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=sel_activated)


cloaca_sel_activated = gene_patterns_over_cell_type_graph %>%
  filter(cell_state == "Cloaca" & interpretation == "Selectively activated") %>%
  dplyr::sample_n(5) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=cloaca_sel_activated)

de_sel_activated = gene_patterns_over_cell_type_graph %>%
  filter(cell_state == "Distal Early" & interpretation == "Selectively activated") %>%
  dplyr::sample_n(5) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=de_sel_activated)

pod_sel_activated = gene_patterns_over_cell_type_graph %>%
  filter(cell_state == "Podocyte" & interpretation == "Selectively activated") %>%
  dplyr::sample_n(5) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=pod_sel_activated)

de_sel_deactivated = gene_patterns_over_cell_type_graph %>%
  filter(cell_state == "Distal Early" & interpretation == "Selectively deactivated") %>%
  dplyr::sample_n(5) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=de_sel_deactivated)

mlp = gene_patterns_over_cell_type_graph %>%
  filter(interpretation == "MLP") %>%
  dplyr::sample_n(5) %>% pull(gene_short_name)
plot_cells(kidney_cds, genes=mlp)

pattern_summary = gene_patterns_over_cell_type_graph %>% dplyr::group_by(cell_state, interpretation) %>% tally()
ggplot(aes(cell_state, n, fill=interpretation), data=pattern_summary %>%
         filter(interpretation != "Absent")) +
  geom_bar(stat="identity") + coord_flip()

#gene_set = msigdbr(species = "Homo sapiens", subcategory = "GO:MF")
#gene_set = msigdbr(species = "Homo sapiens", category = "H")
gene_set = msigdbr(species = "Danio rerio", subcategory = "GO:BP")
gene_set_list = split(x = gene_set$gene_symbol, f = gene_set$gs_name)
gene_set_df = gene_set %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

calc_pathway_enrichment_on_state_specific_genes <- function(gene_df, msigdbr_t2g, sig_thresh = 0.1, ...){
  #gene_set_list = split(x = pathways$gene_symbol, f = pathways$gs_name)
  #gene_ranking = gene_lfc_tbl %>% pull(gene_state_score)
  #names(gene_ranking) = gene_lfc_tbl %>% pull(gene_short_name)
  #pb$tick()
  #gsea_res = fgsea(pathways=gene_set_list, stats=gene_ranking, ...) %>% as_tibble()
  #gsea_res = gsea_res %>% filter(padj < sig_thresh)
  #msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
  gene_symbols_vector = gene_df$gene_short_name
  enrich_res = clusterProfiler::enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...) %>% as_tibble()
  gc()
  return(enrich_res)
}
undebug(calc_pathway_enrichment_on_state_specific_genes)

gene_universe = gene_patterns_over_cell_type_graph %>% filter(interpretation != "Absent") %>% pull(gene_short_name) %>% unique()
state_pattern_pathways = gene_patterns_over_cell_type_graph %>%
  filter(interpretation %in%  c("Selectively activated", "MLP")) %>%
  select(cell_state, interpretation, gene_short_name) %>%
  group_by(cell_state, interpretation) %>% tidyr::nest(data=gene_short_name) %>%
  dplyr::mutate(pathways = purrr::map(.f = purrr::possibly(
    calc_pathway_enrichment_on_state_specific_genes, NA_real_), .x = data, gene_set_df, sig_thresh=1, universe=gene_universe)) %>%
  tidyr::unnest(pathways)
#tidyr::unnest(degs)


gene_set_mf = msigdbr(species = "Danio rerio", subcategory = "GO:MF")
transcription_regulators = gene_set_mf %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("Transcription", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

signaling_genes = gene_set_mf %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("Signaling", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

specific_tfs = gene_patterns_over_cell_type_graph %>%
  filter(interpretation %in%  c("Specifically activated",
                                "Specifically upregulated",
                                "Specifically maintained"
                                )
         & gene_short_name %in% transcription_regulators) %>%
  pull(gene_short_name) %>% unique
plot_genes_by_group(wt_ccm_wl@ccs@cds, markers =specific_tfs, group_cells_by = "cell_type",ordering_type="maximal_on_diag", color_by_group=TRUE)

selective_signals = gene_patterns_over_cell_type_graph %>%
  filter(interpretation %in%  c("Selectively activated", "MLP") & gene_short_name %in% signaling_genes) %>%
  pull(gene_short_name) %>% unique
plot_genes_by_group(wt_ccm_wl@ccs@cds, markers =selective_signals, group_cells_by = "cell_type",ordering_type="maximal_on_diag", color_by_group=TRUE)


selective_genes = gene_patterns_over_cell_type_graph %>%
  filter(interpretation %in%  c("Selectively activated")) %>%
  pull(gene_short_name) %>% unique
selective_genes = selective_genes[grepl("si:", selective_genes) == FALSE]
plot_genes_by_group(wt_ccm_wl@ccs@cds, markers =selective_genes[sample(length(selective_genes), 10)], group_cells_by = "cell_type",ordering_type="maximal_on_diag", color_by_group=TRUE)

# ----------------------------------------------------------------------------
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
                            prior_state_transition_graph=NULL,
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

  if (is.null(prior_state_transition_graph) == FALSE){
    wt_prior_whitelist = prior_state_transition_graph %>% igraph::as_data_frame()
  }else{
    wt_prior_whitelist = NULL
  }

  # FIXME: Should prob use splines
  #  compound_ccm = new_cell_count_model(subset_ccs,
  #                                      main_model_formula_str = "~ splines::ns(compound_dose, knots=c(1,3))")
  genotype_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
                                                       main_model_formula_str = main_model_formula_str,
                                                       #main_model_formula_str = "~ splines::ns(timepoint, knots=c(24, 30, 36)) + knockout",
                                                       #main_model_formula_str = "~ as.factor(timepoint) + knockout",

                                                       nuisance_model_formula_str = nuisance_model_formula_str,
                                                       whitelist = wt_prior_whitelist
                                                       ))
  # FIXME: we should maybe be pulling out the sparsity_factor used for the WT prior model and using that here rather
  # than hardcoding?
  genotype_ccm = select_model(genotype_ccm, sparsity_factor = 0.01)
  return(genotype_ccm)
}
undebug(fit_genotype_ccm)


genotype_df = colData(ccs@cds) %>% as_tibble() %>% dplyr::select(gene_target) %>% distinct()

crispant_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("smo", "noto", "tbxta", "cdx4", "egr2b", "mafba", "epha4a"))
# Fit a model for each genotype:
#

crispant_models_tbl = crispant_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, wt_state_transition_graph, ctrl_ids=c("wt", "ctrl-inj"))) #%>%
#tidyr::unnest(states)

noto_mut_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("noto-mut"))
# Fit a model for each genotype:

noto_models_tbl = noto_mut_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, wt_state_transition_graph, ctrl_ids=c("wt", "ctrl-noto"))) #%>%
#tidyr::unnest(states)


genotype_models_tbl = rbind(noto_models_tbl, crispant_models_tbl)


collect_genotype_effects = function(ccm, timepoint=24, experiment="GAP16"){
  control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=timepoint, experiment=experiment))
  knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=timepoint, experiment=experiment))
  genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
}
#debug(collect_genotype_effects)

# 18 hpf

genotype_models_tbl = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=18)) #%>%
#tidyr::unnest(states)

#plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1) + ggtitle(genotype_models_tbl$gene_target[[1]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[2]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[3]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[4]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[5]], genotype_models_tbl$genotype_eff[[5]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[5]])


genotype_models_tbl = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=24)) #%>%
#tidyr::unnest(states)

#plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1) + ggtitle(genotype_models_tbl$gene_target[[1]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[2]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[3]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[4]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[5]], genotype_models_tbl$genotype_eff[[5]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[5]])


# noto-mut
plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto-mut vs. ctrl 18hpf")
ggsave("noto-mut_pronephros_lfc_18hpf.png", width=4, height=4)


# noto
plot_contrast(genotype_models_tbl$genotype_ccm[[8]], genotype_models_tbl$genotype_eff[[8]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto vs. ctrl 18hpf")
ggsave("noto_pronephros_lfc_18hpf.png", width=4, height=4)


# smo
plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("smo vs. ctrl 18hpf")
ggsave("smo_pronephros_lfc_18hpf.png", width=4, height=4)


# egr2b
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("egr2b vs. ctrl 18hpf")
ggsave("egr2b_pronephros_lfc_18hpf.png", width=4, height=4)



# tbxta
plot_contrast(genotype_models_tbl$genotype_ccm[[6]], genotype_models_tbl$genotype_eff[[6]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("tbxta vs. ctrl 18hpf")
ggsave("tbxta_pronephros_lfc_18hpf.png", width=4, height=4)


# cdx4
plot_contrast(genotype_models_tbl$genotype_ccm[[7]], genotype_models_tbl$genotype_eff[[7]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("cdx4 vs. ctrl 18hpf")
ggsave("cdx4_pronephros_lfc_18hpf.png", width=4, height=4)



# 24 hpf

genotype_models_tbl = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=24)) #%>%
#tidyr::unnest(states)

#plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1) + ggtitle(genotype_models_tbl$gene_target[[1]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[2]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[3]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[4]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[5]], genotype_models_tbl$genotype_eff[[5]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[5]])


# noto-mut
plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto-mut vs. ctrl 24hpf")
ggsave("noto-mut_pronephros_lfc_24hpf.png", width=4, height=4)


# noto
plot_contrast(genotype_models_tbl$genotype_ccm[[8]], genotype_models_tbl$genotype_eff[[8]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto vs. ctrl 24hpf")
ggsave("noto_pronephros_lfc_24hpf.png", width=4, height=4)


# smo
plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("smo vs. ctrl 24hpf")
ggsave("smo_pronephros_lfc_24hpf.png", width=4, height=4)


# egr2b
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("egr2b vs. ctrl 24hpf")
ggsave("egr2b_pronephros_lfc_24hpf.png", width=4, height=4)



# tbxta
plot_contrast(genotype_models_tbl$genotype_ccm[[6]], genotype_models_tbl$genotype_eff[[6]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("tbxta vs. ctrl 24hpf")
ggsave("tbxta_pronephros_lfc_24hpf.png", width=4, height=4)


# cdx4
plot_contrast(genotype_models_tbl$genotype_ccm[[7]], genotype_models_tbl$genotype_eff[[7]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("cdx4 vs. ctrl 24hpf")
ggsave("cdx4_pronephros_lfc_24hpf.png", width=4, height=4)

# 36 hpf

genotype_models_tbl = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=36)) #%>%
#tidyr::unnest(states)

#plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1) + ggtitle(genotype_models_tbl$gene_target[[1]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[2]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[3]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[4]])
#plot_contrast(genotype_models_tbl$genotype_ccm[[5]], genotype_models_tbl$genotype_eff[[5]], q_value_thresh = 0.05) + ggtitle(genotype_models_tbl$gene_target[[5]])


# noto-mut
plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto-mut vs. ctrl 36hpf")
ggsave("noto-mut_pronephros_lfc_36hpf.png", width=4, height=4)

# noto
plot_contrast(genotype_models_tbl$genotype_ccm[[8]], genotype_models_tbl$genotype_eff[[8]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("noto vs. ctrl 36hpf")
ggsave("noto_pronephros_lfc_36hpf.png", width=4, height=4)


# smo
plot_contrast(genotype_models_tbl$genotype_ccm[[4]], genotype_models_tbl$genotype_eff[[4]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("smo vs. ctrl 36hpf")
ggsave("smo_pronephros_lfc_36hpf.png", width=4, height=4)


# egr2b
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("egr2b vs. ctrl 36hpf")
ggsave("egr2b_pronephros_lfc_36hpf.png", width=4, height=4)



# tbxta
plot_contrast(genotype_models_tbl$genotype_ccm[[6]], genotype_models_tbl$genotype_eff[[6]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("tbxta vs. ctrl 36hpf")
ggsave("tbxta_pronephros_lfc_36hpf.png", width=4, height=4)


# cdx4
plot_contrast(genotype_models_tbl$genotype_ccm[[7]], genotype_models_tbl$genotype_eff[[7]], q_value_thresh = 1, plot_edges="none")+
  theme(legend.position="none") +
  ggtitle("cdx4 vs. ctrl 36hpf")
ggsave("cdx4_pronephros_lfc_36hpf.png", width=4, height=4)


## What happens when we run the assembler on the noto crispant model?


ctrl_state_transition_graph = assemble_timeseries_transitions(genotype_models_tbl$genotype_ccm[[8]],
                                                            start=18, stop=48,
                                                            interval_col="timepoint",
                                                            min_interval = 2,
                                                            log_abund_detection_thresh=-2,
                                                            min_dist_vs_time_r_sq=0.0,
                                                            experiment="GAP16",
                                                            knockout=FALSE)

hooke:::plot_path(genotype_models_tbl$genotype_ccm[[8]], path_df = ctrl_state_transition_graph %>% igraph::as_data_frame() , edge_size=0.25)

plot_state_transition_graph(wt_ccm_wl, state_transition_graph %>% igraph::as_data_frame() , color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")



noto_state_transition_graph = assemble_timeseries_transitions(genotype_models_tbl$genotype_ccm[[8]],
                                                              start=18, stop=48,
                                                              interval_col="timepoint",
                                                              min_interval = 2,
                                                              log_abund_detection_thresh=-2,
                                                              min_dist_vs_time_r_sq=0.0,
                                                              experiment="GAP16",
                                                              knockout=TRUE)

hooke:::plot_path(genotype_models_tbl$genotype_ccm[[8]], path_df = noto_state_transition_graph %>% igraph::as_data_frame() , edge_size=0.25)

plot_state_transition_graph(genotype_models_tbl$genotype_ccm[[8]], noto_state_transition_graph %>% igraph::as_data_frame() , color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")


plot_state_transition_graph(genotype_models_tbl$genotype_ccm[[8]], noto_state_transition_graph %>% igraph::as_data_frame() , color_nodes_by = "cell_type", group_nodes_by="cell_type", layer_nodes_by="timepoint")

plot_state_transition_graph(genotype_models_tbl$genotype_ccm[[8]], noto_state_transition_graph %>% igraph::as_data_frame(), cond_b_v_a_tbl=genotype_models_tbl$genotype_eff[[8]], group_nodes_by="cell_type", fc_limits=c(-3,3))

plot_state_transition_graph(genotype_models_tbl$genotype_ccm[[8]], noto_state_transition_graph %>% igraph::as_data_frame(), group_nodes_by="cell_type", genes = c("bcl3"))

###



aberrant_state = "33"
aberrant_state_subgraph_vertices = c(aberrant_state,
                                     igraph::V(noto_state_transition_graph)[get_parents(noto_state_transition_graph, aberrant_state)]$name,
                                     igraph::V(noto_state_transition_graph)[get_children(noto_state_transition_graph, aberrant_state)]$name,
                                     igraph::V(noto_state_transition_graph)[get_siblings(noto_state_transition_graph, aberrant_state)]$name)
aberrant_state_subgraph = igraph::induced_subgraph(noto_state_transition_graph, igraph::V(noto_state_transition_graph)[aberrant_state_subgraph_vertices])

noto_gene_patterns_over_state_graph = classify_genes_over_graph(genotype_models_tbl$genotype_ccm[[8]], aberrant_state_subgraph,  relative_expr_thresh=0.05, abs_expr_thresh=1e-3)
noto_gene_patterns_over_state_graph = noto_gene_patterns_over_state_graph %>% tidyr::unnest(gene_classes) %>% tidyr::unnest(interpretation)
noto_gene_patterns_over_state_graph = left_join(noto_gene_patterns_over_state_graph,
                                           rowData(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds) %>%
                                             as_tibble %>%
                                             select(id, gene_short_name), by=c("gene_id"="id"))

aberrant_state_genes = noto_gene_patterns_over_state_graph %>%
  filter(cell_state == aberrant_state & interpretation %in%  c("Selectively activated", "MLP")) %>%
  pull(gene_short_name) %>% unique
#selective_genes = selective_genes[grepl("si:", selective_genes) == FALSE]
plot_genes_by_group(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds[,clusters(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds) %in% aberrant_state_subgraph_vertices],
                    markers = aberrant_state_genes, group_cells_by = "cluster",ordering_type="maximal_on_diag", color_by_group=TRUE)

aberrant_state_marker_table = top_markers(kidney_cds[rowData(kidney_cds)$gene_short_name %in% aberrant_state_genes,], genes_to_test_per_group = 100)

aberrant_markers = aberrant_state_marker_table %>% filter (cell_group == "33" & marker_test_q_value < 0.01) %>% arrange(desc(marker_score)) %>% pull(gene_short_name)
plot_genes_by_group(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds,
                    markers = aberrant_markers, group_cells_by = "cluster",ordering_type="maximal_on_diag", color_by_group=TRUE)
ggsave("aberrant_state_markers.png", width=12, height=6)

aberrant_state_gene_tfs = noto_gene_patterns_over_state_graph %>%
  filter(cell_state == aberrant_state & interpretation %in%  c("Selectively activated", "MLP") & gene_short_name %in% transcription_regulators) %>%
  pull(gene_short_name) %>% unique
#selective_genes = selective_genes[grepl("si:", selective_genes) == FALSE]
plot_genes_by_group(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds[,clusters(genotype_models_tbl$genotype_ccm[[8]]@ccs@cds) %in% aberrant_state_subgraph_vertices],
                    markers = aberrant_state_gene_tfs, group_cells_by = "cluster",ordering_type="maximal_on_diag", color_by_group=TRUE)

# Plot the aberrant transcription factors:
plot_cells(kidney_cds, genes = aberrant_state_gene_tfs)

#gene_universe = gene_patterns_over_cell_type_graph %>% filter(interpretation != "Absent") %>% pull(gene_short_name) %>% unique()
aberrant_state_pattern_pathways = noto_gene_patterns_over_state_graph %>%
  filter(cell_state == aberrant_state) %>%
  select(cell_state, interpretation, gene_short_name) %>%
  group_by(cell_state, interpretation) %>% tidyr::nest(data=gene_short_name) %>%
  dplyr::mutate(pathways = purrr::map(.f = purrr::possibly(
    calc_pathway_enrichment_on_state_specific_genes, NA_real_), .x = data, gene_set_df, sig_thresh=1)) %>%
  tidyr::unnest(pathways)
#tidyr::unnest(degs)

plot_cells(kidney_cds, color_cells_by="aberant.noto.smo.cells") + scale_color_manual(values=c("lightgrey", "red")) + theme(legend.position="none")
ggsave("highlight_aberrant_state.png", width=4, height=4)

plot_cells(kidney_cds, color_cells_by="aberant.noto.smo.cells", 2, 3) + scale_color_manual(values=c("lightgrey", "red")) + theme(legend.position="none")
ggsave("highlight_aberrant_state_alt_view.png", width=4, height=4)


plot_cells(kidney_cds, genes = c("bcl3","nfkbie", "noxa1", "nox1", "ltb4r2a", "gfpt2"), 2, 3)
ggsave("nfkb_genes_in_aberrant_state.png", width=10, height=6)

plot_cells(kidney_cds, genes = c("aqp3a",  "evx1", "dnase1l4.1"), 2, 3)
ggsave("cloaca_markers_in_aberrant_state.png", width=10, height=3)

plot_cells(kidney_cds, genes = c("twist", "snail", "slug", "chd1", "cdh2"), 2, 3)
ggsave("emt_genes_in_aberrant_state.png", width=10, height=6)


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


