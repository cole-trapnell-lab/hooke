library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
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

set.seed(42)

kidney_cds = cluster_cells(kidney_cds, resolution=1e-3)
colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)

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
                                         marker_file = "./examples/pronephros_cell_types.txt", #marker_file_path,
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
  ggsave("kidney_segment.png", width=3, height=3)


plot_cells(kidney_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_cell_types.png", width=3, height=3)


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
  ggsave("kidney_time.png", width=3, height=3)


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


wt_ccm  = new_cell_count_model(wt_ccs,
                               main_model_formula_str = "~splines::ns(timepoint,df=4)",
                               nuisance_model_formula_str = "~experiment")


timepoint_pred_df = tibble(timepoint=seq(18,48))

timepoint_pred_df = timepoint_pred_df %>%
  dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
    function(tp){ estimate_abundances(wt_ccm, tibble(timepoint=tp, experiment="GAP14"))}, NA_real_), .x = timepoint)) %>% tidyr::unnest()
qplot(timepoint, log_abund, data=timepoint_pred_df, geom="line", color=cell_group) +
  geom_ribbon(aes(timepoint, ymin=log_abund-2*log_abund_se, ymax=log_abund+2*log_abund_se), alpha=I(0.1)) +
  facet_wrap(~cell_group, scales="free_y")


plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                         timepoint_pred_df %>% filter(timepoint == 24),
                                         timepoint_pred_df %>% filter(timepoint == 30)),
              q_value_thresh = 0.01)


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
  return(genotype_ccm)
}
undebug(fit_genotype_ccm)


genotype_df = colData(ccs@cds) %>% as_tibble() %>% dplyr::select(gene_target) %>% distinct()

crispant_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("smo", "noto"))
# Fit a model for each genotype:

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
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=24)) #%>%
#tidyr::unnest(states)

plot_contrast(genotype_models_tbl$genotype_ccm[[1]], genotype_models_tbl$genotype_eff[[1]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[2]], genotype_models_tbl$genotype_eff[[2]], q_value_thresh = 0.05)
plot_contrast(genotype_models_tbl$genotype_ccm[[3]], genotype_models_tbl$genotype_eff[[3]], q_value_thresh = 0.05)

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



### Old stuff:



noto_smo_cds = kidney_cds[,colData(kidney_cds)$genotype %in% c("wt", "noto", "noto-mut", "smo")]


noto_smo_gap16_cds = kidney_cds[,colData(kidney_cds)$genotype %in% c("wt", "noto", "smo", "noto-mut") &
                                  colData(kidney_cds)$experiment == "GAP16" &
                                  colData(kidney_cds)$timepoint %in% c(18,24,36)]

gap16_cds = kidney_cds[,colData(kidney_cds)$experiment == "GAP16"]



plot_cells(noto_smo_gap16_cds, color_cells_by="timepoint.1", show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_time.png", width=3, height=3)



#subset_kidney_cds = kidney_cds[,colData(kidney_cds)$genotype]

ccs = new_cell_count_set(kidney_cds,
                         sample_group = "Oligo",
                         cell_group = "cell_type")


ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~splines::ns(timepoint,df=3) + genotype",
                            nuisance_model_formula_str = "experiment")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=0.1)


cond_wt_18 = estimate_abundances(ccm, tibble::tibble(timepoint=18, genotype="wt", experiment="GAP16"))
cond_wt_24 = estimate_abundances(ccm, tibble::tibble(timepoint=24, genotype="wt", experiment="GAP16"))
cond_wt_36 = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="wt", experiment="GAP16"))
cond_18_vs_24_wt = compare_abundances(ccm, cond_wt_18, cond_wt_24)
plot_contrast(ccm, cond_18_vs_24_wt, scale_shifts_by="none", q_value_thresh=0.01)

cond_24_vs_36_wt = compare_abundances(ccm, cond_wt_24, cond_wt_36)
plot_contrast(ccm, cond_24_vs_36_wt, scale_shifts_by="none", q_value_thresh=0.01)


cond_smo_18 = estimate_abundances(ccm, tibble::tibble(timepoint=18, genotype="smo", experiment="GAP16"))

cond_smo = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="smo", experiment="GAP16"))
cond_wt = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="wt", experiment="GAP16"))
cond_noto = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="noto", experiment="GAP16"))
cond_noto_mut = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="noto-mut", experiment="GAP16"))

#cond_tbxta = estimate_abundances(ccm, tibble::tibble(timepoint.1="36", genotype="tbxta", experiment="GAP16"))

cond_met = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="met", experiment="GAP16"))


# cond_ctrl_vs_wt = compare_abundances(ccm, cond_wt, cond_ctrl_inj)
#
# plot_contrast(ccm, cond_ctrl_vs_wt, scale_shifts_by="none", p_value_thresh=0.05, plot_labels="none") +
#   theme_minimal() + monocle3:::monocle_theme_opts() +
#   ggsave("kidney_wt_vs_ctrl_inj.png", width=7, height=6)

cond_smo_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_smo)

plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01)


plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01, plot_labels="none", receiver_cell_groups=c("none")) +
  ggsave("smo_vs_wt_36h_fc.png", width=4.5, height=3)


#plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", p_value_thresh=1, receiver_cell_groups = c("Neck Segment"))

cond_noto_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_noto)

plot_contrast(ccm, cond_noto_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01)


plot_contrast(ccm, cond_noto_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01, plot_labels="none", receiver_cell_groups=c("none"))+
  ggsave("noto_vs_wt_36h_fc.png", width=4.5, height=3)


