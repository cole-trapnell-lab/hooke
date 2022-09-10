library(tidyverse)
library(hooke)
library(garnett)
library(msigdbr)
library(fgsea)
library(googledrive)

plot_cell_covariances <- function(ccm,
                                  x=1,
                                  y=2,
                                  edge_size=2,
                                  cell_size=1,
                                  group_label_size=2,
                                  plot_labels = c("all", "none"),
                                  cell_groups = NULL,
                                  pcor_limits=c(-1,1)){

  umap_centers = centroids(ccm@ccs)


  #cond_b_vs_a_tbl$delta_q_value = p.adjust(cond_b_vs_a_tbl$delta_q_value, method = "BH")
  # umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

  cov_graph = hooke:::return_igraph(model(ccm))
  corr_edge_coords = igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight != 0.00)
  corr_edge_coords = corr_edge_coords %>% dplyr::select(from, to, weight) %>% dplyr::rename(pcor = weight)
  corr_edge_coords = dplyr::left_join(corr_edge_coords, umap_centers %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group"))
  corr_edge_coords = dplyr::left_join(corr_edge_coords, umap_centers %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))

  undirected_edge_df = corr_edge_coords

  if (is.null(cell_groups) == FALSE){
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(from %in% cell_groups | to %in% cell_groups)
  }

  # corr_edge_coords_umap_delta_abund = left_join(corr_edge_coords_umap_delta_abund,
  #                                               umap_centers,
  #                                               by=c("from"="cell_group"))

  undirected_edge_df = undirected_edge_df %>%
    dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))

  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)

  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,y]

  # plot
  gp = ggplot() +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2),
      color = "black",
      size = 1.5 * cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2),
      color = "white",
      size = cell_size,
      stroke = 0
    ) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    monocle3:::monocle_theme_opts()

  gp = gp  +
    geom_segment(data = undirected_edge_df,
                 aes(x = get(paste0("to_umap_", x)),
                     y = get(paste0("to_umap_", y)),
                     color=pcor,
                     xend = get(paste0("from_umap_", x)),
                     yend = get(paste0("from_umap_", y)),
                     size = edge_size * scaled_weight)) +
    #size=edge_size / 4) +
    scale_color_gradient2(
      low = "#122985",
      mid = "white",
      high = "red4",
      na.value = "white",
      limits = pcor_limits
    )  +
    scale_size_identity()

  if (plot_labels != "none") {
    label_df = umap_centers
    gp <- gp + ggrepel::geom_label_repel(data = label_df,
                                         mapping = aes(umap_1, umap_2, label=cell_group),
                                         size=I(group_label_size),
                                         fill = "white")
  }
  return(gp)
}
#debug(plot_cell_covariances)

drive_download("https://drive.google.com/file/d/15TIyrLLVylsJCiiZqufxIJ1MuImGDz2o/")
a549_ccs = readRDS("a549_ccs.rds")

a549_cds = a549_ccs@cds

set.seed(42)

a549_cds = cluster_cells(a549_cds, resolution=1e-4, random_seed=42)
plot_cells(a549_cds)

colData(a549_cds)$sample = NULL
colData(a549_cds)$log_dose = log10(colData(a549_cds)$dose+1)

colData(a549_cds)$grouping = paste(colData(a549_cds)$catalog_number,
                                   colData(a549_cds)$dose,
                                   colData(a549_cds)$replicate,
                                   colData(a549_cds)$top_oligo_W, sep=".")

colData(a549_cds)$cluster = monocle3::clusters(a549_cds)



state_markers = c("CDK1",  #High in G2M, Low in M-G1 transition (https://elifesciences.org/articles/71356)
                  "TOP2A", #High in G2M, Low in M-G1 transition (https://elifesciences.org/articles/71356)
                  "UBE2C", #High in G2M, Low in M-G1 transition (https://elifesciences.org/articles/71356)
                  "FBXO5", #High in G2M, Low in M-G1 transition (https://elifesciences.org/articles/71356)
                  "FZR1", #High in G2M, Low in M-G1 transition (https://elifesciences.org/articles/71356)
                  "MKI67",
                  "CENPA", #Expressed in both G1M and M-G1 transition (https://elifesciences.org/articles/71356)
                  "PSD3",  #Expressed in both G1M and M-G1 transition (https://elifesciences.org/articles/71356)
                  "SRGAP1", #Expressed in both G1M and M-G1 transition (https://elifesciences.org/articles/71356)
                  "HIST1H2AE",
                  "CCNB1", # G2M
                  "CCNB2", # G2M
                  "CCNA2", # G2M
                  "PCNA", #S phase > G2M (Seurat)
                  "MCM6", #S phase > G2M (Seurat)
                  "ACSS2", # acetate starvation
                  "ACLY",  # acetate starvation
                  "CDKN1A",# p53 DNA damage & arrest. Very high in p53 apoptosis
                  "BAX", # High in p53 apoptosis,
                  "GATA3" # contaminant from other cell lines (e.g. cryptic doublet)
)

plot_cells(a549_cds,
           color_cells_by="viability",
           scale_to_range=TRUE,
           show_trajectory_graph = FALSE)
ggsave("viability.png", width=12, height=12)

plot_cells(a549_cds,
           genes=state_markers,
           scale_to_range=TRUE,
           show_trajectory_graph = FALSE)
ggsave("state_markers.png", width=12, height=12)

plot_genes_by_group(a549_cds, markers=state_markers)
ggsave("state_markers_dotplot.png", width=12, height=12)

# cell_group_df <- tibble::tibble(cell=row.names(colData(a549_cds)),
#                                 cell_group=clusters(a549_cds)[colnames(a549_cds)])
# marker_expr_in_states = aggregate_gene_expression(a549_cds[rowData(a549_cds)$gene_short_name %in% state_markers],
#                                                   cell_group_df=cell_group_df,
#                                                   scale_agg_values=FALSE)
# row.names(marker_expr_in_states) = rowData(a549_cds)[row.names(marker_expr_in_states),]$gene_short_name
# pheatmap::pheatmap(marker_expr_in_states)

# xxx = state_transition_pathways %>% filter(state == 38) %>% pull(pathways)
# ribo_genes = xxx[[1]] %>% filter(padj < 0.01 & pathway == "GO_RIBOSOME_BIOGENESIS") %>% pull(leadingEdge) %>% unlist()
# marker_expr_in_states = aggregate_gene_expression(a549_cds[rowData(a549_cds)$gene_short_name %in% ribo_genes], cell_group_df=cell_group_df)
# row.names(marker_expr_in_states) = rowData(a549_cds)[row.names(marker_expr_in_states),]$gene_short_name
# pheatmap::pheatmap(marker_expr_in_states, fontsize=6)


# See Seurat docs and this paper: https://elifesciences.org/articles/71356
colData(a549_cds)$cell_state_annotation = case_when(
  colData(a549_cds)$cluster %in% c(29) ~ "M-G1 transition",
  colData(a549_cds)$cluster %in% c(11,15) ~ "S",
  colData(a549_cds)$cluster %in% c(23,28,12,5) ~ "G2M",
  colData(a549_cds)$cluster %in% c(17, 21, 27, 10, 22, 3) ~ "G1",
  colData(a549_cds)$cluster %in% c(14) ~ "G0",
  colData(a549_cds)$cluster %in% c(33) ~ "p53 apoptosis", #
  colData(a549_cds)$cluster %in% c(2,4,25,19,20,26) ~ "p53-dependent arrest", #
  colData(a549_cds)$cluster %in% c(8,30) ~ "p53-independent arrest", #
  colData(a549_cds)$cluster %in% c(32) ~ "Apoptosis", # Annotated by looking at the viability plot
  colData(a549_cds)$cluster %in% c(31) ~ "Acetate starvation",
  colData(a549_cds)$cluster %in% c(34) ~ "Glucocorticoid response",
  colData(a549_cds)$cluster %in% c(35) ~ "Contaminant",
  #colData(kidney_cds)$cluster == 16 ~ "Unknown",
  TRUE ~ "Unknown"
)
plot_cells(a549_cds, color_cells_by="cell_state_annotation", show_trajectory_graph = FALSE)

sciplex_classifier <- train_cell_classifier(cds = a549_cds,
                                            marker_file = "./examples/a549_cell_states.txt",
                                            db="none",#org.Dr.eg.db::org.Dr.eg.db,
                                            cds_gene_id_type = "SYMBOL",
                                            num_unknown = 50,
                                            marker_file_gene_id_type = "SYMBOL")
#saveRDS(sciplex_classifier, "sciplex_classifier.rds")

a549_cds = classify_cells(a549_cds, sciplex_classifier, db="none")
prop.table(table(colData(a549_cds[,colData(a549_cds)$catalog_number == "S0000" & colData(a549_cds)$cell_type != "Unknown"])$cell_type))

plot_cells(a549_cds, color_cells_by="cell_type", show_trajectory_graph = FALSE)

a549_cds = a549_cds[,colData(a549_cds)$cell_type != "Contaminant"]

# FIXME: Weird thing that should be maybe fixed?
# Once we drop the contaminant cells we need to re-cluster in order for assembly to work. We
# Can't have cluster levels that don't exist anymore in CDS.

a549_cds = cluster_cells(a549_cds, resolution=1e-4, random_seed=42)
colData(a549_cds)$cluster = monocle3::clusters(a549_cds)

a549_ccs = new_cell_count_set(a549_cds,
                              sample_group = "grouping",
                              cell_group = "cluster")

a549_ccs_coldata = colData(a549_ccs) %>% as.data.frame()

product_names = unique(colData(a549_cds)$product_name)
catalog_numbers = unique(colData(a549_cds)$catalog_number)

for (cn in catalog_numbers) {
  if (cn != "S0000") {
    a549_ccs_coldata = a549_ccs_coldata %>%
      mutate("{cn}":= ifelse(catalog_number==cn, log_dose, 0))
    colData(a549_ccs)[[cn]] = a549_ccs_coldata[[cn]]
  }
}


compound_df = colData(a549_ccs@cds) %>% as_tibble() %>% dplyr::select(catalog_number, pathway_level_1, pathway_level_2, product_name, target) %>% distinct()

#test_compound_df = dplyr::filter(compound_df, catalog_number %in% c("S2693", "S1192"))
# Fit a model for each drug:

fit_drug_ccm = function(compound, ccs, vehicle_compound_id="S0000", num_dose_breaks=3){
  subset_ccs = ccs[,colData(ccs)[,compound] > 0 | colData(ccs)$treatment == vehicle_compound_id]

  colData(subset_ccs)$compound_dose = colData(subset_ccs)[,compound]
  # FIXME: Should prob use splines
  dose_start = min(colData(subset_ccs)$compound_dose)
  dose_stop = max(colData(subset_ccs)$compound_dose)

  dose_breakpoints = c()
  if (num_dose_breaks > 2 & dose_stop > dose_start){
    dose_breakpoints = seq(dose_start, dose_stop, length.out=num_dose_breaks)
    dose_breakpoints = dose_breakpoints[2:(length(dose_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
    main_model_formula_str = paste("~ splines::ns(compound_dose, knots=", paste("c(",paste(dose_breakpoints, collapse=","), ")", sep=""), ")")
  }else{
    main_model_formula_str = "~ 1"
  }
  compound_ccm = new_cell_count_model(subset_ccs,
                                      main_model_formula_str = main_model_formula_str)
  #compound_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
  #                                    main_model_formula_str = "~ compound_dose"))
  return(compound_ccm)
}
#undebug(fit_drug_ccm)

compound_models_tbl = compound_df %>%
  dplyr::mutate(drug_ccm = purrr::map(.f = purrr::possibly(
    fit_drug_ccm, NA_real_), .x = catalog_number, a549_ccs)) #%>%
#tidyr::unnest(states)

# collect_dose_dep_effects = function(ccm, low_dose, high_dose){
#   low_dose_abund = estimate_abundances(ccm, tibble(compound_dose = low_dose))
#   high_dose_abund = estimate_abundances(ccm, tibble(compound_dose = high_dose))
#   dose_comparison_tbl = compare_abundances(ccm, low_dose_abund, high_dose_abund)
# }
# compound_models_tbl = compound_models_tbl %>%
#   dplyr::mutate(dose_eff = purrr::map(.f = purrr::possibly(
#     collect_dose_dep_effects, NA_real_), .x = drug_ccm, low_dose=0, high_dose=3)) #%>%
# #tidyr::unnest(states)

#############
# Test out assembly of a state transition graph over the doses:

#drug_ccm = compound_models_tbl %>% filter(product_name == "Pracinostat (SB939)") %>% pull(drug_ccm)
drug_ccm = compound_models_tbl %>% filter(grepl("Aurora A Inhibitor I", product_name)) %>% pull(drug_ccm)
drug_ccm = drug_ccm[[1]]
drug_ccm = select_model(drug_ccm, criterion = "EBIC", sparsity_factor=0.01)


drug_state_transition_graph = assemble_timeseries_transitions(drug_ccm,
                                                            start=0, stop=4,
                                                            interval_col="compound_dose",
                                                            min_interval = 0.5,
                                                            log_abund_detection_thresh=-2,
                                                            min_dist_vs_time_r_sq=0.0)


hooke:::plot_path(drug_ccm, path_df = drug_state_transition_graph %>% igraph::as_data_frame(), edge_size=1)

plot_state_transition_graph(drug_ccm, drug_state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_state_annotation", group_nodes_by="cell_state_annotation")

#############

collect_dose_dep_effects = function(ccm, vehicle_dose=0, max_dose=4, dose_step=0.25, tcX_sig_thresh=0.1, tcX=0.5){
  #vehicle_abundances = estimate_abundances(ccm, tibble(compound_dose = vehicle_dose))

  dose_pred_df = tibble(dose=seq(vehicle_dose, max_dose, by=dose_step))
  dose_pred_df = dose_pred_df %>%
    dplyr::mutate(dose_eff = purrr::map(.f = purrr::possibly(
      function(dose){
        vehicle_abund = estimate_abundances(ccm, tibble(compound_dose=0))
        dose_abund = estimate_abundances(ccm, tibble(compound_dose=dose))
        dose_comparison_tbl = compare_abundances(ccm, vehicle_abund, dose_abund)
        return(dose_comparison_tbl)
      }, NA_real_), .x = dose)) %>% tidyr::unnest()
  if (is.null(tcX)){
   max_dose_df = dose_pred_df %>% filter(compound_dose_y == max_dose)
   return(max_dose_df)
  }else{
    sig_doses = dose_pred_df %>% filter(delta_q_value < tcX_sig_thresh & delta_log_abund < log(tcX))
    if (nrow(sig_doses) > 0) {
      tcX_dose = sig_doses %>% arrange (compound_dose_y) %>% slice_head(n=1) %>% pull(compound_dose_y)
      tcX_dose = dose_pred_df %>% filter(compound_dose_y == tcX_dose)
      return(tcX_dose)
    } else {
      return(NA_real_)
    }
  }

  #return(dose_pred_df)
}
#debug(collect_dose_dep_effects)

compound_models_tbl = compound_models_tbl %>%
  dplyr::mutate(tc50_eff = purrr::map(.f = purrr::possibly(
                            collect_dose_dep_effects, NA_real_), .x = drug_ccm, tcX = 0.5),
                tc90_eff = purrr::map(.f = purrr::possibly(
                            collect_dose_dep_effects, NA_real_), .x = drug_ccm, tcX = 0.1),
                max_dose_eff = purrr::map(.f = purrr::possibly(
                  collect_dose_dep_effects, NA_real_), .x = drug_ccm, tcX = NULL)
    ) #%>%
#tidyr::unnest(states)




#state_transitions = hooke:::collect_pln_graph_edges(ccm, cond_ra_vs_wt_tbl) %>% as_tibble %>%
#  filter(edge_type == "directed_from_to" & to_delta_p_value < 0.05 & from_delta_p_value < 0.05)


collect_state_transitions = function(ccm, dose_contrast){
  transitions = hooke:::collect_pln_graph_edges(ccm, dose_contrast)
  transitions = transitions %>% dplyr::filter(edge_type != "undirected")
  return(transitions)
}

compound_models_tbl = compound_models_tbl %>%
  dplyr::mutate(state_transitions = purrr::map2(.f = purrr::possibly(
    collect_state_transitions, NA_real_), .x = drug_ccm,
    .y = max_dose_eff)) #%>%
#tidyr::unnest(states)

state_transition_effects = compound_models_tbl %>%
  dplyr::select(catalog_number, pathway_level_1, pathway_level_2, product_name, target, state_transitions) %>%
  tidyr::unnest(state_transitions)

sig_state_transition_effects = state_transition_effects %>%
  mutate(from_delta_q_value = p.adjust(from_delta_q_value, method="BH"),
         to_delta_q_value = p.adjust(to_delta_q_value, method="BH")) %>%
  dplyr::filter(from_delta_q_value < 0.01 & to_delta_q_value < 0.01)

qplot(pathway_level_1, data=sig_state_transition_effects) + coord_flip() + facet_grid(from~to) + monocle3:::monocle_theme_opts()


# Some contrast plots with clear expectations:
hdaci_df = compound_models_tbl %>% filter(grepl("Pracinostat", product_name))
plot_contrast(hdaci_df$drug_ccm[[1]], hdaci_df$max_dose_eff[[1]], q_value_thresh = 0.1)

gr_agonist_df = compound_models_tbl %>% filter(grepl("Triamcinolone", product_name))
plot_contrast(gr_agonist_df$drug_ccm[[1]], gr_agonist_df$max_dose_eff[[1]], q_value_thresh = 0.1)

trametinib_df = compound_models_tbl %>% filter(grepl("Trametinib", product_name))
plot_contrast(trametinib_df$drug_ccm[[1]], trametinib_df$max_dose_eff[[1]], q_value_thresh = 0.1)

epothilone_df = compound_models_tbl %>% filter(grepl("Epothilone A", product_name))
plot_contrast(epothilone_df$drug_ccm[[1]], epothilone_df$tc50_eff[[1]], q_value_thresh = 0.1)
plot_contrast(epothilone_df$drug_ccm[[1]], epothilone_df$max_dose_eff[[1]], q_value_thresh = 0.1)

nuc_analog_df = compound_models_tbl %>% filter(grepl("5-FU", product_name))
plot_contrast(nuc_analog_df$drug_ccm[[1]], nuc_analog_df$max_dose_eff[[1]], q_value_thresh = 0.1)

plot_contrast(nuc_analog_df$drug_ccm[[1]],
              compare_abundances(nuc_analog_df$drug_ccm[[1]],
                                 estimate_abundances(nuc_analog_df$drug_ccm[[1]], tibble(compound_dose = 0)),
                                 estimate_abundances(nuc_analog_df$drug_ccm[[1]], tibble(compound_dose = 1))),
              q_value_thresh = 0.1)

plot_contrast(nuc_analog_df$drug_ccm[[1]],
              compare_abundances(nuc_analog_df$drug_ccm[[1]],
                                 estimate_abundances(nuc_analog_df$drug_ccm[[1]], tibble(compound_dose = 2)),
                                 estimate_abundances(nuc_analog_df$drug_ccm[[1]], tibble(compound_dose = 3))),
              q_value_thresh = 0.1)


topoisomerase_inh_df = compound_models_tbl %>% filter(grepl("Pirarubicin", product_name))
plot_contrast(topoisomerase_inh_df$drug_ccm[[1]], topoisomerase_inh_df$max_dose_eff[[1]], q_value_thresh = 0.1)


plot_contrast(topoisomerase_inh_df$drug_ccm[[1]],
              compare_abundances(topoisomerase_inh_df$drug_ccm[[1]],
                                 estimate_abundances(topoisomerase_inh_df$drug_ccm[[1]], tibble(compound_dose = 3)),
                                 estimate_abundances(topoisomerase_inh_df$drug_ccm[[1]], tibble(compound_dose = 4))),
              q_value_thresh = 0.1)

### Characterizing differences between the molecular states:
# Generate a pseudobulk expression matrix
# Note: normally we'd do this using a549_ccs@metadata[["cell_group_assignments"]]
# but I'm going to collapse more aggressively for the sake of speed during testing:
cell_group_df = colData(a549_ccs@cds) %>% as.data.frame %>%
   dplyr::select(cell, culture_plate, cluster) %>%
   group_by(culture_plate, cluster)
cell_group_df$group_id = group_indices(cell_group_df)
cell_group_df = cell_group_df %>% ungroup %>% dplyr::select(cell, group_id, cell_group=cluster)
row.names(cell_group_df) = cell_group_df$cell

agg_expr_mat = monocle3::aggregate_gene_expression(a549_ccs@cds,
                                                   #cell_group_df = tibble::rownames_to_column(a549_ccs@metadata[["cell_group_assignments"]]),
                                                   cell_group_df = cell_group_df,
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0,
                                                   cell_agg_fun="mean")


agg_coldata = cell_group_df %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
row.names(agg_coldata) = colnames(agg_expr_mat)

pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(a549_ccs@cds) %>% as.data.frame)
#pseudobulk_cds = pseudobulk_cds[rowData(a549_ccs@cds)$num_cells_expressed > 100,
#                                colData(pseudobulk_cds)$num_cells_in_group > 50]
pseudobulk_cds = pseudobulk_cds[rowData(a549_ccs@cds)$num_cells_expressed > 100,]
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
pseudobulk_cds = preprocess_cds(pseudobulk_cds)
pseudobulk_cds = reduce_dimension(pseudobulk_cds)



# # This function compares two cell states to find genes that differ between them
# find_degs_between_states = function(from_state,
#                                     to_state,
#                                     cds,
#                                     model_formula_str = "~cell_group",
#                                     gene_whitelist = NULL,
#                                     q_value_thresh=0.05,
#                                     effect_thresh = 2,
#                                     cores=1){
#   cds_states_tested = cds[,as.character(colData(cds)$cell_group) %in% c(from_state, to_state)]
#   models = monocle3::fit_models(cds_states_tested,
#                                 model_formula_str=model_formula_str,
#                                 weights=colData(cds_states_tested)$num_cells_in_group,
#                                 cores=cores)
#   coefs = monocle3::coefficient_table(models) %>%
#     #filter(grepl("cell_group", term) & q_value < q_value_thresh & abs(normalized_effect) > effect_thresh)
#     filter(grepl("cell_group", term)) %>% dplyr::select(id, gene_short_name, estimate, p_value)
#   rm(models)
#   gc()
#   return(coefs)
# }
# #debug(find_degs_between_states)
#
state_transitions = state_transition_effects %>%
  dplyr::filter(from_delta_q_value < 0.01 & to_delta_q_value < 0.01) %>%
  dplyr::select(from, to) %>%
  dplyr::distinct()
#
#
# # For each reciprocal state transition, find DEGs between the endpoints
# state_transition_degs = state_transitions %>% #head(3) %>%
#   dplyr::mutate(degs = purrr::map2(.f = purrr::possibly(
#     find_degs_between_states, NA_real_), .x = from, .y = to, pseudobulk_cds, cores=4)) #%>%
# #tidyr::unnest(degs)

# This function compares two cell states to find genes that differ between them
compute_gene_fcs_between_states = function(from_state,
                                    to_state,
                                    cds,
                                    gene_whitelist = NULL,
                                    max_abs_effect_size = 10,
                                    cores=1){
  from_cells = cds[,as.character(colData(cds)$cell_group) %in% c(from_state)]
  #from_cells = detect_genes(from_cells)

  to_cells = cds[,as.character(colData(cds)$cell_group) %in% c(to_state)]
  #to_cells = detect_genes(to_cells)

  from_cells_mean_expr = Matrix::rowMeans(normalized_counts(from_cells))
  from_cells_mean_expr = from_cells_mean_expr *  (colData(from_cells)$num_cells_in_group / sum(colData(from_cells)$num_cells_in_group))

  to_cells_mean_expr = Matrix::rowMeans(normalized_counts(to_cells))
  to_cells_mean_expr = to_cells_mean_expr *  (colData(to_cells)$num_cells_in_group / sum(colData(to_cells)$num_cells_in_group))

  log2_fcs = log2(to_cells_mean_expr / from_cells_mean_expr)
  log2_fcs[is.nan(log2_fcs)] = 0
  log2_fcs[is.na(log2_fcs)] = 0
  log2_fcs[log2_fcs > max_abs_effect_size] = max_abs_effect_size
   log2_fcs[log2_fcs < -max_abs_effect_size] = -max_abs_effect_size
  fcs_df = tibble(id=names(log2_fcs),
                  gene_short_name=rowData(cds)[names(log2_fcs),]$gene_short_name,
                  estimate=log2_fcs,
                  p_value=1)

  # coefs = monocle3::coefficient_table(models) %>%
  #   #filter(grepl("cell_group", term) & q_value < q_value_thresh & abs(normalized_effect) > effect_thresh)
  #   filter(grepl("cell_group", term)) %>% dplyr::select(id, gene_short_name, estimate, p_value)

  return(fcs_df)
}
#debug(compute_gene_fcs_between_states)

state_transition_fcs = state_transitions %>% #head(3) %>%
  dplyr::mutate(gene_fcs = purrr::map2(.f = purrr::possibly(
    compute_gene_fcs_between_states, NA_real_), .x = from, .y = to, pseudobulk_cds, cores=4))

transition_pathway_fcs = state_transition_fcs %>%
  dplyr::select(from, to, gene_fcs) %>% unnest() %>%
  dplyr::select(from, to, gene_short_name, estimate) %>% ungroup()
state_entry_gene_scores = transition_pathway_fcs %>%
  group_by(gene_short_name, to) %>% summarize(entry_score = mean(estimate)) %>% ungroup() %>% mutate(state=to) %>%
  dplyr::select(state, gene_short_name, entry_score)
state_exit_gene_scores = transition_pathway_fcs %>%
  group_by(gene_short_name, from) %>% summarize(exit_score = mean(estimate)) %>% ungroup() %>% mutate(state=from)  %>%
  dplyr::select(state, gene_short_name, exit_score)
state_gene_scores = full_join(state_entry_gene_scores, state_exit_gene_scores) %>%
  mutate_if(is.numeric,coalesce,0)
state_gene_scores = state_gene_scores %>% mutate(gene_state_score = entry_score - exit_score)

#### NEW WAY:
#gene_set = msigdbr(species = "Homo sapiens", subcategory = "GO:MF")
#gene_set = msigdbr(species = "Homo sapiens", category = "H")
gene_set = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")
gene_set_list = split(x = gene_set$gene_symbol, f = gene_set$gs_name)
gene_set_df = gene_set %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

calc_pathway_enrichment_on_transition_effects = function(gene_lfc_tbl, gene_set_list, sig_thresh = 0.1, ...){
  #gene_set_list = split(x = pathways$gene_symbol, f = pathways$gs_name)
  gene_ranking = gene_lfc_tbl %>% pull(gene_state_score)
  names(gene_ranking) = gene_lfc_tbl %>% pull(gene_short_name)
  #pb$tick()
  gsea_res = fgsea(pathways=gene_set_list, stats=gene_ranking, ...) %>% as_tibble()
  gsea_res = gsea_res %>% filter(padj < sig_thresh)
  gc()
  return(gsea_res)
}
#debug(calc_pathway_enrichment_on_transition_effects)

state_transition_pathways = state_gene_scores %>%
  group_by(state) %>% nest() %>%
  dplyr::mutate(pathways = purrr::map(.f = purrr::possibly(
    calc_pathway_enrichment_on_transition_effects, NA_real_), .x = data, gene_set_list, sig_thresh=1)) #%>%
#tidyr::unnest(degs)

state_pathway_score_mat = state_transition_pathways %>% dplyr::select(state, pathways) %>% unnest() %>%
  dplyr::select(state, pathway, padj, NES) #%>%
  #dplyr::filter(padj < 0.01) %>%
sig_pathway_ids = state_pathway_score_mat %>% dplyr::filter(padj < 0.01 & abs(NES) > 1) %>% pull(pathway) %>% unique
state_pathway_score_mat = state_pathway_score_mat %>% filter(pathway %in% sig_pathway_ids) %>%
  dplyr::select(state, pathway, NES) %>%
  pivot_wider(names_from=state, values_from=NES, values_fill=0) %>% as.data.frame
row.names(state_pathway_score_mat) = state_pathway_score_mat$pathway
state_pathway_score_mat = state_pathway_score_mat[,-1]
state_pathway_score_mat = as.matrix(state_pathway_score_mat)

pheatmap::pheatmap(state_pathway_score_mat, show_rownames=T, scale="row", clustering_method="ward.D2", fontsize=6)

## We want this stuff below, but plot_cells() needs to work to make it look nice and sensible
#
# selected_pathways = c("GO_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION",
#                       "GO_CELL_CYCLE_G1_S_PHASE_TRANSITION",
#                       "GO_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
#                       "GO_CELL_CYCLE_G2_M_PHASE_TRANSITION",
#                       "GO_MITOTIC_G2_M_TRANSITION_CHECKPOINT",
#                       "GO_DNA_REPLICATION_CHECKPOINT",
#                       "GO_DNA_INTEGRITY_CHECKPOINT",
#                       "GO_RIBOSOME_BIOGENESIS",
#                       "GO_TRANSLATIONAL_ELONGATION",
#                       "GO_APOPTOTIC_SIGNALING_PATHWAY",
#                       "GO_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS",
#                       "GO_APOPTOTIC_PROCESS",
#                       "GO_INFLAMMATORY_RESPONSE",
#                       "GO_CELLULAR_RESPIRATION"
#                       )
# selected_pathway_genes_df = gene_set_df %>% filter(gs_name %in% selected_pathways)
# selected_pathway_genes_df =
#   left_join(selected_pathway_genes_df,
#             rowData(a549_cds) %>% as.data.frame %>% dplyr::select(id, gene_short_name),
#             by=c("gene_symbol"="gene_short_name")) %>% dplyr::select(id, gs_name) %>% filter(is.na(id) == FALSE)
# plot_cells(a549_cds, genes=selected_pathway_genes_df, show_trajectory_graph=F, scale_to_range=TRUE, rasterize=TRUE)
# # pheatmap::pheatmap(pathway_expr_mat[selected_pathways,], scale="row", clustering_method="ward.D2")

# This is totally unhelpful (and might be useful as a way of saying it's really hard to figure out what these clusters mean)
# cluster_markers = top_markers(a549_cds)
# cluster_markers_to_plot = cluster_markers %>% filter(marker_test_q_value < 0.05) %>% group_by(cell_group) %>% arrange(specificity) %>% slice_head(n=5) %>% pull(gene_id) %>% unique
# plot_genes_by_group(a549_cds, markers=cluster_markers_to_plot, max.size=5)

#
# gene_set_membership_df = rowData(pseudobulk_cds) %>% as.data.frame
# gene_set_membership_df = left_join(gene_set_membership_df,
#     gene_set_df, by=c("gene_short_name" = "gene_symbol")) %>%
#     filter(gs_name %in% sig_pathway_ids) %>% dplyr::select(id, gs_name)
#
# pathway_expr_mat = monocle3::aggregate_gene_expression(a549_ccs@cds,
#                                                    #cell_group_df = tibble::rownames_to_column(a549_ccs@metadata[["cell_group_assignments"]]),
#                                                    gene_group_df = gene_set_membership_df,
#                                                    cell_group_df = cell_group_df %>% dplyr::select(cell, cell_group),
#                                                    #norm_method="size_only",
#                                                    scale_agg_values = TRUE,
#                                                    #pseudocount=0,
#                                                    gene_agg_fun="sum",
#                                                    cell_agg_fun="mean")
# pheatmap::pheatmap(pathway_expr_mat, scale="row")



plot_drug_transitions <- function(ccs,
                                  state_transition_summary,
                                  x=1,
                                  y=2,
                                  edge_size=2,
                                  cell_size=1,
                                  group_label_size=2,
                                  plot_labels = c("all", "none"),
                                  cell_groups = NULL,
                                  pcor_thresh = 0,
                                  transition_q_value_thesh = 0.01,
                                  plot_cells=FALSE,
                                  color_cells_by=NULL
                                  ){

  umap_centers = centroids(ccs)

  #cond_b_vs_a_tbl$delta_q_value = p.adjust(cond_b_vs_a_tbl$delta_q_value, method = "BH")
  # umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

  #cov_graph = hooke:::return_igraph(model(ccm))
  #corr_edge_coords = igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight != 0.00)
  corr_edge_coords = state_transition_summary %>% dplyr::select(from, to, treatment_class, num_treatments)
  corr_edge_coords = dplyr::left_join(corr_edge_coords, umap_centers %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group"))
  corr_edge_coords = dplyr::left_join(corr_edge_coords, umap_centers %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))

  corr_edge_coords = corr_edge_coords %>%
    dplyr::mutate(scaled_weight  = abs(num_treatments) / max(abs(num_treatments)))

  plot_df = ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)

  plot_df$umap2D_1 <- reducedDim(ccs@cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(ccs@cds, type="UMAP")[plot_df$cell,y]

  if (is.null(color_cells_by) == FALSE & color_cells_by %in% colnames(colData(ccs@cds))){
    plot_df$cell_color <- colData(ccs@cds)[plot_df$cell, color_cells_by]
  }else{
    plot_df$cell_color <- "white"
  }


  gp = ggplot() + monocle3:::monocle_theme_opts()
  if (plot_cells){
    gp = gp +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2),
      color = "black",
      size = 1.5 * cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2, color = cell_color),
      size = cell_size,
      stroke = 0
    )
  }



    gp = gp  +
      geom_segment(data = corr_edge_coords,
                   aes(x = from_umap_1,
                       y = from_umap_2,
                       xend=to_umap_1,
                       yend = to_umap_2,
                       size=edge_size * scaled_weight),
                   color="black"
                   ) +
      geom_segment(data = corr_edge_coords,
                   aes(x = from_umap_1,
                       y = from_umap_2,
                       xend=(from_umap_1+to_umap_1)/2,
                       yend = (from_umap_2+to_umap_2)/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_size_identity()

  if (plot_labels != "none") {
    label_df = umap_centers
    gp <- gp + ggrepel::geom_label_repel(data = label_df,
                                         mapping = aes(umap_1, umap_2, label=cell_group),
                                         size=I(group_label_size),
                                         fill = "white")
  }
  return(gp)
}
#debug(plot_drug_transitions)

drug_transition_summary = sig_state_transition_effects %>%
  filter(from_delta_q_value < 0.1 & to_delta_q_value < 0.1) %>%
  dplyr::select(from, to, treatment_class = pathway_level_2, pcor) %>%
  dplyr::group_by(from, to, treatment_class) %>%
  summarize(num_treatments = n())


rainbow_cell_state_colors =
  c("M-G1 transition" = "#81C3D7",
    "G1" = "#1965B0",
    "S" = "#F6C141",
    "G2M" = "#4EB265",
    "G0" = "#EF476F",
    "p53-dependent arrest" = "#DF4828",
    "p53-independent arrest" = "#BCE784",
    "p53 apoptosis" = "#DFE2CF",
    "Apoptosis" = "#A288E3",
    "Acetate starvation" = "#E78C35",
    "Glucocorticoid response"= "#3C787E",
    "Unknown" = "grey")



plot_drug_transitions(a549_ccs, drug_transition_summary, edge_size=0.25, plot_cells=TRUE, plot_labels="none", color_cells_by="cell_type") +
  facet_wrap(~treatment_class) +
  scale_color_manual(values=rainbow_cell_state_colors) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ggsave("drug_transition_map.png", width=24, height=24)

cluster_state_map = colData(a549_ccs@cds) %>%
  as.data.frame %>%
  dplyr::count(cluster, cell_type) %>%
  group_by(cluster) %>% slice_max(n)%>% dplyr::select(cluster, cell_type)

# drug_sensitivity_map = sig_state_transition_effects %>%
#   filter(from_delta_q_value < 0.1 & to_delta_q_value < 0.1) %>%
#   dplyr::select(from, to, product_name, treatment_class = pathway_level_2, target, pcor)
#
# drug_sensitivity_map = left_join(drug_sensitivity_map, cluster_state_map, by=c("from"="cluster")) %>% dplyr::rename(from_state=cell_type)
# drug_sensitivity_map = left_join(drug_sensitivity_map, cluster_state_map, by=c("to"="cluster")) %>% dplyr::rename(to_state=cell_type)
# #drug_sensitivity_map = drug_sensitivity_map %>% dplyr::select(from, to, treatment_class, from_state, to_state)
#
# drug_sensitivity_map = drug_sensitivity_map %>% dplyr::select(product_name, treatment_class, target, from_state, to_state) %>% distinct()
#
#
# multi_destination_drugs = drug_sensitivity_map %>% group_by(product_name) %>%
#   summarize(num_destination_states = length(unique(to_state))) %>% filter(num_destination_states > 1) %>% pull(product_name)
# multi_destination_drugs = drug_sensitivity_map %>% filter(product_name %in% multi_destination_drugs)
# ggplot(aes(x=treatment_class), data=multi_destination_drugs) + geom_bar() + facet_grid(to_state ~ from_state) + coord_flip()

target_pathways = unique(sig_state_transition_effects$pathway_level_2)
for (i in (1:length(target_pathways))){
  target_pathway = target_pathways[i]
  product_transitions = sig_state_transition_effects %>%
    filter(from_delta_q_value < 0.1 & to_delta_q_value < 0.1 & pathway_level_2 == target_pathway) %>%
    dplyr::select(from, to, treatment_class = product_name, pcor) %>%
    dplyr::group_by(from, to, treatment_class) %>%
    summarize(num_treatments = n())

  fig_size = ceiling(sqrt(3*length(unique(product_transitions$treatment_class))))
  print (fig_size)
  fig_size = max(6, fig_size)
  filename = paste(stringr::str_replace_all(target_pathway,"/", "-"),"_transition_map.png", sep="")
  plot_drug_transitions(a549_ccs, product_transitions, edge_size=0.25, plot_cells=TRUE, plot_labels="none", color_cells_by="cell_type") +
    facet_wrap(~treatment_class) +
    scale_color_manual(values=rainbow_cell_state_colors) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    ggsave(filename, width=fig_size+1, height=fig_size)

}



##### Extra junk

sig_state_transition_effects %>%
  filter(from_delta_q_value < 0.1 & to_delta_q_value < 0.1) %>%
  dplyr::select(from, to, product_name, target, treatment_class = pathway_level_2, pcor) %>%
  dplyr::group_by(from, to, treatment_class) %>% filter(to == 26)

sig_state_transition_effects %>%
  filter(from_delta_q_value < 0.1 & to_delta_q_value < 0.1) %>%
  dplyr::select(from, to, product_name, target, treatment_class = pathway_level_2, pcor) %>%
  dplyr::group_by(from, to, treatment_class) %>% filter(treatment_class == "Bromodomain")


