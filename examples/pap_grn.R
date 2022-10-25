library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

# First, let's just look at the how cells shift between states in PAP vs healthy controls:

pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/pap.al_2021-09-15.cds.RDS")

pap_cds = detect_genes(pap_cds)

# assign best celltype column and reduce dims
colData(pap_cds)$cell_type = colData(pap_cds)$CW_assignedCellType
colData(pap_cds)$cluster = monocle3::clusters(pap_cds)

colData(pap_cds)$Size_Factor = size_factors(pap_cds)


colData(pap_cds)$experiment = colData(pap_cds)$sample
colData(pap_cds)$sample = NULL

plot_cells(pap_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE)

plot_cells(pap_cds, color_cells_by="Genotype", show_trajectory_graph=FALSE) + facet_wrap(~Genotype)


plot_cells(pap_cds, genes=c("Csf2",
                            "Csf2ra",
                            "Csf2rb",
                            "Chil3",
                            "Lpl",
                            "Car4",
                            "Apoe",
                            "Fabp4",
                            "Fabp5",
                            "F13a1"),
           show_trajectory_graph=FALSE)

#plot_cells(cds)

ccs = new_cell_count_set(pap_cds,
                         sample_group = "sampleName",
                         cell_group = "cell_type")


# Fit a Hooke model:
ccm  = new_cell_count_model(ccs,
                            model_formula_str = "~Genotype + batch")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=0.1)

# Now we can compare the various mutant mice to WT to see how cells shift around.

cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_csf2rb = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))


cond_ra_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)

# This plot shows a big shift from healthy to PAP macrophages
plot_contrast(ccm, cond_ra_vs_wt_tbl, scale_shifts_by="none", p_value_thresh=0.05)

cond_rb_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2rb)

plot_contrast(ccm, cond_rb_vs_wt_tbl, scale_shifts_by="none", p_value_thresh=0.05)

pred_abund_mat = cbind(cond_wt$log_abund, cond_csf2ra$log_abund, cond_csf2rb$log_abund)
colnames(pred_abund_mat) = c("WT", "Csfr2ra-/-", "Csfr2rb-/-")

pheatmap::pheatmap(pred_abund_mat, scale="row")

######

############# Workflow for identifying regulators that mediate shifts in cell state ##########

# Let's use PLN networks to find putative regulators that might mediate the shift from healthy to
# PAP macrophage state

# Generate a pseudobulk expression matrix
agg_expr_mat = monocle3::aggregate_gene_expression(ccs@cds,
                                                   cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]]),
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0,
                                                   cell_agg_fun="mean")


agg_coldata = ccs@metadata[["cell_group_assignments"]] %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
row.names(agg_coldata) = colnames(agg_expr_mat)

### Force expression of KO genes to zero (as these are not functional transcripts)
experiment_samples = ccs@metadata[["cell_group_assignments"]] %>% pull(sample) %>% unique()

ablation_df = rbind(
  data.frame(gene_id = "ENSMUSG00000059326", #Csf2ra, alpha chain of GM-CSF receptor
             groups_to_ablate = ccs@metadata[["cell_group_assignments"]] %>%
               dplyr::filter(sample %in% experiment_samples[grepl("Csf2ra", experiment_samples)]) %>% pull(group_id)),
  data.frame(gene_id = "ENSMUSG00000071713", #Csf2ra, beta chain of GM-CSF receptor
             groups_to_ablate = ccs@metadata[["cell_group_assignments"]] %>%
               dplyr::filter(sample %in% experiment_samples[grepl("Csf2rb|RBCKO", experiment_samples)]) %>% pull(group_id)),
  data.frame(gene_id = "ENSMUSG00000018916", #Csf2, GM-CSF itself
             groups_to_ablate = ccs@metadata[["cell_group_assignments"]] %>%
               dplyr::filter(sample %in% experiment_samples[grepl("GMKO", experiment_samples)]) %>% pull(group_id))
)

### This function takes as input and zeros out the specified genes in the specified samples
ablate_expression <- function(pseudobulk_expr_mat, gene_by_condition_df){
  pseudobulk_expr_mat[which(row.names(agg_expr_mat) %in% ablation_df$gene_id),
                      which(colnames(agg_expr_mat) %in% ablation_df$groups_to_ablate)] = 0
  return(pseudobulk_expr_mat)
}

agg_expr_mat = ablate_expression(agg_expr_mat, experiment_samples)


pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
pseudobulk_cds = preprocess_cds(pseudobulk_cds)
pseudobulk_cds = reduce_dimension(pseudobulk_cds)

###################

# Perform contrast between KO and WT

cond_ra_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)

# Identify reciprocal shifts between cell states:

# TODO: We should explore different policies for the step below. Choices to consider:
# - significant fold changes at one endpoint or both (currently both)?
# - Only look at endpoint states for directed edges? Or also at states along undirected paths between directed edges? Latter will might give us way more genes (possibly too many)
state_transitions = hooke:::collect_pln_graph_edges(ccm, cond_ra_vs_wt_tbl) %>% as_tibble %>%
  filter(edge_type == "directed_from_to" & to_delta_p_value < 0.05 & from_delta_p_value < 0.05)


collect_transition_states = function(from_state, to_state, ccm){
  # This could be expanded to look collect all states along a path from from_state to to_state in the CCM's PLN network
  return (as.character(unique(c(from_state, to_state))))
}

St_f = state_transitions %>%
  dplyr::mutate(states = purrr::map2(.f = purrr::possibly(
    collect_transition_states, NA_real_), .x = from,
    .y = to, ccm)) #%>%
  #tidyr::unnest(states)

# This function compares two cell states to find genes that differ between them
find_degs_between_states = function(states,
                                    cds,
                                    model_formula_str = "~cell_group",
                                    gene_whitelist = NULL,
                                    q_value_thresh=0.05,
                                    effect_thresh = 2,
                                    cores=1){
  cds_states_tested = cds[,as.character(colData(cds)$cell_group) %in% states]
  models = monocle3::fit_models(cds_states_tested,
                                model_formula_str=model_formula_str,
                                weights=colData(cds_states_tested)$num_cells_in_group,
                                cores=cores)
  coefs = monocle3::coefficient_table(models) %>%
    filter(grepl("cell_group", term) & q_value < q_value_thresh & abs(normalized_effect) > effect_thresh)
  gene_ids = coefs %>% pull(gene_id)
  gene_ids = c(gene_ids, gene_whitelist) %>% unique
  return(gene_ids)
}
#debug(find_degs_between_states)

# For each reciprocal state transition, find DEGs between the endpoints
St_f = St_f %>%
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(
    find_degs_between_states, NA_real_), .x = states, pseudobulk_cds, gene_whitelist=c("ENSMUSG00000059326","ENSMUSG00000071713","ENSMUSG00000018916"), cores=4)) #%>%
  #tidyr::unnest(degs)

# Build a PLN network on genes from a bunch of pseudobulks
# FIXME: should this use weights?
build_pln_model_on_genes = function(genes,
                                    states,
                                    cds,
                                    regulatory_genes = NULL,
                                    model_formula_str = "~cell_state",
                                    gene_module_df = NULL,
                                    sparsity_factor=1,
                                    pln_min_ratio=0.001,
                                    pln_num_penalties=30){

  #test_module_expr_cds = pseudobulk_cds[all_gene_modules %>% filter(module %in% test_modules) %>% pull(id),
  #                                      colData(pseudobulk_cds)$cell_group %in% test_states]
  pb_cds = cds[genes, as.character(colData(cds)$cell_group) %in% states]

  gene_expr_cds = new_cell_data_set(t(t(counts(pb_cds)) * colData(pb_cds)$num_cells_in_group),
                                    cell_metadata = colData(pb_cds) %>% as.data.frame,
                                    gene_metadata = rowData(pb_cds) %>% as.data.frame)

  colData(gene_expr_cds)$sample = colData(gene_expr_cds)$group_id
  colData(gene_expr_cds)$cell_state  = as.character(colData(gene_expr_cds)$cell_group)
  colData(gene_expr_cds)$cell_group = NULL

  gene_ccs = methods::new("cell_count_set",
                                 gene_expr_cds,
                                 cds=pap_cds)
  gene_ccs = gene_ccs[,colSums(counts(gene_ccs)) > 0]

  init_gene_penalty_matrix = function(gene_module_df, base_penalty = 1, min_penalty=0.01, max_penalty=1e6){

    coord_matrix = gene_module_df %>% dplyr::select(id, dim_1, dim_2)
    centroid_coords = aggregate(.~id, data=coord_matrix, FUN=mean)
    colnames(centroid_coords)[-1] = paste0(tolower("UMAP"), "_", 1:(length(colnames(centroid_coords))-1))
    dist_matrix = as.matrix(dist(centroid_coords[,-1], method = "euclidean", upper=T, diag = T))

    row.names(dist_matrix) <- centroid_coords$id
    colnames(dist_matrix) <- centroid_coords$id

    # TODO: do I need this? Probably the caller can and should do this.
    #dist_matrix = dist_matrix[colnames(data$Abundance), colnames(data$Abundance)]

    get_rho_mat <- function(DM, distance_parameter = 1, s=1, xmin = NULL) {
      if (is.null(xmin)){
        xmin = min(DM[DM > 0]) / 2
      }
      #out <- (1-(xmin/DM)^s) * distance_parameter
      out = min_penalty + (DM / max(DM))^2
      #out =  1 + DM^s
      # penalties have to be > 0
      out[!is.finite(out)] <- min_penalty
      out[out < 0] <- min_penalty
      return(out)
    }

    penalty_matrix = base_penalty * (min_penalty + get_rho_mat(dist_matrix, distance_parameter=1, s=2))

    # TODO: add support for whitelisting and blacklisting
    #qplot(as.numeric(dist_matrix), as.numeric(out))
    return(penalty_matrix)
  }

  # Ok let's make a model of how genes co-vary across pseudobulks
  if (is.null(gene_module_df) == FALSE){
    module_penalty_matrix = matrix(1, nrow(gene_expr_cds), nrow(gene_expr_cds))
    row.names(module_penalty_matrix) = row.names(gene_expr_cds)
    colnames(module_penalty_matrix) = row.names(gene_expr_cds)
    co_express_penalty_matrix = init_gene_penalty_matrix(gene_module_df %>% dplyr::filter(id %in% row.names(gene_expr_cds)))
    module_penalty_matrix[row.names(co_express_penalty_matrix), colnames(co_express_penalty_matrix)] = co_express_penalty_matrix
    #co_express_penalty_matrix = co_express_penalty_matrix[row.names(gene_expr_cds), row.names(gene_expr_cds)]
  }else{
    module_penalty_matrix = matrix(1, nrow(gene_expr_cds), nrow(gene_expr_cds))
    row.names(module_penalty_matrix) = row.names(gene_expr_cds)
    colnames(module_penalty_matrix) = row.names(gene_expr_cds)
  }

  if (is.null(regulatory_genes) == FALSE){
    non_regulator_genes = rowData(gene_expr_cds) %>% as.data.frame %>% pull(id) %>% setdiff(regulatory_genes)
    blacklist = expand.grid(non_regulator_genes, non_regulator_genes)
  } else{
    blacklist = NULL
  }

  nc_ig = colData(gene_ccs)$num_cells_in_group
  model_weights =  nc_ig / exp(mean(log(nc_ig)))
  gene_ccm = new_cell_count_model(gene_ccs,
                                  model_formula_str,
                                  penalty_matrix = module_penalty_matrix,
                                  blacklist = blacklist,
                                  #weights=model_weights,#weights=colData(gene_ccs)$num_cells_in_group,
                                  verbose=TRUE,
                                  pseudocount=1e-2,
                                  pln_min_ratio=pln_min_ratio,
                                  pln_num_penalties=pln_num_penalties)
  gene_ccm = select_model(gene_ccm, sparsity_factor=sparsity_factor)
  return(gene_ccm)
}
#debug(build_pln_model_on_genes)

# Now let's make a whitelist of regulators:
gene_set = msigdbr(species = "Mus musculus", subcategory = "GO:MF")
transcription_regulators = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("_TRANSCRIPTION", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

receptors = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("_RECEPTOR", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

kinases = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("_KINASE", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

# We will allow TFs and kinases to have edges to other genes. No
# Edges between genes that aren't one of these two classes
allowed_regulator_symbols = c(transcription_regulators,
                              #receptors,
                              kinases)
allowed_regulator_ids = rowData(pseudobulk_cds) %>% as.data.frame %>% filter(gene_short_name %in% allowed_regulator_symbols) %>% pull(id)
allowed_regulator_ids = c(allowed_regulator_ids, c("ENSMUSG00000059326","ENSMUSG00000071713","ENSMUSG00000018916"))

# Project genes into 2D UMAP space. We'll use this to set up penalties between genes
all_degs = c(unlist(St_f$degs), allowed_regulator_ids) %>% unique
all_gene_modules = find_gene_modules(pseudobulk_cds[all_degs,], resolution=1e-2)
all_gene_modules = left_join(all_gene_modules, rowData(pseudobulk_cds) %>% as.data.frame, by="id")
qplot(dim_1, dim_2, color=module, data=all_gene_modules)

# OK now fit the PLN for each reciprocal transition:
St_f = St_f %>%
  dplyr::mutate(gene_count_model = purrr::map2(.f = purrr::possibly(
    build_pln_model_on_genes, NA_real_),
    .x = degs,
    .y = states,
    pseudobulk_cds,
    model_formula_str="~cell_state",
    regulatory_genes = allowed_regulator_ids,
    gene_module_df = all_gene_modules,
    sparsity_factor=1,
    pln_min_ratio=1e-3)) #%>%
#tidyr::unnest(states)

# Rank the genes based on their degree in the PLN network, with edges weighted
# by partial correlation
rank_regulators = function(states, gene_model_ccm, cds, log_abundance_thresh=1e-5){
  from_state = states[1]
  to_state = states[2]

  from_cond_est = estimate_abundances(gene_model_ccm, tibble::tibble(cell_state=from_state))
  to_cond_est = estimate_abundances(gene_model_ccm, tibble::tibble(cell_state=to_state))

  from_to_cond_diff = compare_abundances(gene_model_ccm, from_cond_est, to_cond_est)

  # "Activators" are genes that are upregulated in the transition and positively correlated with
  # many other upregulated genes, or negatively correlated with downregulated genes
  # "Repressors" are the opposite. Genes that are downregulated and positively correlated with
  # other downregulated genes, or negatively correlated with upregulated genes

  # get ids of upregulated genes
  activators = from_to_cond_diff %>% dplyr::filter(delta_log_abund > 0) %>% pull(cell_group)

  # get the subgraph induced by upregulated genes
  upreg_igraph = hooke:::return_igraph(model(gene_model_ccm))
  upreg_igraph = igraph::induced_subgraph(upreg_igraph, igraph::V(upreg_igraph)[activators])
  upreg_igraph = igraph::delete_edges(upreg_igraph, igraph::E(upreg_igraph)[weight < 0])

  activator_scores = igraph::strength(upreg_igraph)
  activator_score_df = tibble(gene_id = names(activator_scores), regulator_score = activator_scores)

  # get ids of downregulated genes
  represssors = from_to_cond_diff %>% dplyr::filter(delta_log_abund < 0) %>% pull(cell_group)

  # get the subgraph induced by upregulated genes
  downreg_igraph = hooke:::return_igraph(model(gene_model_ccm))
  downreg_igraph = igraph::induced_subgraph(downreg_igraph, igraph::V(downreg_igraph)[represssors])
  downreg_igraph = igraph::delete_edges(downreg_igraph, igraph::E(downreg_igraph)[weight < 0])

  repressor_scores = -igraph::strength(downreg_igraph)
  repressor_score_df = tibble(gene_id = names(repressor_scores), regulator_score = repressor_scores)

  score_df = from_to_cond_diff %>% dplyr::select("gene_id"=cell_group, delta_log_abund)
  score_df = dplyr::left_join(score_df, rbind(activator_score_df, repressor_score_df))
  return(score_df)
}
#debug(rank_regulators)

St_f = St_f %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(
    rank_regulators, NA_real_), .x = states,
    .y = gene_count_model, pseudobulk_cds)) #%>%
#tidyr::unnest(states)

# Now plot a comparison of the fold change of each gene vs. its "regulator score"
plot_top_regulators = function(regulator_df, cds, group_label_size=2, resid_std_devs=1, abund_lfc_range=c(-5,5)){
  regulator_df$label = rowData(cds)[regulator_df$gene_id,]$gene_short_name
  regulator_df = regulator_df %>% mutate(delta_log_abund = ifelse(delta_log_abund < min(abund_lfc_range),
                                                                  min(abund_lfc_range),
                                                                  delta_log_abund))
  regulator_df = regulator_df %>% mutate(delta_log_abund = ifelse(delta_log_abund > max(abund_lfc_range),
                                                                  max(abund_lfc_range),
                                                                  delta_log_abund))
  reg_score_fit = glm(regulator_score ~ splines::ns(delta_log_abund, df=3), data=regulator_df)
  top_regulators = which(abs(rstandard(reg_score_fit)) > resid_std_devs)

  gp = ggplot(aes(x=delta_log_abund, y=regulator_score), data=regulator_df) + geom_point()
  gp = gp + monocle3:::monocle_theme_opts()
  gp = gp + ggrepel::geom_label_repel(data = regulator_df[top_regulators,],
                                      mapping = aes(x=delta_log_abund, y=regulator_score, label=label),
                                      size=I(group_label_size),
                                      fill = "white")
  return(gp)
}
#debug(plot_top_regulators)
plot_top_regulators(St_f$regulators[[1]], pseudobulk_cds)

get_neighboring_genes = function(gene_model_ccm, pseudobulk_cds, gene_short_name, pos_p_cor_only=TRUE, min_abs_pcor=0.001){
  upreg_igraph = hooke:::return_igraph(model(gene_model_ccm))
  upreg_igraph = igraph::induced_subgraph(upreg_igraph, igraph::V(upreg_igraph))
  upreg_igraph = igraph::delete_edges(upreg_igraph, igraph::E(upreg_igraph)[abs(weight) < min_abs_pcor])
  if (pos_p_cor_only)
    upreg_igraph = igraph::delete_edges(upreg_igraph, igraph::E(upreg_igraph)[weight < 0])
  source_gene_id = rowData(pseudobulk_cds)[rowData(pseudobulk_cds)$gene_short_name == gene_short_name,] %>% as.data.frame %>% pull(id)
  n_ids = igraph::V(upreg_igraph)[igraph::neighbors(upreg_igraph, source_gene_id)]$name
  rowData(pseudobulk_cds)[n_ids,] %>% as.data.frame %>% pull(gene_short_name)
}
#debug(get_neighboring_genes)
get_neighboring_genes(St_f$gene_count_model[[1]], pseudobulk_cds, "Csf2ra", min_abs_pcor=0.01)

get_neighboring_genes(St_f$gene_count_model[[1]], pseudobulk_cds, "Cebpb", min_abs_pcor=0.01)
get_neighboring_genes(St_f$gene_count_model[[1]], pseudobulk_cds, "Zeb2", min_abs_pcor=0.01)

# TODO: after we've done the above for a bunch of transitions, we should be able to pull out the top regulators for each
# and assemble them into a single regulatory network that explains all the shifts!


######## Some diagnostic stuff: ##########

test_gene_model_ccm = select_model(St_f$gene_count_model[[1]], sparsity_factor=1)

#test_gene_model_ccm = select_model(St_f$gene_count_model[[1]], sparsity_factor=100)


test_gene_cond_hpap_am = estimate_abundances(test_gene_model_ccm, tibble::tibble(cell_state="hPAP alveolar macrophages"))
test_gene_cond_healthy_am = estimate_abundances(test_gene_model_ccm, tibble::tibble(cell_state="Healthy alveolar macrophages"))

test_cond_healthy_vs_hpap_am_tbl = compare_abundances(test_gene_model_ccm, test_gene_cond_healthy_am, test_gene_cond_hpap_am)

plot(model(test_gene_model_ccm),  output = "corrplot")


coefficient_path(test_gene_model_ccm@model_family, corr = TRUE) %>%
  ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) +
  geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw()

plot_gene_contrast <- function(ccm,
                               #umap_centers,
                               #gene_module_df,
                               cond_b_vs_a_tbl,
                               log_abundance_thresh = -5,
                               scale_shifts_by=c("receiver", "sender", "none"),
                               #cell_group="cluster",
                               edge_size=2,
                               cell_size=1,
                               p_value_thresh = 1.0,
                               group_label_size=2,
                               max_effect_size = 4,
                               min_effect_size = -4,
                               mask_pcor_range = c(0, 0.05),
                               genes_to_label=NULL){

  #umap_centers = centroids(ccm@ccs)

  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund = ifelse(delta_p_value <= p_value_thresh, delta_log_abund, 0))
  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% mutate(delta_log_abund = ifelse(delta_log_abund > max_effect_size, max_effect_size, delta_log_abund),
                                               delta_log_abund = ifelse(delta_log_abund < min_effect_size, min_effect_size, delta_log_abund))

  ccm_pcor_graph = hooke:::return_igraph(model(ccm))
  ccm_pcor_graph = igraph::delete_edges(ccm_pcor_graph, igraph::E(ccm_pcor_graph)[igraph::E(ccm_pcor_graph)$weight > mask_pcor_range[2] | igraph::E(ccm_pcor_graph)$weight < mask_pcor_range[1]])
  ccm_pcor_graph_laplacian_embed <- igraph::embed_laplacian_matrix(ccm_pcor_graph, 5)
  ccm_pcor_network_umap = uwot::umap(ccm_pcor_graph_laplacian_embed$X)
  #qplot(xxx_network_umap[,1],
  #      xxx_network_umap[,2],
  #      label=rowData(pseudobulk_cds)[igraph::V(xxx)$name,] %>% as.data.frame %>% pull(gene_short_name), geom="text")

  #gene_module_df$cell_group = gene_module_df$module
  ccm_pcor_network_umap = as.data.frame(ccm_pcor_network_umap)
  ccm_pcor_network_umap$cell_group = igraph::V(ccm_pcor_graph)$name#rowData(ccm@ccs)[igraph::V(ccm_pcor_graph)$name,] %>% as.data.frame %>% pull(gene_short_name)
  ccm_pcor_network_umap = ccm_pcor_network_umap[,c(3,1,2)]
  #coord_matrix = gene_module_df %>% dplyr::select("cell_group"=id, dim_1, dim_2)
  #centroid_coords = aggregate(.~cell_group, data=coord_matrix, FUN=mean)
  colnames(ccm_pcor_network_umap)[-1] = paste0(tolower("UMAP"), "_", 1:(length(colnames(ccm_pcor_network_umap))-1))
  umap_centers = ccm_pcor_network_umap

  umap_centers_delta_abund = umap_centers
  umap_centers_delta_abund = dplyr::left_join(umap_centers_delta_abund, cond_b_vs_a_tbl, by=c("cell_group"="cell_group"))
  umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  umap_centers_delta_abund = umap_centers_delta_abund %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))
  umap_centers_delta_abund = left_join(umap_centers_delta_abund, rowData(ccm@ccs) %>% as.data.frame %>% dplyr::select(id, gene_short_name), by=c("cell_group"="id"))
  #cond_b_vs_a_tbl$delta_q_value = p.adjust(cond_b_vs_a_tbl$delta_p_value, method = "BH")

  corr_edge_coords_umap_delta_abund = hooke:::collect_pln_graph_edges(ccm,
                                                                      umap_centers_delta_abund,
                                                                      log_abundance_thresh)

  corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund %>% filter(pcor < min(mask_pcor_range) | pcor > max(mask_pcor_range))
  positive_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(pcor > 0)
  negative_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(pcor < 0)

  positive_edge_df = positive_edge_df %>%
    dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  negative_edge_df = negative_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))


  plot_df = ccm_pcor_network_umap

  #cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% mutate(cluster = stringr::str_split_fixed(cell_group, "\\.", 3)[,3])
  plot_df = dplyr::left_join(plot_df,
                             cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                             by=c("cell_group"="cell_group"))

  #directed_edge_df = directed_edge_df %>% filter(edge_type %in% c("directed_to_from", "directed_from_to"))


  gp = ggplot()  +
    geom_segment(data = positive_edge_df,
                 aes(x = to_umap_1,
                     y = to_umap_2,
                     xend=from_umap_1,
                     yend = from_umap_2,
                     alpha=scaled_weight),
                 color="tomato") +
    geom_segment(data = negative_edge_df,
                 aes(x = to_umap_1,
                     y = to_umap_2,
                     xend=from_umap_1,
                     yend = from_umap_2,
                     alpha=scaled_weight),
                 color="steelblue")

  # plot
  gp = gp +
    geom_point(
      data = plot_df,
      aes(umap_1, umap_2),
      color = "black",
      size = 1.5 * cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df,
      aes(umap_1, umap_2),
      color = "white",
      size = cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>%
        arrange(!is.na(abs(delta_log_abund)),
                abs(delta_log_abund)),
      aes(umap_1, umap_2, color = delta_log_abund),
      size = cell_size,
      stroke = 0
    ) +
    scale_color_gradient2(
      low = "#122985",
      mid = "white",
      high = "red4",
      na.value = "white"
    )  +
    theme_void() +
    #theme(legend.position = "none") +

    #plot_oriented_pln_network(directed_edge_df) +

    #geom_point(data = umap_centers_delta_abund, aes(umap2D_1, umap2D_2, color=delta_log_abund, size=max_log_abund)) +
    # #geom_text(data = umap_centers_delta_abund, aes(umap2D_1, umap2D_2, color=delta_log_abund, label=label)) +
    #scale_color_gradient2(low = 'steelblue', mid = 'grey', high = 'darkred') +
    # #theme(legend.position = "none") +
    monocle3:::monocle_theme_opts() #+ ggtitle(paste(cond_b, "vs",cond_a, "pos partial corr w/ umap penalty"))


  #if(label_cell_groups) {
  if (is.null(genes_to_label)){
    label_df = umap_centers_delta_abund %>% filter(delta_log_abund != 0)
  }
  else{
    label_df = umap_centers_delta_abund %>% filter(delta_log_abund != 0) %>% filter(gene_short_name %in% genes_to_label)
  }
  gp <- gp + ggrepel::geom_label_repel(data = label_df,
                                       mapping = aes(umap_1, umap_2, label=gene_short_name),
                                       size=I(group_label_size),
                                       fill = "white")

  return(gp)
}
debug(plot_gene_contrast)
# gene_callouts = rowData(test_gene_model_ccm@ccs) %>% as.data.frame %>% filter(TF == TRUE) %>% pull(gene_short_name)
# gene_callouts = c(gene_callouts, "Csf2ra", "Cfs2", "Csf2rb")
# plot_gene_contrast(test_gene_model_ccm,
#                           #gene_subset_umap_df %>% filter(id %in% row.names(rowData(test_gene_model_ccm@ccs))),
#                           test_cond_healthy_vs_hpap_am_tbl,
#                           scale_shifts_by="none",
#                           p_value_thresh=1,
#                    edge_size = 0.25,
#                           cell_size=5,
#                    mask_pcor_range = c(0,0.1),
#                    genes_to_label = gene_callouts)






