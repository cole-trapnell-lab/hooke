
# given gene short names, return a dataframe that contains 
# the samples in which to ablate the counts
#' @param ccs a cell_count_set object
#' @param genes a list of genes to ablate
get_ablation_samples <- function(ccs, genes) {
  
  experiment_samples = ccs@metadata[["cell_group_assignments"]] %>% pull(sample) %>% unique()
  
  get_samples = function(gene_short_name, ccs) {
    
    experiment_samples = ccs@metadata[["cell_group_assignments"]] %>% pull(sample) %>% unique()
    groups_to_ablate = ccs@metadata[["cell_group_assignments"]] %>%
      dplyr::filter(sample %in% experiment_samples[grepl(gene_short_name, experiment_samples)]) %>% 
      pull(group_id)
    return(groups_to_ablate)
  }
  
  ablation_df = rowData(ccs@cds) %>% 
    as.data.frame %>% 
    filter(gene_short_name %in% genes) %>% 
    rownames_to_column("gene_id") %>% 
    select(gene_id, gene_short_name) %>%
    mutate(groups_to_ablate = purrr::map(.f = get_samples, 
                                         .x = gene_short_name,ccs)) %>% 
    tidyr::unnest(groups_to_ablate)
  
  return(ablation_df)
}


### This function takes as input and zeros out the specified genes in the specified samples
#' @param agg_expr_mat aggregate gene expression matrix
#' @param ablation_df data frame of specified genes + samples to be zeroed out
ablate_expression <- function(agg_expr_mat, ablation_df){
  agg_expr_mat[which(row.names(agg_expr_mat) %in% ablation_df$gene_id),
               which(colnames(agg_expr_mat) %in% ablation_df$groups_to_ablate)] = 0
  return(agg_expr_mat)
}

#' Returns a pseudobulked cell count dataset object based on cell count set groupings
#' @param ccs a cell count set object
#' @param gene_ids gene ids to ablate
#' @param norm_method
#' @param scale_agg_values
#' @param cell_agg_fun
#' @param cell_group_threshold 
#' @param cell_state_threshold filters cell groups that have low counts
#' @export
pseudobulk <- function(ccs, 
                       gene_ids = NULL,
                       norm_method="size_only",
                       scale_agg_values = FALSE,
                       pseudocount = 0, 
                       cell_agg_fun = "mean",
                       total_cell_threshold = 25,
                       num_cells_in_group_threshold = 5) {
  
  # get rid of samples with low counts
  ccs = ccs[, Matrix::colSums(counts(ccs)) >= total_cell_threshold]
  
  agg_expr_mat = monocle3::aggregate_gene_expression(ccs@cds,
                                                     cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]]),
                                                     norm_method = norm_method,
                                                     scale_agg_values = scale_agg_values,
                                                     pseudocount = pseudocount,
                                                     cell_agg_fun = cell_agg_fun)
  
  
  agg_coldata = ccs@metadata[["cell_group_assignments"]] %>%
    dplyr::group_by(group_id, cell_group) %>%
    dplyr::summarize(num_cells_in_group = n()) %>%
    as.data.frame %>% 
    filter(num_cells_in_group >= num_cells_in_group_threshold)
  
  agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
  row.names(agg_coldata) = colnames(agg_expr_mat)

  if (!is.null(gene_ids)) {
    ablation_df = get_ablation_samples(ccs, gene_ids)
    agg_expr_mat = ablate_expression(agg_expr_mat, ablation_df)
  }
  
  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
  pseudobulk_cds = pseudobulk_cds[,Matrix::colSums(exprs(pseudobulk_cds)) != 0]
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  pseudobulk_cds = preprocess_cds(pseudobulk_cds)
  pseudobulk_cds = reduce_dimension(pseudobulk_cds)

}

#' @param from_state
#' @param to_state
collect_transition_states = function(from_state, to_state){
  return (as.character(unique(c(from_state, to_state))))
}


collect_between_transition_states = function(shortest_path) {
  states = shortest_path %>% select(from,to) %>% t %>% c 
  return(as.character(unique(states)))
}


# This function compares two cell states to find genes that differ between them
#' @param states
#' @param cds
#' @param model_formula_str
#' @param gene_whitelist
#' @param q_value_thresh
#' @param effect_thresh
#' @param cores
find_degs_between_states = function(states,
                                    cds,
                                    model_formula_str = "~cell_group",
                                    gene_whitelist = NULL,
                                    q_value_thresh = 0.05,
                                    effect_thresh = 2,
                                    cores = 1) {
  cds_states_tested = cds[,as.character(colData(cds)$cell_group) %in% states]
  models = monocle3::fit_models(cds_states_tested,
                                model_formula_str=model_formula_str,
                                weights=colData(cds_states_tested)$num_cells_in_group,
                                cores=cores)
  coefs = monocle3::coefficient_table(models) %>%
    filter(!grepl("Intercept", term) & 
           q_value < q_value_thresh & 
           abs(normalized_effect) > effect_thresh)
  gene_ids = coefs %>% pull(gene_id)
  gene_ids = c(gene_ids, gene_whitelist) %>% unique
  return(gene_ids)
}


#' @param genes
#' @param states only build model on specified cell states 
#' @param pb_cds A pseudobulked Monocole cell data set object
#' @param regulatory_genes
#' @param model_str
#' @param gene_module_df
#' @param sparsity_factor
#' @param pln_min_ratio
#' @param pln_num_penalties
build_pln_model_on_genes = function(genes,
                                    states,
                                    cds,
                                    regulatory_genes = NULL,
                                    model_formula_str = "~cell_state",
                                    gene_module_df = NULL,
                                    sparsity_factor=1,
                                    pln_min_ratio=0.001,
                                    pln_num_penalties=30){
  
  pb_cds = cds[genes, as.character(colData(cds)$cell_group) %in% states]
  
  gene_expr_cds = new_cell_data_set(t(t(as.matrix(counts(pb_cds))) * colData(pb_cds)$num_cells_in_group),
                                    cell_metadata = colData(pb_cds) %>% as.data.frame,
                                    gene_metadata = rowData(pb_cds) %>% as.data.frame)
  
  colData(gene_expr_cds)$sample = colData(gene_expr_cds)$group_id
  colData(gene_expr_cds)$cell_state  = as.character(colData(gene_expr_cds)$cell_group)
  colData(gene_expr_cds)$cell_group = NULL
  
  gene_ccs = methods::new("cell_count_set",
                          gene_expr_cds,
                          # cds=pap_cds, 
                          cds = pb_cds) # is this right?  
  gene_ccs = gene_ccs[,Matrix::colSums(counts(gene_ccs)) > 0]
  
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


#' Rank the genes based on their degree in the PLN network, with edges weighted
#' by partial correlation
#' @param states
#' @param gene_model_ccm a cell_count_model object
#' @param cds A Monocle cell data set object
#' @param log_abundance_thresh
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



#' Plots a comparison of the fold change of each gene vs. its "regulator score" 
#' @param regulator_df
#' @param cds
#' @param group_label_size
#' @param resid_std_devs
#' @param abund_lfc_range
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
