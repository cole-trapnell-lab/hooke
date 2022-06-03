
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

  # not sure the best place to put this yet
  colData(pseudobulk_cds)$cell_group = as.character(colData(pseudobulk_cds)$cell_group)

  return(pseudobulk_cds)
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
                                    model_formula_str = "~ cell_group",
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


#'
#' @param pb_cds
#' @param knockout_genes
#'
get_allowed_regulators = function(pb_cds, knockout_genes, allow_regulators = c("transcription")) {

  gene_set = msigdbr(species = "Danio rerio", subcategory = "GO:MF")

  transcription_regulators = list()
  receptors = list()
  kinases = list()
  if ("transcription" %in% allow_regulators){
    transcription_regulators = gene_set %>%
      dplyr::select(gs_id, gene_symbol, gs_name) %>%
      dplyr::filter(grepl("_TRANSCRIPTION", gs_name, ignore.case=TRUE)) %>%
      pull(gene_symbol) %>% unique %>% sort
  }

  if ("receptors" %in% allow_regulators) {
    receptors = gene_set %>%
      dplyr::select(gs_id, gene_symbol, gs_name) %>%
      dplyr::filter(grepl("_RECEPTOR", gs_name, ignore.case=TRUE)) %>%
      pull(gene_symbol) %>% unique %>% sort

  }
  if ("kinases" %in% allow_regulators) {
    kinases = gene_set %>%
      dplyr::select(gs_id, gene_symbol, gs_name) %>%
      dplyr::filter(grepl("_KINASE", gs_name, ignore.case=TRUE)) %>%
      pull(gene_symbol) %>% unique %>% sort
  }

  allowed_regulator_symbols = c(transcription_regulators,
                                kinases, receptors)

  allowed_regulator_ids = rowData(pb_cds) %>%
    as.data.frame %>%
    filter(gene_short_name %in% allowed_regulator_symbols) %>%
    pull(id)

  allowed_regulator_ids = c(allowed_regulator_ids, knockout_genes)
  return(allowed_regulator_ids)

}

#'
#' @param pb_cds
#' @param degs
#' @param allowed_regulator_ids
#'
get_gene_modules <- function(pb_cds, degs, allowed_regulator_ids) {

  all_degs = c(degs, allowed_regulator_ids) %>% unique
  all_gene_modules = monocle3::find_gene_modules(pb_cds[all_degs,], resolution=1e-2)
  all_gene_modules = left_join(all_gene_modules, rowData(pb_cds) %>% as.data.frame, by="id")
  return(all_gene_modules)
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
                                    knockout_genes = NULL,
                                    regulatory_genes = NULL,
                                    model_formula_str = "~cell_state",
                                    gene_module_df = NULL,
                                    sparsity_factor=1,
                                    pln_min_ratio=0.001,
                                    pln_num_penalties=30){

  # to do: handle knockouts?
  if (is.null(knockout_genes) == FALSE) {
    genes = union(genes, knockout_genes)
  }

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
  repressors = from_to_cond_diff %>% dplyr::filter(delta_log_abund < 0) %>% pull(cell_group)

  # get the subgraph induced by upregulated genes
  downreg_igraph = hooke:::return_igraph(model(gene_model_ccm))
  downreg_igraph = igraph::induced_subgraph(downreg_igraph, igraph::V(downreg_igraph)[repressors])
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



classify_genes_in_cell_state <- function(cell_state, state_graph, estimate_matrix, stderr_matrix, state_term="cell_group", log_fc_thresh=1, abs_expr_thresh = 1e-3, sig_thresh=0.05, cores=1){
  #expr_self = expr_mat[,cell_state]

  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))

  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))

  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))

  states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()

  expr_df = tibble(gene_id=row.names(estimate_matrix))

  message("      examining coeffficients", cell_state)

  expr_df$expr_self = pnorm(estimate_matrix[,cell_state] - log(abs_expr_thresh), sd = stderr_matrix[,cell_state], lower.tail=FALSE)
  expr_df$expr_self = p.adjust(expr_df$expr_self, method="BH") < sig_thresh

  expr_df$expressed_in_parents = NA
  expr_df$expressed_in_siblings = NA
  expr_df$higher_than_parents = NA
  expr_df$lower_than_parents = NA
  expr_df$higher_than_all_siblings = NA
  expr_df$lower_than_all_siblings = NA
  expr_df$higher_than_siblings = NA
  expr_df$lower_than_siblings = NA

  if (length(parents) > 0){
    expressed_in_parents_mat = pnorm(estimate_matrix[,parents, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,parents, drop=F], lower.tail=FALSE)
    expressed_in_parents_mat = apply(expressed_in_parents_mat, 2, p.adjust, method="BH")

    expressed_in_parents_mat = expressed_in_parents_mat < sig_thresh
    expr_df$expressed_in_parents = Matrix::rowSums(expressed_in_parents_mat) > 0

    higher_than_parents_stat = -t(sweep(t(estimate_matrix[,parents, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_parents_pval = pnorm(higher_than_parents_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_parents_pval = apply(higher_than_parents_pval, 2, p.adjust, method="BH")

    higher_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & higher_than_parents_pval < sig_thresh
    expr_df$higher_than_parents = Matrix::rowSums(higher_than_parents_mat) > 0

    lower_than_parents_pval = pnorm(-higher_than_parents_stat,
                                    sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_parents_pval = apply(lower_than_parents_pval, 2, p.adjust, method="BH")

    lower_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & lower_than_parents_pval < sig_thresh
    expr_df$lower_than_parents = Matrix::rowSums(lower_than_parents_mat) > 0
  }else{
    expr_df$expressed_in_parents = NA
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_parents = NA
    expr_df$lower_than_parents = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA

  }

  if (length(siblings) > 0){
    expressed_in_siblings_mat = pnorm(estimate_matrix[,siblings, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,siblings, drop=F], lower.tail=FALSE)
    expressed_in_siblings_mat = apply(expressed_in_siblings_mat, 2, p.adjust, method="BH")

    expressed_in_siblings_mat = expressed_in_siblings_mat < sig_thresh
    expr_df$expressed_in_siblings = Matrix::rowSums(expressed_in_siblings_mat) > 0

    higher_than_siblings_stat = -t(sweep(t(estimate_matrix[,siblings, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_siblings_pval = pnorm(higher_than_siblings_stat,
                                      sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_siblings_pval = apply(higher_than_siblings_pval, 2, p.adjust, method="BH")

    higher_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & higher_than_siblings_pval < sig_thresh
    expr_df$higher_than_all_siblings = Matrix::rowSums(higher_than_siblings_mat) == ncol(higher_than_siblings_pval)
    expr_df$higher_than_siblings = Matrix::rowSums(higher_than_siblings_mat) > 0

    lower_than_siblings_pval = pnorm(-higher_than_siblings_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_siblings_pval = apply(lower_than_siblings_pval, 2, p.adjust, method="BH")

    lower_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & lower_than_siblings_pval < sig_thresh
    expr_df$lower_than_all_siblings = Matrix::rowSums(lower_than_siblings_mat) == ncol(lower_than_siblings_mat)
    expr_df$lower_than_siblings = Matrix::rowSums(lower_than_siblings_mat) > 0


  }else{
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
  }

  expr_df = expr_df %>% tidyr::nest(data = !gene_id)

  message("      interpreting patterns")
  interpret_expression_pattern = function(pat_df){
    if (pat_df$expr_self){
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Maintained")
      }else if (pat_df$expressed_in_parents){
        # Expressed in self and parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self and parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Upregulated")
          else if(pat_df$lower_than_parents)
            return("Downregulated")
          else
            return("Maintained")
        } else {
          # Expressed in self and parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than siblings
              return("Specifically upregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than siblings
              return("Selectively upregulated")
            }
            else if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than siblings
              return("Upregulated")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
          else if(pat_df$lower_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than siblings
              return("Selectively downregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Selectively downregulated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than  all siblings
              return("Specifically downregulated")
            }
            else { # same as parent, same as siblings
              return("Downregulated")
            }
          }
          else { # same as parent
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Specifically maintained")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Selectively maintained")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than siblings
              return("Maintained")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
        }

      }else{
        # Expressed in self but not in parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self, not in parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Activated")
          else if(pat_df$lower_than_parents)
            return("Downregulated") # shouldn't happen
          else
            return("Maintained") # shouldn't happen
        } else {
          # Expressed in self and not in parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than all siblings
              return("Specifically activated")
            } else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than some siblings
              return("Selectively activated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Higher than parent, but lower than all siblings
              return("Activated")
            }
            if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than some siblings
              return("Activated")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
          else if(pat_df$lower_than_parents){
            # if the gene is lower in the parent, which is off, just mark the gene absent
            if (pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than all siblings
              return("Absent") # shouldn't happen
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Absent") # shouldn't happen
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
          else { # same as parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Absent")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Absent")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Same as parent, but lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
        }
      }
      return ("Expressed")
    }else{
      # Not expressed in self
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Absent")
      }else if (pat_df$expressed_in_parents){
        # Not expressed in self, but expressed in parents
        if (is.na(pat_df$expressed_in_siblings)){
          # Not expressed in self, expressed parent and there are no siblings
          if(pat_df$lower_than_parents)
            return("Deactivated")
          else
            return("Absent") # shouldn't happen
        } else {
          # Not expressed in self, expressed in parent and there are siblings
          if(pat_df$lower_than_parents){
            # Lower than parent
            if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than siblings
              return("Specifically deactivated")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than siblings
              return("Selectively deactivated")
            }
            return("Deactivated")
          }
          else {
            #Not expressed in self, not lower than parent
            return ("Absent")
          }
        }
      }else{
        # Not expressed in self or parents
        return ("Absent")
      }
      return ("Absent")
    }
    return ("Absent")
    #match_row = match(data.frame(t(pat_df)), data.frame(t(interp_table)))
    #interpetation[match_row]
  }
  #debug(interpret_expression_pattern)
  expr_df = expr_df %>% mutate(interpretation = purrr::map(.f = purrr::possibly(
    interpret_expression_pattern, NA_character_), .x = data))
  message("      completed", cell_state)
  return(expr_df)
}
#debug(classify_genes_in_cell_state)


#' Classify each gene's pattern of expresison in each state in a state transition graph
#' @export
classify_genes_over_graph <- function(ccm,
                                      state_graph,
                                      gene_ids = NULL,
                                      group_nodes_by=NULL,
                                      log_fc_thresh=1,
                                      abs_expr_thresh = 1e-3,
                                      sig_thresh=0.05,
                                      min_samples_detected = 2,
                                      min_cells_per_pseudobulk = 3,
                                      cores=1,
                                      ...){
  if (is.null(group_nodes_by)){
    pb_cds = pseudobulk_cds_for_states(wt_ccm_wl)
    state_term = "cell_group"
  }else{
    pb_cds = pseudobulk_cds_for_states(wt_ccm_wl, state_col = group_nodes_by)
    state_term = group_nodes_by
  }

  #cds_to_test = pb_cds[,as.character(colData(pb_cds)[,state_term]) %in% states_in_model]

  #colData(cds_to_test)[,state_term] = factor(as.character(colData(cds_to_test)[,state_term]), levels=states_in_model) # set the "self" state as the reference level

  #norm_expr_mat = normalized_counts(pb_cds, "size_only", pseudocount = 0)

  if (is.null(gene_ids) == FALSE){
    pb_cds = pb_cds[gene_ids,]
  }

  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds = pb_cds[genes_to_test,]

  pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group > min_cells_per_pseudobulk)

  message("fitting regression models")
  pb_cds = pb_cds[,pseudobulks_to_test]
  pb_group_models = fit_models(pb_cds,
                               model_formula_str=paste("~ 0 + ", state_term),
                               weights=colData(pb_cds)$num_cells_in_group,
                               cores=cores) %>% dplyr::select(gene_short_name, id, model, model_summary)

  message("      collecting coeffficients")
  pb_group_models = coefficient_table(pb_group_models) %>%
    dplyr::select(gene_short_name, id, term, estimate, std_err) %>%
    mutate(term = stringr::str_replace_all(term, state_term, ""))
  estimate_matrix = pb_group_models %>% dplyr::select(id, term, estimate)
  estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)

  gene_ids = estimate_matrix$id
  estimate_matrix$id = NULL
  estimate_matrix = as.matrix(estimate_matrix)
  row.names(estimate_matrix) = gene_ids
  colnames(estimate_matrix) = as.character(colnames(estimate_matrix))

  stderr_matrix = pb_group_models %>% dplyr::select(id, term, std_err)
  stderr_matrix = stderr_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)

  gene_ids = stderr_matrix$id
  stderr_matrix$id = NULL
  stderr_matrix = as.matrix(stderr_matrix)
  row.names(stderr_matrix) = gene_ids
  colnames(stderr_matrix) = as.character(colnames(stderr_matrix))

  #p_val_matrix = pnorm(estimate_matrix - log(abs_expr_thresh), sd = stderr_matrix, lower.tail=FALSE)

  #expr_thresh_mat = p_val_matrix < sig_thresh

  #cell_states = tibble(cell_state = unlist(igraph::V(state_graph)$name))
  states_to_assess = intersect(as.character(unique(colData(pb_cds)[,state_term])), unlist(igraph::V(state_graph)$name))
  cell_states = tibble(cell_state = states_to_assess)

  cell_states = cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
      classify_genes_in_cell_state, NA_real_), .x = cell_state,
      state_graph, estimate_matrix, stderr_matrix, state_term,
      log_fc_thresh=log_fc_thresh,
      abs_expr_thresh = abs_expr_thresh,
      sig_thresh=sig_thresh,
      cores=cores))

}

#' get the parent(s) of a state in a state transition graph
get_parents = function(state_graph, cell_state){
  parents = igraph::neighbors(state_graph, cell_state, mode="in")
  if (length(parents) > 0)
    return (parents$name)
  else
    return (c())
}

#' get the children of a state in a state transition graph
get_children = function(state_graph, cell_state){
  children = igraph::neighbors(state_graph, cell_state, mode="out")
  if (length(children) > 0)
    return (children$name)
  else
    return (c())
}

#' get the siblings of a state in a state transition graph
get_siblings = function(state_graph, cell_state){
  parents = get_parents(state_graph, cell_state)
  siblings = igraph::neighbors(state_graph, parents, mode="out")
  siblings = setdiff(siblings$name, cell_state) #exclude self
  return(siblings)
}

#' Compute a pseudobulk expression matrix for a model
pseudobulk_cds_for_states <- function(ccm, state_col=NULL, collapse_samples=FALSE){

  if (is.null(state_col)){
    cell_group_df = tibble::rownames_to_column(ccm@ccs@metadata[["cell_group_assignments"]])
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = cell_group)
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, "cell_group", sep="_")) %>% dplyr::select(rowname, pseudobulk_id, cell_group)
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, cell_group) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
    #%>% select(rowname, cell_group)
  }else{
    cell_group_df = tibble::rownames_to_column(ccm@ccs@metadata[["cell_group_assignments"]])
    cds_group_df = colData(ccm@ccs@cds) %>%
      as.data.frame %>% tibble::rownames_to_column() %>% dplyr::select(rowname, !!sym(state_col))
    cell_group_df = left_join(cell_group_df, cds_group_df, by=c("rowname"))
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = !!sym(state_col))
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, !!sym(state_col), sep="_")) %>% dplyr::select(rowname, pseudobulk_id, !!sym(state_col))
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, !!sym(state_col)) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
  }

  agg_expr_mat = monocle3::aggregate_gene_expression(wt_ccs@cds,
                                                     cell_group_df=cell_group_df,
                                                     norm_method="size_only",
                                                     scale_agg_values = FALSE,
                                                     pseudocount=0,
                                                     cell_agg_fun="mean")

  agg_expr_mat = agg_expr_mat[,agg_coldata$pseudobulk_id]

  row.names(agg_coldata) = agg_coldata$pseudobulk_id
  agg_coldata = agg_coldata[colnames(agg_expr_mat),]

  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccm@ccs@cds) %>% as.data.frame)
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  return(pseudobulk_cds)
}
