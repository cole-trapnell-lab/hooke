#' add umap coords to a data frame
#' @noRd
#' @import dplyr
add_umap_coords <- function(df, umap_centers) {

  from_df = left_join(df, umap_centers, by = c("from"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_from_", .x))

  to_df = left_join(df, umap_centers, by = c("to"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_to_", .x))

  return(full_join(from_df, to_df))

}

#' return the complementing edges
#' @noRd
#' @importFrom igraph graph_from_data_frame
blacklist <- function(edges) {
  igraph::graph_from_data_frame(edges) %>%
    igraph::complementer() %>%
    igraph::as_data_frame()
}

#' @noRd
#' returns edges based on pcor values
get_pcor_edges <- function(ccm, selected_model=c("reduced", "full")) {
  model(ccm, selected_model)$latent_network() %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    pivot_longer(-c("from"), names_to = "to") %>%
    filter(from!=to, value!=0)
}

#' @noRd
#' returns pairwise distances
get_distances <- function(ccs, method="euclidean", matrix=T) {
  umap_centers = centroids(ccs)
  row.names(umap_centers) <- umap_centers$cell_group

  dist_matrix = dist(umap_centers[,-1], method = method, upper=T, diag = T)

  if (matrix) {
    return(dist_matrix)
  } else {
    dist_df =  dist_matrix %>%
      as.matrix() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("from") %>%
      tidyr::pivot_longer(-from, names_to = "to", values_to = "dist")
    return(dist_df)
  }
}






#' calculates the fraction expression, mean expression and specificity
#' for each gene split by the specified group
#' @param cds
#' @param group_cells_by
#' @noRd
#'
aggregated_expr_data <- function(cds, group_cells_by = "cell_type_broad"){

  cds = cds[, !is.na(colData(cds)$timepoint)]
  cds = cds[, !is.na(colData(cds)[[group_cells_by]])]
  cds = cds[, colData(cds)[[group_cells_by]] != ""]

  cell_group_df <- data.frame(row.names = row.names(colData(cds)),
                              cell_id = row.names(colData(cds)))

  cell_group_df$cell_group <- colData(cds)[, group_cells_by]
  cell_group_df$cell_group <- as.character(cell_group_df$cell_group)
  cluster_binary_exprs = as.matrix(aggregate_gene_expression(cds,
                                                             cell_group_df = cell_group_df,
                                                             norm_method = "binary",
                                                             scale_agg_values=FALSE))

  cluster_fraction_expressing_table = tibble::rownames_to_column(as.data.frame(cluster_binary_exprs))

  cluster_fraction_expressing_table = tidyr::gather(cluster_fraction_expressing_table,
                                                    "cell_group", "fraction_expressing", -rowname)

  cluster_mean_exprs = as.matrix(aggregate_gene_expression(cds,
                                                           cell_group_df = cell_group_df,
                                                           norm_method = "size_only",
                                                           scale_agg_values=FALSE))

  cluster_expr_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs))
  cluster_expr_table = tidyr::gather(cluster_expr_table, "cell_group",
                                     "mean_expression", -rowname)

  cluster_fraction_expressing_table$mean_expression = cluster_expr_table$mean_expression

  cluster_spec_mat = monocle3:::specificity_matrix(cluster_mean_exprs, cores = 4)
  cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
  cluster_spec_table = tidyr::gather(cluster_spec_table, "cell_group",
                                     "specificity", -rowname)

  cluster_fraction_expressing_table$specificity = cluster_spec_table$specificity
  cluster_fraction_expressing_table = cluster_fraction_expressing_table %>%
    dplyr::rename("gene_id" = rowname) %>%
    dplyr::left_join(rowData(cds) %>%
                       as.data.frame() %>%
                       dplyr::select("gene_id" = "id", gene_short_name),
                     by = "gene_id") %>%
    dplyr::select(cell_group, gene_id, gene_short_name, everything())

  return(cluster_fraction_expressing_table)

}


#' this switches the umap space to the sub umap space
#' @param cds
#' @param sub_space
#' @noRd
switch_umap_space <- function(cds, sub_space = TRUE) {
  curr_umap_matrix = reducedDims(cds)[["UMAP"]]

  if (sub_space) {
    umap_coords = colData(cds) %>%
      as.data.frame() %>%
      select(subumap3d_1, subumap3d_2, subumap3d_3) %>%
      as.matrix()
  } else {
    umap_coords = colData(cds) %>%
      as.data.frame() %>%
      select(umap3d_1, umap3d_2, umap3d_3) %>%
      as.matrix()
  }
  reducedDims(cds)[["UMAP"]] = umap_coords[rownames(curr_umap_matrix),]
  return(cds)
}




#' wrapper for plot contrast that allows you to facet the plot
#' you can also switch between sub + full umap space
#' @param ccm
#' @param cond_b_vs_a_tbl
#' @param cell_group
#' @param facet_group
#' @noRd
#'
plot_sub_contrast <- function (ccm,
                               cond_b_vs_a_tbl,
                               facet_group = "major_group",
                               select_group = NULL,
                               log_abundance_thresh = -5,
                               edge_size=2,
                               cell_size=1,
                               q_value_thresh = 1.0,
                               group_label_size=2,
                               plot_labels = c("significant", "all", "none"),
                               plot_edges = c("all", "directed", "undirected", "none"),
                               fc_limits=c(-3,3),
                               sub_space = T,
                               ...) {


  ccm@ccs@cds = switch_umap_space(ccm@ccs@cds, sub_space = sub_space)

  colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
  colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]

  cg_to_mg = as.data.frame(colData(ccm@ccs@cds)) %>%
    select("cell_group", "facet_group") %>%
    distinct()

  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>%
    left_join(cg_to_mg, by = "cell_group")

  if (is.null(select_group) == FALSE) {
    ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[["facet_group"]] == select_group]
    cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% filter(facet_group == select_group)

    sub_cell_groups = unique(colData(ccm@ccs@cds)$cell_group)
    ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]

  }


  plot_contrast(ccm,
                cond_b_vs_a_tbl,
                edge_size=edge_size,
                cell_size=cell_size,
                q_value_thresh = q_value_thresh,
                group_label_size=group_label_size,
                plot_labels = plot_labels,
                fc_limits=fc_limits,
                plot_edges = plot_edges,
                ...) +
    facet_wrap(~facet_group)

}




#' @noRd
convert_cluster_to_cell_group = function(ccs, cell_group) {

  colData(ccs@cds)[["cell_group"]] = colData(ccs@cds)[[cell_group]]

  cluster_to_cellgroup = colData(ccs@cds) %>%
    as.data.frame() %>%
    group_by(cluster, cell_group) %>%
    tally() %>%
    top_n(1)

  return(cluster_to_cellgroup)

}

#' @noRd
convert_to_col = function(ccs, df, colname) {

  df %>%
    left_join(convert_cluster_to_cell_group(ccs, colname),
              by = c("cell_group" = "cluster")) %>%
    select(-cell_group, -n) %>%
    rename("cell_group" = cell_group.y)

}


#'
#' @param perturbation of interest
#' @param ccs a cell count set object
#' @param ctrl_ids
#' @param col_name column name to select for perturbations
#' @param interaction boolean whether to include an additive model or interaction term
#' @noRd
#'
fit_perturb_ccm = function(perturbation,
                                ccs,
                                ctrl_ids = c("ctrl-uninj", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                                col_name = "gene_target",
                                num_time_breaks = 3,
                                sparsity_factor = 0.2,
                                main_model_formula_string = NULL,
                                nuisance_model_formula_string = NULL,
                                whitelist = NULL,
                                interaction = F){

  subset_ccs = ccs[, !is.na(colData(ccs)[[col_name]]) ]
  subset_ccs = subset_ccs[,colData(subset_ccs)[[col_name]] == perturbation | colData(subset_ccs)[[col_name]] %in% ctrl_ids]
  colData(subset_ccs)$knockout = colData(subset_ccs)[[col_name]] == perturbation

  # subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]
  # colData(subset_ccs)$knockout = colData(subset_ccs)$gene_target == genotype

  knockout_time_start = min(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  knockout_time_stop = max(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  subset_ccs = subset_ccs[,colData(subset_ccs)$timepoint >= knockout_time_start & colData(subset_ccs)$timepoint <= knockout_time_stop]
  time_breakpoints = c()

  if (num_time_breaks > 2 & knockout_time_stop > knockout_time_start){
    time_breakpoints = seq(knockout_time_start, knockout_time_stop, length.out=num_time_breaks)
    time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
    main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")
  }
  else if (num_time_breaks == 2){
    colData(subset_ccs)$timepoint = as.factor(colData(subset_ccs)$timepoint )
    main_model_formula_str = "~ timepoint"
  }
  else{
    main_model_formula_str = ""
  }

  if (length(unique(colData(subset_ccs)$expt)) > 1)
      nuisance_model_formula_str = "~ expt"
  else
      nuisance_model_formula_str = "~ 1"

  if (interaction) {
    main_model_formula_str = paste0(main_model_formula_str, " * knockout")
  } else {
    main_model_formula_str = paste0(main_model_formula_str, " + knockout")
  }

  # print(main_model_formula_str)

  # if predefined overwrite
  if (!is.null(main_model_formula_string)) {
    main_model_formula_str = main_model_formula_string
  }

  if (!is.null(nuisance_model_formula_string)) {
    nuisance_model_formula_str = nuisance_model_formula_string
  }

  # print(main_model_formula_str)
  # print(nuisance_model_formula_str)

  perturb_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
                                                       main_model_formula_str = main_model_formula_str,
                                                       nuisance_model_formula_str = nuisance_model_formula_str,
                                                       whitelist = whitelist
  ))

  perturb_ccm = select_model(perturb_ccm, sparsity_factor = sparsity_factor)
  return(perturb_ccm)
}

#' @noRd
collect_genotype_effects = function(ccm, timepoint=24, expt="GAP16"){
  control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=timepoint, expt=expt))
  knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=timepoint, expt=expt))
  genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
}



#' filters a cds
#' @param cds
#' @param ... expressions that return a logical value
#' @noRd
#' @return a cell data set object
filter_cds <- function(cds, ...) {
  cell_names = as.data.frame(cds@colData) %>% filter(...) %>% rownames()
  new_cds = cds[, cell_names] %>%
            estimate_size_factors() %>%
            detect_genes()
  return(cds[, cell_names])
}

#'
#' @param ccs
#' @param ... expressions that return a logical value
#' @noRd
#' @return a cell count set object
filter_ccs <- function(ccs,
                       # recompute_sf = TRUE,
                       ...) {

  ccs@cds = filter_cds(ccs@cds, ...)

  ccs = new_cell_count_set(ccs@cds,
                           sample_group = ccs@info$sample_group,
                           cell_group = ccs@info$cell_group)

  # if (recompute_sf) {
  #   ccs = new_cell_count_set(ccs@cds,
  #                            sample_group = ccs@info$sample_group,
  #                            cell_group = ccs@info$cell_group)
  # } else {
  #   cds_version = ccs@metadata$cds_version
  #   cell_group_assignments = ccs@metadata[["cell_group_assignments"]][colnames(ccs@cds),]
  #   ccs@metadata = list(cds_version = cds_version,
  #                       cell_group_assignments = cell_group_assignments)
  # }

  return(ccs)
}

#'
#' @param ccm
#' @param ... expressions that return a logical value
#' @noRd
#' @return a cell count model object
filter_ccm <- function(ccm, ...) {
  ccm@ccs = filter_ccs(ccm@ccs, ...)
  return(ccm )
}

#' change cell count set grouping
#' @param ccm a cell count model
#' @param cell_group string specifying how to aggregate the cell counts
#' @noRd
contract_ccm <- function(ccm, cell_group = NULL, sample_group = NULL) {

  if (is.null(sample_group)) {
    sample_group = ccm@ccs@info$sample_group
  }

  if (is.null(cell_group)) {
    cell_group = ccm@ccs@info$cell_group
  }

  ccm@ccs = new_cell_count_set(ccm@ccs@cds,
                               sample_group = sample_group,
                               cell_group = cell_group)
  return(ccm)
}






#' score DEGs on specificity
#' @param gene_patterns_over_cell_graph
#' @param pattern
#' @param cell_group
#' @noRd
#'
top_gene_pattern <- function(ccm,
                             gene_patterns_over_cell_graph,
                             pattern = NULL,
                             grep_pattern = NULL,
                             cell_group = "cell_type_sub",
                             n = 5,
                             return_genes = T) {

  if (is.null(pattern) == FALSE) {
    # use 1 type of pattern
    all_genes = gene_patterns_over_cell_graph %>%
      filter(interpretation == pattern) %>%
      pull(gene_short_name) %>%
      unique

  }
  else if (is.null(grep_pattern) == FALSE) {
    all_genes = gene_patterns_over_cell_graph %>%
      filter(grepl(grep_pattern, interpretation, ignore.case = T)) %>%
      pull(gene_short_name) %>%
      unique
  }
  else {
    # use all genes
    all_genes = gene_patterns_over_cell_graph %>%
      pull(gene_short_name) %>%
      unique
  }

  if (length(all_genes) <= 5) {
    return(all_genes)
  }

  all_genes_marker_scores = top_markers(ccm@ccs@cds[rowData(ccm@ccs@cds)$gene_short_name %in% all_genes,],
                                        group_cells_by = cell_group)

  top_genes = all_genes_marker_scores %>%
    as_tibble() %>% filter(marker_test_q_value < 0.01) %>%
    group_by(cell_group) %>% slice_max(marker_score, n=5)

  if (return_genes) {
    top_genes = top_genes %>% pull(gene_short_name)
    return(unique(top_genes))
  } else {
    # top_genes = top_genes %>% select(cell_group, gene_short_name)
    return(top_genes %>% distinct)
  }

}


#' @noRd
single_thread_omp = function(){
  old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
  if (is.na(old_omp_num_threads)){
    old_omp_num_threads = 1
  }
  RhpcBLASctl::omp_set_num_threads(1)
  return (old_omp_num_threads)
}

single_thread_blas = function(){
  old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))
  if (is.na(old_blas_num_threads)){
    old_blas_num_threads = 1
  }
  RhpcBLASctl::blas_set_num_threads(1)
  return (old_blas_num_threads)
}






