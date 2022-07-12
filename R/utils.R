#' adds umap coords to a data frame
#' @import dplyr
add_umap_coords <- function(df, umap_centers) {

  from_df = left_join(df, umap_centers, by = c("from"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_from_", .x))

  to_df = left_join(df, umap_centers, by = c("to"= "cell_group")) %>%
    rename_with(~gsub(pattern="_",replacement = "_to_", .x))

  return(full_join(from_df, to_df))

}

#' return the complementing edges
#' @importFrom igraph graph_from_data_frame
blacklist <- function(edges) {
  igraph::graph_from_data_frame(edges) %>%
    igraph::complementer() %>%
    igraph::as_data_frame()
}

#' returns edges based on pcor values
get_pcor_edges <- function(ccm, selected_model=c("reduced", "full")) {
  model(ccm, selected_model)$latent_network() %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(-c("from"), names_to = "to") %>%
    filter(from!=to, value!=0)
}

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
      rownames_to_column("from") %>%
      tidyr::pivot_longer(-from, names_to = "to", values_to = "dist")
    return(dist_df)
  }
}


add_covariate <- function(ccs, pb_cds, covariate) {
  assertthat::assert_that(
    tryCatch(expr = covariate %in% colnames(colData(ccs@cds)),
             error = function(e) FALSE),
    msg = "covariate not in colnames")

  assertthat::assert_that(
    tryCatch(expr = !covariate %in% colnames(colData(pb_cds)),
             error = function(e) FALSE),
    msg = "covariate already in colnames")

  group_to_covariate = colData(ccs@cds) %>%
    as.data.frame %>%
    select(group_id, all_of(covariate)) %>%
    distinct()

  pb_coldata = colData(pb_cds) %>%
    as.data.frame %>%
    left_join(group_to_covariate, by = "group_id")

  colData(pb_cds)[[covariate]] =  pb_coldata[[covariate]]

  return(pb_cds)
}

find_edges <- function(states, gene_model_ccm) {

  from_cond_est = estimate_abundances(gene_model_ccm, tibble("cell_state" = states[1]))
  to_cond_est = estimate_abundances(gene_model_ccm, tibble("cell_state" = states[2]))
  from_to_cond_diff = compare_abundances(gene_model_ccm, from_cond_est, to_cond_est) %>%
    filter(delta_p_value < 0.05)

  gene_edges = collect_pln_graph_edges(gene_model_ccm, from_to_cond_diff)
  gene_edges

}

#' calculates the fraction expression, mean expression and specificity
#' for each gene split by the specified group
#' @param cds
#' @param group_cells_by
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



#' projection wrap up
#' @param query_cds
#' @param ref_cds
#' @param directory_path a string giving the name of the directory from which to read the model files
#' @param transfer_type
run_projection <- function(query_cds,
                           ref_cds,
                           directory_path,
                           transfer_type = "cluster") {

  query_cds <- load_transform_models(query_cds, directory_path)
  query_cds <- preprocess_transform(query_cds, method="PCA")
  query_cds <- align_beta_transform(query_cds)
  query_cds <- reduce_dimension_transform(query_cds, method="UMAP")

  umap_dims = reducedDims(query_cds)[["UMAP"]]
  annoy_res = uwot:::annoy_search(umap_dims,
                                  k = 10,
                                  ann = query_cds@reduce_dim_aux$UMAP$nn_index$annoy_index)

  labels  = get_nn_labels(umap_dims,
                          annoy_res,
                          as.data.frame(colData(ref_cds)),
                          transfer_type = transfer_type)

  query_cds = query_cds[,!is.na(colData(query_cds)[[transfer_type]])]

  return(query_cds)
}

convert_cluster_to_cell_group = function(ccs, cell_group) {

  colData(ccs@cds)[["cell_group"]] = colData(ccs@cds)[[cell_group]]

  cluster_to_cellgroup = colData(ccs@cds) %>%
    as.data.frame() %>%
    group_by(cluster, cell_group) %>%
    tally() %>%
    top_n(1)

  return(cluster_to_cellgroup)

}

convert_to_col = function(ccs, df, colname) {

  df %>%
    left_join(convert_cluster_to_cell_group(ccs, colname),
              by = c("cell_group" = "cluster")) %>%
    select(-cell_group, -n) %>%
    rename("cell_group" = cell_group.y)

}


fit_genotype_ccm = function(genotype,
                            ccs,
                            # ctrl_ids=c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                            ctrl_ids = c("ctrl-uninj", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                            num_time_breaks = 3,
                            sparsity_factor = 0.2,
                            main_model_formula_string = NULL,
                            nuisance_model_formula_string = NULL,
                            whitelist = NULL,
                            multiply = F){

  subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]

  colData(subset_ccs)$knockout = colData(subset_ccs)$gene_target == genotype
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
    main_model_formula_str = "~ as.factor(timepoint)"
  }
  else{
    main_model_formula_str = ""
  }

  if (multiply) {
    nuisance_model_formula_str = main_model_formula_str

    if (length(unique(colData(subset_ccs)$expt)) > 1) {
      nuisance_model_formula_str = paste0(nuisance_model_formula_str, "+ expt")
    }
  } else {
    if (length(unique(colData(subset_ccs)$expt)) > 1)
      nuisance_model_formula_str = "~ expt"
    else
      nuisance_model_formula_str = "~ 1"
  }

  if (multiply) {
    main_model_formula_str = paste0(main_model_formula_str, "*knockout")
  }
  main_model_formula_str = paste0(main_model_formula_str, " + knockout")


  # if predefined overwrite
  if (!is.null(main_model_formula_string)) {
    main_model_formula_str = main_model_formula_string
  }

  if (!is.null(nuisance_model_formula_string)) {
    nuisance_model_formula_str = nuisance_model_formula_string
  }

  # print(main_model_formula_str)
  # print(nuisance_model_formula_str)

  genotype_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
                                                       main_model_formula_str = main_model_formula_str,
                                                       nuisance_model_formula_str = nuisance_model_formula_str,
                                                       whitelist = whitelist
  ))

  genotype_ccm = select_model(genotype_ccm, sparsity_factor = sparsity_factor)
  return(genotype_ccm)
}

collect_genotype_effects = function(ccm, timepoint=24, expt="GAP16"){
  control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=timepoint, expt=expt))
  knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=timepoint, expt=expt))
  genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
}


# contract_ccm <- function(ccm, group_nodes_by = "cell_type_broad") {
#
#   cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
#
#   cell_group_metadata = colData(ccm@ccs@cds) %>%
#     as.data.frame %>% select(!!sym(group_nodes_by))
#
#   cell_group_metadata$cell_group = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)
#
#   group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
#     as.data.frame %>%
#     dplyr::count(cell_group, !!sym(group_nodes_by)) %>%
#     dplyr::group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
#   colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
#
#   new_cell_groups = unique(group_by_metadata$group_nodes_by)
#
#   ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]] %>%
#     left_join(group_by_metadata, by = "cell_group") %>%
#     select(-cell_group) %>%
#     dplyr::rename("cell_group" = group_nodes_by)
#
#   return(ccm)
#
# }

#' filters a cds
#' @param cds
filter_cds <- function(cds, ...) {
  cell_names = as.data.frame(cds@colData) %>% filter(...) %>% rownames()
  return(cds[, cell_names])
}

#'
#' @param ccs
filter_ccs <- function(ccs, ...) {
  ccs@cds = filter_cds(ccs@cds, ...)
  cell_groups = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  ccs@metadata = ccs@metadata[["cell_group_assignments"]][colnames(ccs@cds),]
  return(ccs)
}

#'
#' @param ccm
filter_ccm <- function(ccm, ...) {
  ccm@ccs = filter_ccs(ccm@ccs, ...)
  return(ccm )
}

#' change cell group
#' @param ccm
#' @param cell_group
contract_ccm <- function(ccm, cell_group) {
  ccm@ccs = new_cell_count_set(ccm@ccs@cds,
                                                sample_group = ccm@ccs@info$sample_group,
                                                cell_group = cell_group)
  return(ccm)
}




subset_gap <- function(cds, major_group) {

  sub_cds = cds[, colData(cds)$major_group == major_group]

  reducedDims(sub_cds)[["UMAP"]] = sub_cds@colData %>% as.data.frame %>%
    select(subumap3d_1, subumap3d_2, subumap3d_3) %>%
    as.matrix

  return(sub_cds)
}

# subset_ccs = function(ccs, col_name, col_value) {
#
#   ccs = ccs[, colData(ccs)[[col_name]] == col_value]
#   ccs@cds = ccs@cds[, ccs@cds@colData[[col_name]] == col_value]
#
#   ccs@metadata$cell_group_assignments = ccs@metadata$cell_group_assignments[colnames(ccs@cds),]
#   return(ccs)
# }
#
# subset_ccm = function(ccm, col_name, col_values) {
#
#   ccs = ccm@ccs
#   ccs = ccs[, colData(ccs)[[col_name]] %in% col_values]
#   ccs@cds = ccs@cds[, ccs@cds@colData[[col_name]] %in% col_values]
#
#   ccs@metadata$cell_group_assignments = ccs@metadata$cell_group_assignments[colnames(ccs@cds),]
#   ccm@ccs = ccs
#   return(ccm)
# }


threshold_expression_matrix <- function(norm_expr_mat, relative_expr_thresh = 0.25, abs_expr_thresh = 1e-3, scale_tpc=1e6){
  norm_expr_mat = norm_expr_mat %>% as.matrix
  norm_expr_mat = norm_expr_mat #* scale_tpc
  expression_max_values = norm_expr_mat %>% matrixStats::rowMaxs()
  names(expression_max_values) = row.names(norm_expr_mat)
  dynamic_range = expression_max_values # size of dynamic range of expression for each gene, in logs (assuming min is zero)

  dynamic_range_thresh = relative_expr_thresh * dynamic_range
  #abs_thresh = rep(abs_expr_thresh, nrow(norm_expr_mat))

  expression_thresh_vec = dynamic_range_thresh
  expression_thresh_vec[expression_thresh_vec < abs_expr_thresh] = Inf

  names(expression_thresh_vec) = row.names(norm_expr_mat)
  expr_over_thresh = norm_expr_mat > expression_thresh_vec
  return(expr_over_thresh)
}


#' score DEGs on specificity
#' @param gene_patterns_over_cell_graph
#' @param pattern
#' @param cell_group
#'
top_gene_pattern <- function(ccm,
                             gene_patterns_over_cell_graph,
                             pattern = NULL,
                             grep_pattern = NULL,
                             cell_group = "cell_type_sub",
                             n = 5) {

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
    group_by(cell_group) %>% slice_max(marker_score, n=5) %>% pull(gene_short_name)

  return(unique(top_genes))
}


plot_single_gene <- function(cds, gene, file_path = "./", x=1, y=3) {

  plot_cells(cds, x = x, y = y,  genes = c(gene), label_cell_groups = F, show_trajectory_graph = F,
             cell_size = 0.5, cell_stroke = 0, alpha = 0.8) +
    scale_color_viridis_c() +
    theme_void() +
    theme(legend.position = "none")

  ggsave(paste0(file_path, gene, ".png"),
         dpi = 750,
         height = 2,
         width = 2.2,
         bg = "transparent")
}
