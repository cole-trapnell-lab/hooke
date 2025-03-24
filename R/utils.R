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
denylist <- function(edges) {
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
aggregated_expr_data <- function(cds, group_cells_by = "cell_type", gene_group_df=NULL, 
                                 gene_agg_fun = "sum", cell_agg_fun = "mean"){

  cds = cds[, !is.na(colData(cds)$timepoint)]
  cds = cds[, !is.na(colData(cds)[[group_cells_by]])]
  cds = cds[, colData(cds)[[group_cells_by]] != ""]

  cell_group_df <- data.frame(row.names = row.names(colData(cds)),
                              cell_id = row.names(colData(cds)))

  cell_group_df$cell_group <- colData(cds)[, group_cells_by]
  cell_group_df$cell_group <- as.character(cell_group_df$cell_group)
  cluster_binary_exprs = as.matrix(aggregate_gene_expression(cds,
                                                             cell_group_df = cell_group_df,
                                                             gene_group_df = gene_group_df,
                                                             norm_method = "binary",
                                                             scale_agg_values=FALSE))

  cluster_fraction_expressing_table = tibble::rownames_to_column(as.data.frame(cluster_binary_exprs))

  cluster_fraction_expressing_table = tidyr::gather(cluster_fraction_expressing_table,
                                                    "cell_group", "fraction_expressing", -rowname)

  cluster_mean_exprs = as.matrix(aggregate_gene_expression(cds,
                                                           cell_group_df = cell_group_df,
                                                           gene_group_df= gene_group_df,
                                                           norm_method = "size_only",
                                                           scale_agg_values=FALSE))
  rownames(cluster_mean_exprs) = rownames(cluster_binary_exprs)

  cluster_expr_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs))
  cluster_expr_table = tidyr::gather(cluster_expr_table, "cell_group",
                                     "mean_expression", -rowname)

  cluster_fraction_expressing_table$mean_expression = cluster_expr_table$mean_expression

  cluster_spec_mat = monocle3:::specificity_matrix(cluster_mean_exprs)
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
switch_umap_space <- function(cds, prefix = "subumap3d") {

  curr_umap_matrix = reducedDims(cds)[["UMAP"]]
  col_names = colnames(colData(cds))
  col_names = col_names[grepl(pattern = prefix, col_names)]

  umap_coords = colData(cds) %>%
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(col_names)) %>%
    as.matrix()

  reducedDims(cds)[["UMAP"]] = umap_coords[rownames(curr_umap_matrix),]
  return(cds)
}

#' this switches between global and sub umap space
#' @param cds
#' @param sub_space
#' @noRd
switch_reducedDims = function(cds, umap_space = c("global_UMAP", "sub_UMAP", "coemb_UMAP")) {
  # must put as as.matrix() otherwise can't build a nn index
  reducedDims(cds)[["PCA"]] = as.matrix(reducedDims(cds)[[gsub("UMAP", "PCA", umap_space)]])
  reducedDims(cds)[["Aligned"]] = as.matrix(reducedDims(cds)[[gsub("UMAP", "Aligned", umap_space)]])
  reducedDims(cds)[["UMAP"]] = as.matrix(reducedDims(cds)[[umap_space]])
  cds@metadata$umap_space = umap_space
  print(paste0("switching to ", umap_space, " space"))
  return(cds)
}

# check which umap space the cds currently is in
get_umap_space = function(cds) {
  return(cds@metadata$umap_space)
}


#' this switches the ccm umap space allowing you to plot contrast in
#' global and local spaces
#' @param ccm a cell count model object
#' @param umap_space
#' @noRd
switch_ccs_space <- function(ccs, umap_space = c("sub_UMAP", "global_UMAP")) {
  umap_space <- match.arg(umap_space)
  ccs@cds <- switch_reducedDims(ccs@cds, umap_space)
  ccs@cds_reduced_dims[["UMAP"]] <- ccs@cds_reduced_dims[[umap_space]]
  return(ccs)
}

#' this switches the ccm umap space allowing you to plot contrast in
#' global and local spaces
#' @param ccm a cell count model object
#' @param umap_space
#' @noRd
switch_ccm_space <- function(ccm, umap_space = c("sub_UMAP", "global_UMAP")) {
  umap_space <- match.arg(umap_space)
  ccm@ccs@cds <- switch_reducedDims(ccm@ccs@cds, umap_space)
  ccm@ccs@cds_reduced_dims[["UMAP"]] <- ccm@ccs@cds_reduced_dims[[umap_space]]
  return(ccm)
}


#' wrapper for plot contrast that allows you to facet the plot
#' you can also switch between sub + full umap space
#' @param ccm
#' @param cond_b_vs_a_tbl
#' @param cell_group
#' @param facet_group
#' @noRd
plot_sub_contrast = function(ccm,
                             cond_a_v_b_tbl,
                             selected_groups = NULL,
                             umap_space = "sub_UMAP",
                             facet_group = "assembly_group",
                             log_abundance_thresh = -5,
                             edge_size=2,
                             cell_size=1,
                             q_value_thresh = 1.0,
                             group_label_size=2,
                             plot_labels = c("significant", "all", "none"),
                             plot_edges = c("none", "all", "directed", "undirected"),
                             fc_limits=c(-3,3),
                             ...) {

  ccm = switch_ccm_space(ccm, umap_space = umap_space)

  # ccm@ccs = subset_ccs(ccm@ccs, ...)

  colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
  colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]
  ccm@ccs@cds_coldata$cell_group = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]

  cg_to_mg = colData(ccm@ccs@cds) %>%
              as.data.frame() %>%
              select(cell_group, facet_group) %>%
              distinct()

  cond_a_v_b_tbl = left_join(cond_a_v_b_tbl,
                             cg_to_mg, by = "cell_group")

  if (!is.null(selected_groups)){

    partition_cell_groups = cg_to_mg %>% filter(facet_group %in% selected_groups) %>% pull(cell_group)

    # ccm@ccs <- hooke:::subset_ccs(ccm@ccs, partition_cell_groups)
    ccm@ccs <- subset_ccs(ccm@ccs, cell_group %in% partition_cell_groups)
    cond_a_v_b_tbl = cond_a_v_b_tbl[cond_a_v_b_tbl$facet_group %in% selected_groups,]
  }


  plot_contrast(ccm,
                cond_a_v_b_tbl,
                edge_size = edge_size,
                cell_size = cell_size,
                q_value_thresh = q_value_thresh,
                group_label_size = group_label_size,
                plot_labels = plot_labels,
                fc_limits = fc_limits,
                plot_edges = plot_edges,
                ...) +
    facet_wrap(~facet_group)

}



plot_sub_contrast_3d = function(ccm,
                             cond_a_v_b_tbl,
                             selected_groups = NULL,
                             umap_space = "sub_UMAP",
                             facet_group = "assembly_group",
                             log_abundance_thresh = -5,
                             edge_size=2,
                             # cell_size=1,
                             q_value_thresh = 1.0,
                             group_label_size=2,
                             fc_limits=c(-3,3),
                             ...) {

  ccm = switch_ccm_space(ccm, umap_space = umap_space)

  # ccm@ccs = subset_ccs(ccm@ccs, ...)

  colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
  colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]
  ccm@ccs@cds_coldata$cell_group = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]

  cg_to_mg = colData(ccm@ccs@cds) %>%
    as.data.frame() %>%
    select(cell_group, facet_group) %>%
    distinct()

  cond_a_v_b_tbl = left_join(cond_a_v_b_tbl,
                             cg_to_mg, by = "cell_group")

  if (!is.null(selected_groups)){

    partition_cell_groups = cg_to_mg %>% filter(facet_group %in% selected_groups) %>% pull(cell_group)

    # ccm@ccs <- hooke:::subset_ccs(ccm@ccs, partition_cell_groups)
    ccm@ccs <- subset_ccs(ccm@ccs, cell_group %in% partition_cell_groups)
    cond_a_v_b_tbl = cond_a_v_b_tbl[cond_a_v_b_tbl$facet_group %in% selected_groups,]
  }


  plot_contrast_3d(ccm,
                  cond_a_v_b_tbl,
                  edge_size = edge_size,
                  # cell_size = cell_size,
                  q_value_thresh = q_value_thresh,
                  group_label_size = group_label_size,
                  fc_limits = fc_limits,
                  ...)

}


plot_sub_abundance = function(ccs,
                             cond_a_v_b_tbl,
                             selected_groups = NULL,
                             umap_space = "sub_UMAP",
                             facet_group = "assembly_group",
                             log_abundance_thresh = -5,
                             edge_size=2,
                             cell_size=1,
                             q_value_thresh = 1.0,
                             group_label_size=2,
                             plot_labels = c("significant", "all", "none"),
                             plot_edges = c("none", "all", "directed", "undirected"),
                             fc_limits=c(-3,3),
                             nrow = NULL,
                             ncol = NULL,
                             ...) {

  ccs = switch_ccs_space(ccs, umap_space = umap_space)
  plot_labels <- match.arg(plot_labels)

  # ccm@ccs = subset_ccs(ccm@ccs, ...)

  colData(ccs@cds)[["cell_group"]] = colData(ccs@cds)[[ccs@info$cell_group]]
  colData(ccs@cds)[["facet_group"]] = colData(ccs@cds)[[facet_group]]
  ccs@cds_coldata$cell_group = colData(ccs@cds)[[ccs@info$cell_group]]

  cg_to_mg = colData(ccs@cds) %>%
    as.data.frame() %>%
    select(cell_group, facet_group) %>%
    distinct()

  cond_a_v_b_tbl = left_join(cond_a_v_b_tbl,
                             cg_to_mg, by = "cell_group")

  if (!is.null(selected_groups)){

    partition_cell_groups = cg_to_mg %>% filter(facet_group %in% selected_groups) %>% pull(cell_group)

    # ccm@ccs <- hooke:::subset_ccs(ccm@ccs, partition_cell_groups)
    ccs <- subset_ccs(ccs, cell_group %in% partition_cell_groups)
    cond_a_v_b_tbl = cond_a_v_b_tbl[cond_a_v_b_tbl$facet_group %in% selected_groups,]
  }

  plot_abundance(ccs,
                cond_a_v_b_tbl,
                #edge_size = edge_size,
                cell_size = cell_size,
                q_value_thresh = q_value_thresh,
                group_label_size = group_label_size,
                plot_labels = plot_labels,
                fc_limits = fc_limits,
                #plot_edges = plot_edges,
                ...) +
    facet_wrap(~facet_group, nrow = nrow, ncol = ncol)

}


plot_sub_contrast_2 = function(ccm,
                               cond_a_v_b_tbl,
                               umap_space = "sub_UMAP",
                               # log_abundance_thresh = -5,
                               # edge_size=2,
                               # cell_size=1,
                               q_value_thresh = 1.0,
                               # group_label_size=2,
                               plot_labels = c("significant", "all", "none"),
                               plot_edges = c("none", "all", "directed", "undirected"),
                               # fc_limits=c(-3,3),
                               # downsample = 1e5,
                               ...) {

  # plot_labels = match.arg(plot_labels)
  # plot_edges = match.arg(plot_edges)

  ccm = switch_ccm_space(ccm, umap_space = umap_space)

  ccm@ccs = subset_ccs(ccm@ccs, ...)

  plot_contrast(ccm,
                cond_a_v_b_tbl,
                # edge_size = edge_size,
                # cell_size = cell_size,
                q_value_thresh = q_value_thresh,
                # group_label_size = group_label_size,
                # plot_labels = plot_labels,
                # fc_limits = fc_limits,
                # plot_edges = plot_edges,
                downsample = downsample
                )

}


#' #' wrapper for plot contrast that allows you to facet the plot
#' #' you can also switch between sub + full umap space
#' #' @param ccm
#' #' @param cond_b_vs_a_tbl
#' #' @param cell_group
#' #' @param facet_group
#' #' @noRd
#' #'
#' plot_sub_contrast <- function (ccm,
#'                                cond_b_vs_a_tbl,
#'                                facet_group = "major_group",
#'                                select_group = NULL,
#'                                log_abundance_thresh = -5,
#'                                edge_size=2,
#'                                cell_size=1,
#'                                q_value_thresh = 1.0,
#'                                group_label_size=2,
#'                                plot_labels = c("significant", "all", "none"),
#'                                plot_edges = c("all", "directed", "undirected", "none"),
#'                                fc_limits=c(-3,3),
#'                                prefix = "subumap3d",
#'                                ...) {
#'
#'
#'   ccm@ccs@cds = switch_umap_space(ccm@ccs@cds, prefix = prefix)
#'
#'   colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
#'   colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]
#'
#'   cg_to_mg = as.data.frame(colData(ccm@ccs@cds)) %>%
#'     select("cell_group", "facet_group") %>%
#'     distinct()
#'
#'   cond_b_vs_a_tbl = cond_b_vs_a_tbl %>%
#'     left_join(cg_to_mg, by = "cell_group")
#'
#'   if (is.null(select_group) == FALSE) {
#'     ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[["facet_group"]] == select_group]
#'     cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% filter(facet_group == select_group)
#'
#'     sub_cell_groups = unique(colData(ccm@ccs@cds)$cell_group)
#'     ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]
#'
#'   }
#'
#'
#'   plot_contrast(ccm,
#'                 cond_b_vs_a_tbl,
#'                 edge_size=edge_size,
#'                 cell_size=cell_size,
#'                 q_value_thresh = q_value_thresh,
#'                 group_label_size=group_label_size,
#'                 plot_labels = plot_labels,
#'                 fc_limits=fc_limits,
#'                 plot_edges = plot_edges,
#'                 ...) +
#'     facet_wrap(~facet_group)
#'
#' }


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
                                allowlist = NULL,
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
                                                       allowlist = allowlist
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

#' function to easily plot to easily feed in normalized counts from ccs
#' into boxplots
#' @param ccs
get_norm_df = function(ccs) {
  
  cell_type_total <- Matrix::colSums(counts(ccs))
  geometric_mean = exp(mean(log(cell_type_total)))
  
  get_norm_counts(ccs) %>%
    as.data.frame() %>%
    rownames_to_column("cell_group") %>%
    pivot_longer(-cell_group,
                 names_to = "sample",
                 values_to = "count") %>%
    left_join(colData(ccs) %>% as.data.frame, by = "sample") %>% 
    mutate(count_per_1000 = count*1000/geometric_mean)

}



fill_missing_terms_with_default_values = function(ccm, newdata, pln_model = c("full", "reduced"), verbose=FALSE){
  pln_model <- match.arg(pln_model)

  # check that all terms in new data have been specified
  if (pln_model == "reduced")
    missing_terms = setdiff(names(ccm@model_aux[["reduced_model_xlevels"]]), names(newdata))
  else if (pln_model == "full")
    missing_terms = setdiff(names(ccm@model_aux[["full_model_xlevels"]]), names(newdata))

  if (length(missing_terms) >= 1) {

    default_df = lapply(missing_terms, function(term){
      df = data.frame(t = levels(factor(colData(ccm@ccs)[[term]]))[1])
      names(df) = term
      df
    }) %>% bind_cols()

    newdata = cbind(newdata, tibble(default_df))

    if (verbose){
      print( paste0(paste(missing_terms,collapse = ", "),
                    " missing from specified newdata columns. Assuming default values: ",
                    paste(default_df[1,],collapse = ", ")))
    }
  }
  return (newdata)
}

