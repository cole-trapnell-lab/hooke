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
get_pcor_edges <- function(ccm) {
  ccm@best_model$latent_network() %>%
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
      pivot_longer(-from, names_to = "to", values_to = "dist")
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
                               # cell_group = "cell_type_broad",
                               facet_group = "major_group",
                               select_group = NULL,
                               log_abundance_thresh = -5,
                               scale_shifts_by=c("receiver", "sender", "none"),
                               edge_size=2,
                               cell_size=1,
                               q_value_thresh = 1.0,
                               group_label_size=2,
                               plot_labels = c("significant", "all", "none"),
                               fc_limits=c(-3,3),
                               sender_cell_groups=NULL,
                               receiver_cell_groups=NULL,
                               plot_edges = TRUE,
                               sub_space = T) {


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

  gp = plot_contrast(ccm,
                cond_b_vs_a_tbl,
                scale_shifts_by=scale_shifts_by,
                edge_size=edge_size,
                cell_size=cell_size,
                q_value_thresh = q_value_thresh,
                group_label_size=group_label_size,
                plot_labels = plot_labels,
                fc_limits=fc_limits,
                sender_cell_groups=sender_cell_groups,
                receiver_cell_groups=receiver_cell_groups,
                plot_edges = plot_edges) +
    facet_wrap(~facet_group)

  return(gp)

}

#' this one subsets the contrast plot to a sub_cds dimensionality space
#' @param
#'
plot_sub_contrast_2 <- function(ccm,
                              cond_b_vs_a_tbl,
                              sub_cds,
                              cell_group = "cell_type",
                              scale_shifts_by=c("receiver", "sender", "none"),
                              edge_size=2,
                              cell_size=1,
                              q_value_thresh = 1.0,
                              group_label_size=2,
                              plot_labels = c("significant", "all", "none"),
                              fc_limits=c(-3,3),
                              sender_cell_groups=NULL,
                              receiver_cell_groups=NULL,
                              plot_edges = TRUE) {

  # colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
  colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[["cell_type"]]
  colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]


  # get cell groups in sub
  sub_cell_group = colData(sub_cds)[[ccm@ccs@info$cell_group]]

  # subset to select group
  ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[["cell_group"]] %in% sub_cell_group]
  ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]

  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% filter(cell_group %in% sub_cell_group)

  # switch coords to the new ones

  sub_umap = reducedDims(sub_cds)[["UMAP"]]
  reducedDims(ccm@ccs@cds)[["UMAP"]] = sub_umap[colnames(ccm@ccs@cds),]

  gp = plot_contrast(ccm,
                     cond_b_vs_a_tbl,
                     scale_shifts_by=scale_shifts_by,
                     edge_size=edge_size,
                     cell_size=cell_size,
                     q_value_thresh = q_value_thresh,
                     group_label_size=group_label_size,
                     plot_labels = plot_labels,
                     fc_limits=fc_limits,
                     sender_cell_groups=sender_cell_groups,
                     receiver_cell_groups=receiver_cell_groups,
                     plot_edges = plot_edges)
  return(gp)
}


subset_ccm <- function(ccm,
                       cond_b_vs_a_tbl,
                       sub_space = F,
                       facet_group = "major_group",
                       select_group = NULL) {

  ccm@ccs@cds = switch_umap_space(ccm@ccs@cds, sub_space = sub_space)

  colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
  colData(ccm@ccs@cds)[["facet_group"]] = colData(ccm@ccs@cds)[[facet_group]]

  cg_to_mg = as.data.frame(colData(ccm@ccs@cds)) %>%
    select("cell_group", "facet_group") %>%
    distinct()

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

  return(list(ccm = ccm, cond_b_vs_a_tbl = cond_b_vs_a_tbl))

}



#' projection wrap up
#' @param query_cds
#' @param ref_cds
#' @param directory_path a string giving the name of the directory from which to read the model files
#' @param transfer_type
run_projection <- function(query_cds,
                           ref_cds,
                           directory_path,
                           transfer_type = "cluster",
                           prefix = "") {

  query_cds <- load_transform_models(query_cds, directory_path)
  query_cds <- preprocess_transform(query_cds, method="PCA")
  query_cds <- align_beta_transform(query_cds)
  query_cds <- reduce_dimension_transform(query_cds, method="UMAP")

  append_umap_coords <- function(cds, prefix = "") {
    colData(cds)[[paste0(prefix,"umap3d_1")]] = reducedDims(cds)[["UMAP"]][,1]
    colData(cds)[[paste0(prefix,"umap3d_2")]] = reducedDims(cds)[["UMAP"]][,2]
    colData(cds)[[paste0(prefix,"umap3d_3")]] = reducedDims(cds)[["UMAP"]][,3]
    return(cds)
  }

  query_cds <- append_umap_coords(query_cds, prefix = prefix)

  ref_coldata = colData(ref_cds) %>% as.data.frame

  query_cds = transfer_cell_labels(query_cds,
                                   reduction_method = "UMAP",
                                   transform_models_dir = directory_path,
                                   ref_coldata = ref_coldata,
                                   ref_column_name = transfer_label)


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



comp_abund <- function(to, from, ccm) {
  to_abund = estimate_abundances(ccm, tibble("cell_state" = to))
  from_abund = estimate_abundances(ccm, tibble("cell_state" = from))
  return(compare_abundances(ccm, to_abund, from_abund))
}

plot_sub_path <- function(ccm,
                          cond_b_vs_a_tbl,
                          path_df,
                          sub_space = T,
                          cell_group = "cell_type_broad",
                          facet_group = "major_group",
                          select_group = NULL,
                          scale_shifts_by = "none",
                          plot_labels = "significant",
                          plot_edges = F,
                          q_value_thresh = 0.05,
                          edge_size = 1) {


  sub_ccm_tbl = subset_ccm(ccm,
                           cond_b_vs_a_tbl,
                           sub_space = sub_space,
                           facet_group = facet_group,
                           select_group = select_group)

  cell_groups = sub_ccm_tbl$cond_b_vs_a_tbl %>% pull(cell_group)

  # subset the path
  path_df = path_df %>%
    filter(to %in% cell_groups,
           from %in% cell_groups)

  path_df = add_umap_coords(path_df, centroids(sub_ccm_tbl$ccm@ccs))



  sub_contrast = plot_sub_contrast(ccm,
                                   cond_b_vs_a_tbl,
                                   scale_shifts_by = scale_shifts_by,
                                   select_group = select_group,
                                   plot_labels = plot_labels,
                                   plot_edges = plot_edges,
                                   q_value_thresh = q_value_thresh)


  gp =sub_contrast +
    geom_segment(data = path_df,
                 aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=umap_from_1,
                     yend = umap_from_2),
                 size=edge_size) +
    geom_segment(data = path_df,
                 aes(x = umap_from_1,
                     y = umap_from_2,
                     xend = (umap_to_1+umap_from_1)/2,
                     yend = (umap_to_2+umap_from_2)/2),
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm")),
                 size=edge_size)

  return(gp)

}



transitive.closure <- function(g,mat=FALSE,loops=TRUE){
  g <- as(g, "matrix")

  n <- ncol(g)

  matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
    while(i < pow) {
      z <- z %*% x
      y <- y+z
      i <- i+1
    }
    return(y)
  }

  h <- matExpIterativ(g,n)
  h <- (h>0)*1
  dimnames(h) <- dimnames(g)
  if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
  if (!mat) h <- as(h,"graphNEL")
  return(h)
}
transitive.reduction <- function(g){

  g = as(g, "matrix")
  # modified algorithm from Sedgewick book: just remove transitive edges instead of inserting them
  g = transitive.closure(g, mat=TRUE) # bug fix: algorithm only works for transitively closed graphs!
  g = g - diag(diag(g))
  type = (g > 1)*1 - (g < 0)*1
  for(y in 1:nrow(g)){
    for(x in 1:nrow(g)){
      if(g[x,y] != 0){
        for(j in 1:nrow(g)){
          if((g[y,j] != 0) & sign(type[x,j])*sign(type[x,y])*sign(type[y,j]) != -1){
            g[x,j] = 0
          }
        }
      }
    }
  }
  g
}


reduce_graph <- function(edges) {

  # transitive.reduction takes in an adjacency matrix
  G <- igraph::graph_from_data_frame(edges %>% dplyr::select(from, to), directed=TRUE )
  adj_matrix <- igraph::as_adjacency_matrix(G) %>% as.matrix()

  reduced_adj_matrix <- transitive.reduction(adj_matrix)
  reduced_graph <- igraph::graph_from_adjacency_matrix(reduced_adj_matrix)
  reduced_graph

}


# convert cluster to cell group

convert_cluster_to_cell_group <- function(ccs, cell_group) {

  colData(ccs@cds)[["cell_group"]] = colData(ccs@cds)[[cell_group]]

  cluster_to_cellgroup = colData(ccs@cds) %>%
    as.data.frame() %>%
    group_by(cluster, cell_group) %>%
    tally() %>%
    top_n(1)

  return(cluster_to_cellgroup)

}


#'
#' @param size_factor_group
#'
#'
compute_sizefactors <- function(ccs, size_factor_group) {


  colData(ccs)[["size_factor_group"]]= colData(ccs)[[size_factor_group]]

  size_factors_list = lapply(unique(colData(ccs)$size_factor_group), function(sfg) {
    ccs_sfg = ccs[,colData(ccs)$size_factor_group == sfg]
    ccs_sfg = estimate_size_factors(ccs_sfg)
    data.frame(size_factors(ccs_sfg))
  })

  size_factors_all = do.call(rbind, size_factors_list)

  sfs = size_factors_all[names(colData(ccs)$Size_Factor),]
  names(sfs) = names(colData(ccs)$Size_Factor)

  # sfs = size_factors_all[,1]
  # names(sfs) = rownames(size_factors_all)
  return(sfs)
}




