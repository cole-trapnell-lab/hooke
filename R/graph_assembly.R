#' Orient graph edges from a PLNnetwork using a contrast between conditions
#'
#'@param ccm A cell_count_model
#'@param cond_b_vs_a_tbl A contrast between two conditions as returned by compare_abundances()
collect_pln_graph_edges <- function(ccm,
                                    cond_b_vs_a_tbl,
                                    log_abundance_thresh = 1-5){
  pln_model = model(ccm)

  abundance_corr_tbl = correlate_abundance_changes(pln_model, cond_b_vs_a_tbl)

  abundance_corr_tbl = abundance_corr_tbl %>% dplyr::filter(
    (to_log_abund_x > log_abundance_thresh | to_log_abund_y > log_abundance_thresh) &  # Keep if the "to" node is above abundance thresh in at least one condition
      (from_log_abund_x > log_abundance_thresh | from_log_abund_y > log_abundance_thresh) # Keep if the "to" node is above abundance thresh in at least one condition
  )

  corr_edge_delta_abund = abundance_corr_tbl

  corr_edge_delta_abund = corr_edge_delta_abund %>% dplyr::mutate(scaled_weight = -pcor)

  corr_edge_delta_abund = corr_edge_delta_abund %>%
    dplyr::mutate(edge_type = dplyr::case_when(
      from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
      from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor < 0 ~ "directed_to_from",
      from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor < 0 ~ "directed_from_to",
      from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor < 0 ~ "undirected",
      TRUE ~ "hidden"
    ))

  # Fix the edges so all directed edges are directed_from_to:
  backwards_edges = corr_edge_delta_abund %>% dplyr::filter(edge_type == "directed_to_from")
  ok_edges = corr_edge_delta_abund %>% dplyr::filter(edge_type != "directed_to_from")

  be_colnames = colnames(backwards_edges)
  be_colnames = unlist(lapply(be_colnames, stringr::str_replace, "to", "from"))
  from_cols = grepl("from", colnames(backwards_edges))
  be_colnames[from_cols] = unlist(lapply(be_colnames[from_cols], stringr::str_replace, "from", "to"))
  colnames(backwards_edges) = be_colnames
  if (nrow(backwards_edges) != 0) {
    backwards_edges$edge_type = "directed_from_to"
  }
  backwards_edges = backwards_edges[colnames(ok_edges)]
  corr_edge_coords_umap_delta_abund = rbind(ok_edges, backwards_edges)
  return (corr_edge_coords_umap_delta_abund)
}

return_igraph <- function(model, type = "partial_cor", remove.isolated=FALSE,
                          edge.color = c("#F8766D", "#00BFC4"),
                          output = "igraph"){

  net <- model$latent_network(type = type)

  if (output == "igraph") {

    G <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)
    igraph::V(G)$label <- colnames(net)

    ## Nice nodes
    V.deg <- igraph::degree(G)/sum(igraph::degree(G))
    igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
    igraph::V(G)$size <- V.deg * 100
    igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
    igraph::V(G)$frame.color <- NA

    ## Nice edges
    igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
    if (type == "support")
      igraph::E(G)$width <- abs(igraph::E(G)$weight)
    else
      igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)


    if (remove.isolated) {
      G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
    }
    G$layout <- igraph::layout_in_circle
  }
  if (output == "corrplot") {

    if (ncol(net) > 100)
      colnames(net) <- rownames(net) <- rep(" ", ncol(net))
    G <- net
    diag(net) <- 0
    corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
  }
  G
}


#' Extract a partitioned abstract graph from a Monocle cell_data_set object
#'
#' QUESTION: Do the cells need to be grouped by monocle cluster? Or can they
#' be grouped arbitrarily (e.g. by cell type? Would be good to generalize if
#' possible. For now, only works when cells are grouped by Monocle cluster
#'
#' @param cds A cell_data_set object. cluster_cells() must have been called.
#' @param reduction_method The coordinate space in which to build the graph
get_paga_graph <- function(cds, reduction_method = "UMAP") {

  cluster_result <- cds@clusters[[reduction_method]]$cluster_result
  #clusters <- ccm@ccs@cds@clusters[[reduction_method]]$clusters
  #partitions <- ccm@ccs@cds@clusters[[reduction_method]]$partitions

  #colData(cds)$cluster = clusters
  #colData(cds)$partition = partitions

  #coldata = colData(ccm@ccs@cds) %>% as_tibble() %>% dplyr::rename(embryo = Oligo) %>%
  #  mutate(timepoint=as.numeric(stringr::str_remove(timepoint, "hpf")))
  #coldata = ccm@ccs@metadata[["cell_group_assignments"]]

  cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g,
                                                     cluster_result$optim_res,
                                                     0.05, FALSE)


  #cluster_time = coldata %>% group_by(cluster) %>% summarise(avg_time = mean(as.numeric(timepoint)))
  #cluster_info = merge(cluster_maj, cluster_time, by="cluster") %>%
  #  mutate("label" = paste0(cell_type, ".", round(avg_time), ".", cluster)) %>%
  #  arrange(as.numeric(cluster))

  cluster_g <- cluster_graph_res$cluster_g
  cluster_g <- igraph::set_vertex_attr(cluster_g, "name", value = stringr::str_replace(igraph::V(cluster_g)$name, "cell_membership", ""))
  cluster_g
}





