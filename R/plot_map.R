
# first plot map with cytoscape coords
# G <- readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/MLP/lineage_map_igraph_210824.rds")
# load a layout
# cytoscape_pos <- readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/in progress/cytoscape/cytoscape_positions_20210824.rds")



get_time_levels <- function(ccs,
                            G,
                            num_levels = 10) {


  colData(ccs@cds)$cell_group = colData(ccs@cds)[[ccs@info$cell_group]]

  timepoint_df = colData(ccs@cds) %>%
    as.data.frame %>%
    group_by(cell_group) %>%
    summarize(avg_time = mean(as.numeric(timepoint), na.rm=T))

  # separate the unattached and put them on the bottom
  isolated_nodes = tidygraph::as_tbl_graph(G) %>%
    tidygraph::activate(nodes) %>%
    mutate(isolated = tidygraph::node_is_isolated()) %>%
    as.data.frame()

  timepoint_df = timepoint_df %>%
    left_join(isolated_nodes, c("cell_group" = "name")) %>%
    mutate(min_time = ifelse(isolated, 100, avg_time))

  level_df = timepoint_df %>%
    arrange(avg_time) %>%
    mutate(rn = row_number()) %>%
    mutate(group = cut(rn, num_levels, labels=F)) %>%
    mutate(group = ifelse(avg_time == 100, group + 1, group)) %>%
    tibble::column_to_rownames("cell_group")

  return(level_df)
}

get_layout <- function(G, level_df) {

  lay1 <- igraph::layout_with_sugiyama(G, layers=level_df[igraph::V(G)$name,][["group"]])

  return(lay1)
}


# this only works if you have cell type broad

# ccs = ref_ccs
# ccm = ref_ccm
# cond_b_vs_a_tbl = cond_18_vs_30_tbl
# path_df = green_edges %>% select(from, to) %>% distinct()
# #
# # subset_col = "major_group"
# # subset_by = "mesoderm"
#
#
# tmp_function(ref_ccm,
#              cond_18_vs_30_tbl,
#              green_edges)


get_sub_layout <- function(ccm,
                           path_df,
                           subset_col,
                           subset_by,
                           arrow.gap = 0.02) {

  # subset to select group
  if ( (is.null(subset_col) & is.null(subset_by)) == FALSE) {
    ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[[subset_col]] == subset_by]
    ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]

  }

  cell_groups = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]] %>% unique()
  path_df = path_df %>% select(from, to) %>%filter(from %in% cell_groups, to %in% cell_groups)

  # are there any in path_df that don't make it into
  lone_nodes = setdiff(cell_groups, unique(union(path_df$from, path_df$to)))
  # alone_df = data.frame(from = alone, to = alone)
  # path_df = rbind(path_df, alone_df) %>% distinct()

  G = igraph::graph_from_data_frame(path_df, directed = T)
  for (n in lone_nodes) {
    G <- igraph::add_vertices(G, 1, name = n)
  }

  level_df = get_time_levels(ccm@ccs, G)
  lay1 = get_layout(G, level_df)

  g = ggnetwork::ggnetwork(G, layout = lay1$layout, arrow.gap = arrow.gap, multiple = T, loops=T)

  return(list(lay1 = lay1, g = g))
}



plot_edges <- function(ccm,
                       cond_b_vs_a_tbl,
                       path_df,
                       subset_col = NULL,
                       subset_by = NULL,
                       arrow.gap = 0.02,
                       q_value_thresold = 0.05,
                       legend.position = "right",
                       cluster_to_cellgroup = NULL) {

  # subset to select group
  if ( (is.null(subset_col) & is.null(subset_by)) == FALSE) {
    ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[[subset_col]] == subset_by]
    ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]

  }

  cell_groups = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]] %>% unique()
  path_df = path_df %>% filter(from %in% cell_groups, to %in% cell_groups)

  G = igraph::graph_from_data_frame(path_df, directed = T)

  # then can plot by cytoscape
  if (ccm@ccs@info$cell_group == "cell_type_broad") {

    cytoscape_pos = cytoscape_pos %>% tibble::column_to_rownames("id")
    geo = cytoscape_pos[igraph::V(G)$name,] %>% select(x,y) %>% as.matrix
    g = ggnetwork::ggnetwork(G, layout = geo, arrow.gap = arrow.gap)

  }

  level_df = get_time_levels(ccm@ccs, G)
  lay1 = get_layout(G, level_df)
  g = ggnetwork::ggnetwork(G, layout = lay1$layout, arrow.gap = arrow.gap, multiple = T, loops=T)


  if (is.null(cond_b_vs_a_tbl) == FALSE) {
    g = g %>%
      left_join(cond_b_vs_a_tbl, by = c("name" = "cell_group")) %>%
      mutate(delta_log_abund = ifelse(delta_q_value < q_value_thresold, delta_log_abund, 0)) %>%
      mutate(delta_log_abund = tidyr::replace_na(delta_log_abund, 0))

  }

  # g replace names with cluster
  # cluster_to_cellgroup

  if (!is.null(cluster_to_cellgroup)) {
    g = g %>%
      left_join(cluster_to_cellgroup, by = c("name" = "cluster")) %>%
      select(-name) %>% rename(name = "cell_group")
  }


  # if it has edge weights, scale by weight
  if ("weight" %in% colnames(path_df)) {

    g = g %>% mutate(scaled_weight = weight / max(abs(weight), na.rm = T))

    p = ggplot(g, aes(x, y, xend = xend, yend = yend)) +
      ggnetwork::geom_edges(aes(size = 0.35*scaled_weight),
                                arrow = arrow(length = unit(6, "pt"), type="closed")) +
      scale_size_identity()

  } else {

    p = ggplot(g, aes(x, y, xend = xend, yend = yend)) +
      ggnetwork::geom_edges(arrow = arrow(length = unit(6, "pt"), type="closed"))

  }

  p = p + ggnetwork::geom_nodes(aes(fill = delta_log_abund),
                              size = 7,colour="black",shape=21) +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "red4") +
        ggnetwork::geom_nodetext_repel(aes(label = name), size = 3) +
        monocle3:::monocle_theme_opts() +
        ggnetwork::theme_blank() +
        scale_size_identity() +
        theme(legend.position = legend.position)

  return(p)

}


# this organizes things into 4 quarters
#
plot_everything <- function(ccm = ccm,
                            path = path,
                            cond_b_vs_a_tbl = cond_b_vs_a_tbl,
                            q_value_thresold = 0.05,
                            arrow.gap = 0.02,
                            legend.position = "right") {

  lay_mesoderm = get_sub_layout(ccm,
                                path_df = path,
                                subset_col = "major_group",
                                subset_by = "mesoderm")

  lay_cns = get_sub_layout(ccm,
                           path_df = path,
                           subset_col = "major_group",
                           subset_by = "CNS")

  lay_peri = get_sub_layout(ccm,
                            path_df = path,
                            subset_col = "major_group",
                            subset_by = "periderm-other")

  lay_mesen = get_sub_layout(ccm,
                             path_df = path,
                             subset_col = "major_group",
                             subset_by = "mesenchyme-fin")



  q1 = lay_mesoderm$g %>% select(x,y, name) %>% distinct()
  q2 = lay_cns$g %>% select(x,y, name) %>% distinct() %>% mutate(x=x-1)
  q3 = lay_peri$g %>% select(x,y, name) %>% distinct() %>% mutate(x=x-2)
  q4 = lay_mesen$g %>% select(x,y, name) %>% distinct() %>% mutate(x=x+1)


  coords = rbind(q1, q2, q3, q4)
  rownames(coords) = NULL
  coords = tibble::column_to_rownames(coords, "name")

  G = igraph::graph_from_data_frame(path, directed = T)

  # what names arent in coords

  geo = coords[igraph::V(G)$name,] %>% select(x,y) %>% as.matrix
  g = ggnetwork::ggnetwork(G, layout = geo, arrow.gap = arrow.gap, multiple = T, loops=T)

  g = g %>%
    left_join(cond_b_vs_a_tbl, by = c("name" = "cell_group")) %>%
    mutate(delta_log_abund = ifelse(delta_q_value < q_value_thresold, delta_log_abund, 0)) %>%
    mutate(delta_log_abund = tidyr::replace_na(delta_log_abund, 0))

  if ("weight" %in% colnames(g)) {

    g = g %>% mutate(scaled_weight = weight / max(abs(weight), na.rm = T))

    p = ggplot(g, aes(x, y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(aes(size = 0.35*scaled_weight),
                              arrow = arrow(length = unit(6, "pt"), type="closed")) +
        scale_size_identity()

  } else {

    p = ggplot(g, aes(x, y, xend = xend, yend = yend)) +
      ggnetwork::geom_edges(arrow = arrow(length = unit(6, "pt"), type="closed"))

  }

  p = p +
      ggnetwork::geom_nodes(aes(fill = delta_log_abund),
                            size = 7,colour="black",shape=21) +
      scale_fill_gradient2(low = "darkblue", mid = "white", high = "red4") +
      ggnetwork::geom_nodetext_repel(aes(label = name), size = 3) +
      monocle3:::monocle_theme_opts() +
      ggnetwork::theme_blank()+
      theme(legend.position = legend.position)

  return(p)

}


compare_to_cytoscape = function(ccm,
                                path_df,
                                legend.position = "right",
                                arrow.gap = 0.02) {

  # if not already cell type broad,
  # collapse edges to match that
  if (ccm@ccs@info$cell_group != "cell_type_broad") {
    cts_to_ctb = colData(ccm@ccs@cds) %>%
      as.data.frame() %>%
      group_by(cell_type_sub, cell_type_broad) %>%
      tally()

    path_df = path_df %>%
      left_join(cts_to_ctb, by = c("from" = "cell_type_sub")) %>%
      left_join(cts_to_ctb, by = c("to" = "cell_type_sub"), suffix = c(".from", ".to")) %>%
      select(-c(from, to)) %>%
      select("from" = cell_type_broad.from, "to" = cell_type_broad.to) %>% distinct()

  }


  cytoscape_pos <- readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/in progress/cytoscape/cytoscape_positions_20210824.rds")
  G_lm <- readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/MLP/lineage_map_igraph_210824.rds")
  # get the coords of the
  cytoscape_pos = cytoscape_pos %>% tibble::column_to_rownames("id")
  geo = cytoscape_pos[igraph::V(G_lm)$name,] %>% select(x,y) %>% as.matrix


  # path df add null nodes to path_df
  lone_nodes = setdiff(igraph::V(G_lm)$name, unique(union(path_df$from, path_df$to)))
  G = igraph::graph_from_data_frame(path_df, directed = T)
  for (n in lone_nodes) {
    G <- igraph::add_vertices(G, 1, name = n)
  }

  g1 = ggnetwork::ggnetwork(G_lm, layout = geo, arrow.gap = arrow.gap)
  g2 = ggnetwork::ggnetwork(G, layout = geo, arrow.gap = arrow.gap)

  # what edges are shared
  # should make those green
  # missing ones are gray
  # wrong are red?

  p = ggplot(mapping=aes(x, y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(data = g1, arrow = arrow(length = unit(6, "pt"), type="closed"), color="red") +
    ggnetwork::geom_nodes(data = g1, size = 7,colour="black",shape=21) +
    ggnetwork::geom_nodetext_repel(data = g1, aes(label = name), size = 3) +

    ggnetwork::geom_edges(data = g2, arrow = arrow(length = unit(6, "pt"), type="closed")) +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank()+
    theme(legend.position = legend.position)

  return(p)

}







