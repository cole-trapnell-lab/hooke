#' Plot a UMAP colored by how cells shift in a given contrast
#'
#' @param ccm A cell_count_model object.
#' @param criterion a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".
#' @param scale_shifts_by
#' @param q_value_thresh
#' @param plot_labels
#' @return an updated cell_count_model object
#' @export
#' @import ggplot2
#' @import dplyr
plot_contrast <- function(ccm,
                          #umap_centers,
                          cond_b_vs_a_tbl,
                          log_abundance_thresh = -5,
                          scale_shifts_by=c("receiver", "sender", "none"),
                          #cell_group="cluster",
                          edge_size=2,
                          cell_size=1,
                          q_value_thresh = 1.0,
                          group_label_size=2,
                          plot_labels = c("significant", "all", "none"),
                          fc_limits=c(-3,3),
                          sender_cell_groups=NULL,
                          receiver_cell_groups=NULL,
                          plot_edges = TRUE){

  umap_centers = centroids(ccm@ccs)

  umap_centers_delta_abund = umap_centers
  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))
  umap_centers_delta_abund = dplyr::left_join(umap_centers_delta_abund, cond_b_vs_a_tbl, by=c("cell_group"="cell_group"))
  umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  umap_centers_delta_abund = umap_centers_delta_abund %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))

  #cond_b_vs_a_tbl$delta_q_value = p.adjust(cond_b_vs_a_tbl$delta_q_value, method = "BH")
  # umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

  corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                              umap_centers_delta_abund,
                                                              log_abundance_thresh)
  directed_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("directed_to_from", "directed_from_to"))
  undirected_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("undirected"))

  if (is.null(sender_cell_groups) == FALSE){
    directed_edge_df = directed_edge_df %>% dplyr::filter(from %in% sender_cell_groups)
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(from %in% sender_cell_groups | to %in% sender_cell_groups)
  }

  if (is.null(receiver_cell_groups) == FALSE){
    directed_edge_df = directed_edge_df %>% dplyr::filter(to %in% receiver_cell_groups)
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(from %in% receiver_cell_groups | to %in% receiver_cell_groups)
  }

  # corr_edge_coords_umap_delta_abund = left_join(corr_edge_coords_umap_delta_abund,
  #                                               umap_centers,
  #                                               by=c("from"="cell_group"))

  if (scale_shifts_by == "sender"){
    directed_edge_df = directed_edge_df %>%
      dplyr::group_by(to) %>%
      dplyr::mutate(flow_factor = -pmin(0, pcor),
             total_weight = sum(flow_factor),
             scaled_weight  = flow_factor / total_weight)
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }else if (scale_shifts_by == "receiver"){
    directed_edge_df = directed_edge_df %>%
      dplyr::group_by(from) %>%
      dplyr::mutate(flow_factor = -pmin(0, pcor),
             total_weight = sum(flow_factor),
             scaled_weight  = flow_factor / total_weight)
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }else{
    directed_edge_df = directed_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
    undirected_edge_df = undirected_edge_df %>%
      dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  }

  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)

  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]

  plot_df = dplyr::left_join(plot_df,
                             cond_b_vs_a_tbl,
                              # cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                              by=c("cell_group"="cell_group"))


  if (is.null(fc_limits)) {
    fc_limits = range(plot_df$delta_log_abund)
  } else {
    min = fc_limits[1]
    max = fc_limits[2]
    plot_df = plot_df %>%
      mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
      mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  }

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
    geom_point(
      data = plot_df %>%
        arrange(!is.na(abs(delta_log_abund)),
                abs(delta_log_abund)),
      aes(umap2D_1, umap2D_2, color = delta_log_abund),
      size = cell_size,
      stroke = 0
    ) +
    scale_color_gradient2(
      low = "#122985",
      mid = "white",
      high = "red4",
      na.value = "white",
      limits = fc_limits
    )  +
    #theme_void() +
    #theme(legend.position = "none") +
    monocle3:::monocle_theme_opts() +
    xlab("UMAP 1") + ylab("UMAP 2")


  if (plot_edges) {
    gp = gp  +
      geom_segment(data = undirected_edge_df,
                   aes(x = to_umap_1,
                       y = to_umap_2,
                       xend=from_umap_1,
                       yend = from_umap_2,
                       size=edge_size * scaled_weight),
                   #size=edge_size / 4,
                   color="lightgray") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = to_umap_1,
                       y = to_umap_2,
                       xend=from_umap_1,
                       yend = from_umap_2,
                       size=edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = to_umap_1,
                       y = to_umap_2,
                       xend=(to_umap_1+from_umap_1)/2,
                       yend = (to_umap_2+from_umap_2)/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = from_umap_1,
                       y = from_umap_2,
                       xend=to_umap_1,
                       yend = to_umap_2,
                       size=edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = from_umap_1,
                       y = from_umap_2,
                       xend=(from_umap_1+to_umap_1)/2,
                       yend = (from_umap_2+to_umap_2)/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_size_identity()
  }

  if (plot_labels != "none") {
    label_df = umap_centers_delta_abund
    if (plot_labels == "significant")
      label_df = label_df %>% filter(delta_log_abund != 0)
    gp <- gp + ggrepel::geom_label_repel(data = label_df,
                                         mapping = aes(umap_1, umap_2, label=cell_group),
                                         size=I(group_label_size),
                                         fill = "white")
  }
  return(gp)
}

#' returns different color palettes
#' @param num_colors the number of colors needed
#' @param
get_colors <- function(num_colors, type = "rainbow") {


  if (type == "rainbow") {

    colors =
      c("18h" = "#DF4828",
        "24h" = "#E78C35",
        "36h" = "#F6C141",
        "48h" = "#4EB265",
        "72h" = "#1965B0")

  } else if (type == "vibrant")  {

    colors =
      c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
        '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')

  } else {

    colors =
        c('#4477AA',
          '#EE6677',
          '#228833',
          '#CCBB44',
          '#66CCEE',
          '#AA3377',
          '#BBBBBB')

  }

  full_spectrum =  colorRampPalette(colors)(num_colors)

  return(full_spectrum)


}


#' rename this later?
#' probably needs fixing
my_plot_cells <- function(data,
                      color_cells_by = "cluster",
                      cell_size=1,
                      legend_position="none",
                      residuals = NULL,
                      effect_sizes = NULL,
                      cond_b_vs_a_tbl = NULL,
                      cell_group = NULL,
                      q_value_thresh = 1.0) {

  if (class(data) == "cell_count_set") {
    cds = data@cds
    plot_df = as.data.frame(colData(cds))
  } else if (class(data) == "cell_data_set") {
    cds = data
    plot_df = as.data.frame(colData(cds))
    # plot_df = data@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else if (class(data) == "cell_count_model") {
    cds = data@ccs@cds
    plot_df = as.data.frame(colData(cds))
    plot_df$cell_group = data@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else {
    print("some error message")
  }


  plot_df$cell = row.names(plot_df)
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[plot_df$cell,2]


  # could be an vector that corresponds to a label
  if (!is.null(residuals)) {
    res_df = as.data.frame(residuals) %>% rownames_to_column("cell_group")
    colnames(res_df)[1] = "residuals"
    plot_df = plot_df %>% left_join(res_df, by="cell_group")
    color_cells_by = "residuals"

  }

  else if (!is.null(effect_sizes)) {

    plot_df = plot_df %>% left_join(effect_sizes, by="cell_group")
    color_cells_by = "delta_log_needed"

  }

  else if (!is.null(cond_b_vs_a_tbl)) {

    cond_b_vs_a_tbl = cond_b_vs_a_tbl %>%
      dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

    plot_df = dplyr::left_join(plot_df,
                               cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                               by=c("cell_group"="cell_group"))

    color_cells_by = "delta_log_abund"

  }

  plot_df$color_cells_by = plot_df[[color_cells_by]]

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
    theme_void() +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()


  if (color_cells_by == "timepoint") {

    num_colors = unique(plot_df$timepoint) %>% sort() %>% length()
    full_spectrum_timepoint = get_time_colors(num_colors)
    names(full_spectrum_timepoint) = unique(plot_df$timepoint) %>% sort()

    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = as.character(timepoint)),
        size = cell_size,
        stroke = 0
      ) +
      scale_color_manual(values = full_spectrum_timepoint)

  } else if (color_cells_by %in% c("residuals", "delta_log_needed")) {
    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = color_cells_by),
        size = cell_size,
        stroke = 0
      ) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = c(-3,3)
      )

  } else if (color_cells_by == "delta_log_abund") {
    gp = gp  +
      geom_point(
        data = plot_df %>%
          arrange(!is.na(abs(delta_log_abund)),
                  abs(delta_log_abund)),
        aes(umap2D_1, umap2D_2, color = delta_log_abund),
        size = cell_size,
        stroke = 0
      ) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits=c(-3,3)
      )
  } else if (color_cells_by == "none") {
    # do nothing
  }
  else {
    num_colors = unique(plot_df[[color_cells_by]]) %>% sort() %>% length()
    full_spectrum = get_colors(num_colors, "vibrant")
    names(full_spectrum) = unique(plot_df[[color_cells_by]]) %>% sort()

    gp = gp + geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2,
          color = color_cells_by),
      size = cell_size,
      stroke = 0
    ) +
      scale_color_manual(values = full_spectrum)
  }

  return(gp)

}


#' plots a path on top of
#' @param data
#' @param path_df
#' @param edge_size
#' @param path_color
#' @param color_cells_by
#' @param residuals
#' @param cond_b_vs_a_tbl
plot_path <- function(data,
                      path_df = path_df,
                      edge_size=2,
                      path_color = "black",
                      color_cells_by = "cluster",
                      color_path_by = NULL,
                      residuals = NULL,
                      cond_b_vs_a_tbl = NULL) {


  gp = my_plot_cells(data, color_cells_by = color_cells_by, residuals = residuals, cond_b_vs_a_tbl = cond_b_vs_a_tbl)

  if (class(data) == "cell_count_set") {
    umap_centers = centroids(data)
  } else if (class(data) == "cell_count_model") {
    umap_centers = centroids(data@ccs)
  }

  path_df = add_umap_coords(path_df, umap_centers)

  if (is.null(color_path_by) == FALSE) {

    gp = gp +
      ggnewscale::new_scale_color() +
      geom_segment(data = path_df,
                           aes(x = umap_to_1,
                               y = umap_to_2,
                               xend=umap_from_1,
                               yend = umap_from_2,
                               color = get(color_path_by)),
                           size=2 ) +
      geom_segment(data = path_df,
                   aes(x = umap_from_1,
                       y = umap_from_2,
                       xend = (umap_to_1+umap_from_1)/2,
                       yend = (umap_to_2+umap_from_2)/2,
                       color = get(color_path_by)),
                   size=2,
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white"
      )

  } else {
    gp = gp + geom_segment(data = path_df,
                           aes(x = umap_to_1,
                               y = umap_to_2,
                               xend=umap_from_1,
                               yend = umap_from_2),
                           color = path_color,
                           size = 2 ) +
      geom_segment(data = path_df,
                   aes(x = umap_from_1,
                       y = umap_from_2,
                       xend = (umap_to_1+umap_from_1)/2,
                       yend = (umap_to_2+umap_from_2)/2),
                   size=2,
                   color=path_color,
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))
  }

  return(gp)

}


add_edge <- function(gp,
                     path_df,
                     umap_centers,
                     size = 1,
                     color = "black") {
    path_df = add_umap_coords(path_df, umap_centers)

    gp = gp +
      geom_segment(data = path_df,
              aes(x = umap_to_1,
                  y = umap_to_2,
                  xend = umap_from_1,
                  yend = umap_from_2),
              size = size,
              color = color) +
      geom_segment(data = path_df,
              aes(x = umap_from_1,
                  y = umap_from_2,
                  xend = (umap_to_1+umap_from_1)/2,
                  yend = (umap_to_2+umap_from_2)/2),
              size = size,
              color = color,
              linejoin='mitre',
              arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))

    return(gp)
}


#' this currently just plots the mappy thing with the abundance results
#' to be completed later to handle automatic plotting
#' @param G igraph object of mappy thing
#' @param cytoscape_pos coordinates from cytoscape
#' @param cond_b_vs_a_tbl output of compare_abundances
#'
plot_graph <- function(G,
                       cytoscape_pos,
                       cond_b_vs_a_tbl,
                       arrow.gap = 0.02,
                       scale = F,
                       q_value_thresold = 0.05,
                       color_nodes_by = "delta_log_abund",
                       legend.position = "right") {


  cytoscape_pos = cytoscape_pos %>% column_to_rownames("id")
  geo = cytoscape_pos[igraph::V(G)$name,] %>% select(x,y) %>% as.matrix
  g = ggnetwork::ggnetwork(G, layout = geo, arrow.gap=arrow.gap) %>%
    left_join(cond_b_vs_a_tbl, by = c("name" = "cell_group")) %>%
    mutate(delta_log_abund = ifelse(delta_q_value < q_value_thresold, delta_log_abund, 0)) %>%
    mutate(delta_log_abund = replace_na(delta_log_abund, 0))

  p = ggplot(g, aes(x, y, xend = xend, yend = yend)) +
      ggnetwork::geom_edges(arrow = arrow(length = unit(6, "pt"), type="closed")) +
      ggnetwork::geom_nodes(aes(fill = delta_log_abund),
                            size = 7,colour="black",shape=21) +
      scale_fill_gradient2(low = "darkblue", mid = "white", high = "red4") +
      ggnetwork::geom_nodetext_repel(aes(label = name), size = 3) +
      monocle3:::monocle_theme_opts()+
      theme_blank()+
      theme(legend.position = legend.position)
  p
  return(p)

}



# plot igraph
#' plots the resulting path as an igraph
#' instead of as a UMAP
#' @param data
#' @param edges
#' @param color_nodes_by
#' @param arrow.gap
#' @param scale
plot_map <- function(data, edges, color_nodes_by = "", arrow.gap = 0.02, scale = F) {

  if (class(data) == "cell_count_set") {
    ccs = data
    plot_df = as.data.frame(colData(cds))
  } else if (class(data) == "cell_count_model") {
    ccs = data@ccs
  } else {
    print("some error message")
  }

  umap_centers = centroids(ccs)
  cell_groups = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  nodes = data.frame(id = cell_groups)
  if ("source" %in% colnames(edges)) {
    edges = edges %>% rename("from" = source, "to" = target)
  }
  edges = edges %>% select(from,to)

  no_edge = setdiff(nodes$id, union(edges$from, edges$to))
  edges = rbind(edges, data.frame("from" = no_edge, "to" = no_edge))
  n = network(edges %>% dplyr::select(from,to), directed = T, loops = T)
  nodes = nodes[match(network.vertex.names(n), nodes$id),]
  n %v% "id" = network.vertex.names(n)

  merged_coords = ggnetwork(n) %>% select(id, vertex.names) %>% unique() %>%
    left_join(umap_centers, by = c("id"="cell_group"))
  rownames(merged_coords) = merged_coords$id
  coords = merged_coords[network.vertex.names(n),] %>% select(umap_1,umap_2)
  geo = as.matrix(sapply(coords, as.numeric))

  g <- ggnetwork(x = n, layout = geo, arrow.gap=arrow.gap, scale = scale)

  show(ggplot(g, aes(x, y, xend = xend, yend = yend)) +
          geom_edges(arrow = arrow(length = unit(6, "pt"), type="closed")) +
          geom_nodes(size = 7,colour="black",shape=21) +
          geom_nodetext_repel(aes(label = id), size=3) +
          theme_blank())

}


#' plots a grn
#' @param g a ggnetwork dataframe object
#' @color_nodes_by
plot_grn <- function(g,
                     color_nodes_by,
                     regulator_score_min = NULL,
                     regulator_score_max = NULL,
                     legend.position = "right") {


  # colors = c( "#735CDD", "#1B9AAA","#EF476F", "#101D42")
  # names(colors) = c("mesodermal progenitor cells (contains PSM)",
  #                   "myoblast/fast-committed myocyte",
  #                   "fast-committed myocyte",
  #                   "mature fast muscle")

  num_colors = unique(g[[color_nodes_by]]) %>% length()
  full_spectrum = get_colors(num_colors, "vibrant")
  names(full_spectrum_timepoint) = unique(plot_df$timepoint) %>% sort()

  if (is.null(regulator_score_min)) {
    g = g %>% mutate(regulator_score = ifelse(regulator_score < regulator_score_min, regulator_score_min, regulator_score))

  }
  if (is.null(regulator_score_max)) {
    g = g %>% mutate(regulator_score = ifelse(regulator_score > regulator_score_max, regulator_score_max, regulator_score))

  }

  p = ggplot(mapping = aes(x, y, xend = xend, yend = yend)) +
    # draw activator edges
    geom_nodes(data = g, size=0.5) +
    geom_edges(data = g %>% filter(edge_type == "activator"),
               aes(size = 0.35*scaled_weight),
               arrow = arrow(length = unit(1.5, "pt"), type="closed"),
               curvature = 0.1) +
    # draw repressor edges
    geom_edges(data = g %>% filter(edge_type == "repressor"),
               aes(size = 0.35*scaled_weight),
               arrow = arrow(angle = 90, length = unit(.075, "cm")),
               curvature = 0) +
    geom_nodelabel_repel(data = g ,
                         aes(fill = cell_type_broad,
                             label = gene_short_name,
                             size = 1.25*abs(regulator_score)),
                         color="white",
                         label.size = 0.2,
                         box.padding = unit(0.1, "lines"),
                         label.padding = 0.15,
                         force = 0.1) +
    scale_fill_manual(values = colors) +
    scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    theme_blank() +
    theme(legend.position = legend.position)


}

# make sure ccs is made on clusters

plot_paga <- function(ccs,
                      color_cells_by,
                      legend_position = "right") {


  assertthat::assert_that(
    tryCatch(expr = ccs@info$cell_group == "cluster",
             error = function(e) FALSE),
    msg = "cell group must be cluster")

  umap_centers = centroids(ccs) %>% tibble::column_to_rownames("cell_group")

  paga_graph = get_paga_graph(ccs@cds)

  colData(ccs@cds)[["color_cells_by"]] = colData(ccs@cds)[[color_cells_by]]

  cluster_to_cellgroup = colData(ccs@cds) %>%
                          as.data.frame() %>%
                          group_by(cluster, color_cells_by) %>%
                          tally() %>%
                          top_n(1) %>%
                          tibble::column_to_rownames("cluster")

  igraph::V(paga_graph)$cell_group = cluster_to_cellgroup[igraph::V(paga_graph)$name, "color_cells_by"]

  coords = umap_centers[igraph::V(paga_graph),1:2] %>% as.matrix()

  p = ggnetwork::ggnetwork(paga_graph, layout=coords, arrow.gap= 0.02) %>%
      ggplot(aes(x, y, xend = xend, yend = yend)) +
      ggnetwork::geom_edges() +
      ggnetwork::geom_nodes(aes(fill = cell_group), size = 7, color="black",shape=21) +
      ggnetwork::theme_blank() +
      theme(legend.position = legend_position)

  return(p)

}




#' helper function for adding attributes to network node
#' assumes gene id is in the dataframe
# add_node_attribute <- function(n, df, colname) {
#   df = df[match(network.vertex.names(n), df$gene_id),]
#   n %v% colname = df[[colname]]
#   return(n)
# }
#
#
# plot_regulators <- function(edges) {
#
#   n = network(edges, directed = F, loops = F)
#
#   # add regulator score to the nodes
#
#   add_node_attribute <- function(n, df, colname) {
#     df = df[match(network.vertex.names(n), df$gene_id),]
#     n %v% colname = reg_score_df[[colname]]
#     return(n)
#   }
#
#   n = add_node_attribute(n, reg_score_df, "regulator_score")
#
#   g = ggnetwork(n, arrow.gap=0.02)
#
#   ggplot(mapping=aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
#     geom_edges(data = g %>% filter(edge_type == "undirected"), color="gray") +
#
#     # geom_edges(data = g %>% filter(edge_type != "undirected"),
#     #            arrow = arrow(length = unit(3, "pt"), type="closed")) +
#
#     geom_edges(data = g %>% filter(edge_type != "undirected"),
#                arrow = arrow(angle = 90, length = unit(.075, "cm"))) +
#
#     geom_nodes(data = g, aes(fill=regulator_score),
#                # fill = "red4",
#                size = 3, colour="black",shape=21) +
#     scale_size_identity() +
#     scale_fill_gradient2(low = "darkblue", mid = "white", high="red4") +
#     monocle3:::monocle_theme_opts() + theme_blank()
#
# }

# To add:
#' #' plot a similar heatmap to lauren's BB testing
#' #'
#' plot_heatmap <- function(ccm) {
#'
#' }
#'
#' plot_heatmap_2 <- function() {
#'
#' }
#'
#' #' plot normalized counts as boxplots
#' plot_boxplots <- function() {
#'
#' }
