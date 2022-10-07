#' Plot a UMAP colored by how cells shift in a given contrast
#'
#' @param ccm A cell_count_model object.
#' @param criterion a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".
#' @return an updated cell_count_model object
#' @export
#' @import ggplot2
#' @import dplyr
plot_contrast <- function(ccm,
                          cond_b_vs_a_tbl,
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
                          plot_edges = c("all", "directed", "undirected", "none"),
                          label_cell_groups = list(),
                          repel_labels = TRUE,
                          model_for_pcors="reduced",
                          switch_label = NULL,
                          sub_cds = NULL,
                          alpha = 1.0,
                          x=1,
                          y=2){

  scale_shifts_by = match.arg(scale_shifts_by)
  plot_labels = match.arg(plot_labels)
  plot_edges = match.arg(plot_edges)

  if (!is.null(sub_cds)) {

    colData(ccm@ccs@cds)[["cell_group"]] = colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
    # get cell groups in sub
    sub_cell_group = colData(sub_cds)[[ccm@ccs@info$cell_group]] %>% unique()

    # subset to select group
    ccm@ccs@cds = ccm@ccs@cds[, colnames(sub_cds)]
    # ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[["cell_group"]] %in% sub_cell_group]
    ccm@ccs@metadata[["cell_group_assignments"]] = ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds),]

    cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% filter(cell_group %in% sub_cell_group)

    # switch coords to the new ones
    sub_umap = reducedDims(sub_cds)[["UMAP"]]
    reducedDims(ccm@ccs@cds)[["UMAP"]] = sub_umap[colnames(ccm@ccs@cds),]

  }

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
                                                              log_abundance_thresh,
                                                              model_for_pcors=model_for_pcors)
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

  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,y]

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
      alpha = alpha,
      stroke = 0
    ) +
    scale_color_gradient2(
      low = "royalblue3",
      mid = "white",
      high = "orangered3",
      na.value = "white",
      limits = fc_limits
    )  +
    #theme_void() +
    #theme(legend.position = "none") +
    monocle3:::monocle_theme_opts()

  # force to be empty
  if (plot_edges == "directed") {
    undirected_edge_df = undirected_edge_df %>% dplyr::filter(!edge_type %in% c("undirected"))
  } else if (plot_edges == "undirected") {
    directed_edge_df = directed_edge_df %>% dplyr::filter(edge_type %in% c("undirected"))
  }

  if (plot_edges != "none") {
    gp = gp  +
      geom_segment(data = undirected_edge_df,
                   aes(x = get(paste0("to_umap_", x)),
                       y = get(paste0("to_umap_", y)),
                       xend = get(paste0("from_umap_", x)),
                       yend = get(paste0("from_umap_", y)),
                       size = edge_size * scaled_weight),
                   color="lightgray") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = get(paste0("to_umap_", x)),
                       y = get(paste0("to_umap_", y)),
                       xend = get(paste0("from_umap_", x)),
                       yend = get(paste0("from_umap_", y)),
                       size = edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = get(paste0("to_umap_", x)),
                       y = get(paste0("to_umap_", y)),
                       xend=(get(paste0("to_umap_", x))+get(paste0("from_umap_", x)))/2,
                       yend = (get(paste0("to_umap_", y))+get(paste0("from_umap_", x)))/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = get(paste0("from_umap_", x)),
                       y = get(paste0("from_umap_", y)),
                       xend = get(paste0("to_umap_", x)),
                       yend = get(paste0("to_umap_", y)),
                       size=edge_size * scaled_weight),
                   color="black") +
      geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = get(paste0("from_umap_", x)),
                       y = get(paste0("from_umap_", y)),
                       xend=(get(paste0("from_umap_", x))+get(paste0("to_umap_", x)))/2,
                       yend = (get(paste0("from_umap_", y))+get(paste0("to_umap_", y)))/2,
                       size=edge_size * scaled_weight),
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_size_identity()

  }

  if (is.null(switch_label) == FALSE) {
    label_df = centroids(ccs, switch_group = switch_label)

    if (length(label_cell_groups) > 0)
      label_df = label_df %>% filter(cell_group %in% label_cell_groups)

    # can't have both
    plot_labels = "none"

    if (repel_labels) {
      gp = gp + ggrepel::geom_label_repel(data = label_df,
                                          mapping = aes(get(paste0("umap_", x)),
                                                        get(paste0("umap_", y)),
                                                        label=cell_group),
                                          size=I(group_label_size),
                                          fill = "white")
    } else {

      gp <- gp + geom_text(data = label_df,
                           mapping = aes(get(paste0("umap_", x)),
                                         get(paste0("umap_", y)),
                                         label=cell_group),
                           size=I(group_label_size))
    }

  }


  if (plot_labels != "none") {
    label_df = umap_centers_delta_abund

    if (plot_labels == "significant")
      label_df = label_df %>% filter(delta_log_abund != 0)

    if (length(label_cell_groups) > 0)
      label_df = label_df %>% filter(cell_group %in% label_cell_groups)

    if (repel_labels) {
      gp = gp + ggrepel::geom_label_repel(data = label_df,
                                          mapping = aes(get(paste0("umap_", x)),
                                                        get(paste0("umap_", y)),
                                                        label=cell_group),
                                          size=I(group_label_size),
                                          fill = "white")
    } else {

      gp <- gp + geom_text(data = label_df,
                           mapping = aes(get(paste0("umap_", x)),
                                         get(paste0("umap_", y)),
                                         label=cell_group),
                           size=I(group_label_size))
    }
  }

  gp = gp + xlab(paste0("UMAP",x)) + ylab(paste0("UMAP",y))
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


my_plot_cells <- function(data,
                          color_cells_by = NULL,
                          cell_size=1,
                          legend_position="none",
                          coefficients = NULL,
                          effect_sizes = NULL,
                          cond_b_vs_a_tbl = NULL,
                          cell_group = NULL,
                          q_value_thresh = 1.0,
                          fc_limits = c(-3,3),
                          plot_labels = TRUE,
                          label_cells_by = NULL,
                          group_label_size=2,
                          repel_labels = TRUE,
                          lab_title = NULL,
                          color_values = NULL,
                          alpha = 1.0,
                          x = 1,
                          y = 2) {

  if (class(data) == "cell_count_set") {
    cds = data@cds
    ccs = data
    plot_df = as.data.frame(colData(cds))
    plot_df$cell_group = ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else if (class(data) == "cell_data_set") {
    cds = data
    plot_df = as.data.frame(colData(cds))
    if (is.null(color_cells_by)) {
      color_cells_by = "cluster"
    }
    plot_df$cell_group = plot_df[[color_cells_by]]
  } else if (class(data) == "cell_count_model") {
    ccs = data@ccs
    cds = data@ccs@cds
    plot_df = as.data.frame(colData(cds))
    plot_df$cell_group = data@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else {
    print("some error message")
  }

  if (is.null(color_cells_by)) {
    color_cells_by = ccs@info$cell_group
  }

  plot_df$cell = row.names(plot_df)
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[plot_df$cell,y]


  # could be an vector that corresponds to a label
  if (!is.null(coefficients)) {
    res_df = data.frame("coefficients" = coefficients) %>% rownames_to_column("cell_group")
    plot_df = plot_df %>% left_join(res_df, by="cell_group")

    if (is.null(color_cells_by))
      color_cells_by = "coefficients"

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
    # theme_void() +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()


  if (color_cells_by == "timepoint") {

    num_colors = unique(plot_df$timepoint) %>% sort() %>% length()
    full_spectrum_timepoint = get_colors(num_colors, type="rainbow")
    names(full_spectrum_timepoint) = unique(plot_df$timepoint) %>% sort()

    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = as.character(timepoint)),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      scale_color_manual(values = full_spectrum_timepoint)

  } else if (color_cells_by %in% c("coefficients", "delta_log_needed")) {
    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = color_cells_by),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = fc_limits
      )


  } else if (grepl("dose", color_cells_by)) {

    # plot_df$color_cells_by


    gp =  gp  +
      # geom_point(
      #   data = plot_df %>% filter(is.na(color_cells_by)),
      #   aes(umap2D_1, umap2D_2),
      #   color = "gray",
      #   size = cell_size,
      #   alpha = 0.2,
      #   stroke = 0) +
      geom_point(
        data = plot_df %>% filter(!is.na(color_cells_by),
                                  color_cells_by > 0),
        aes(umap2D_1, umap2D_2,
            color = log_dose),
        size = cell_size,
        alpha = alpha,
        stroke = 0) +
      # viridis::scale_color_viridis(option = "C")
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = fc_limits
      )


  } else if (color_cells_by %in% c("viridis", "inferno", "C")) {
    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = coefficients),
        size = cell_size,
        alpha = alpha,
        stroke = 0) +
      viridis::scale_color_viridis(option = color_cells_by)

  } else if (color_cells_by == "delta_log_abund") {
    gp = gp  +
      geom_point(
        data = plot_df %>%
          arrange(!is.na(abs(delta_log_abund)),
                  abs(delta_log_abund)),
        aes(umap2D_1, umap2D_2, color = delta_log_abund),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = fc_limits
      )
  } else if (color_cells_by == "none") {
    # do nothing
  }
  else if (is.numeric(plot_df[[color_cells_by]])) {

    # num_colors = unique(plot_df[[color_cells_by]]) %>% sort() %>% length()
    # full_spectrum_timepoint = get_colors(num_colors, type="rainbow")
    # names(full_spectrum_timepoint) = unique(plot_df$timepoint) %>% sort()

    gp = gp  +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
            color = color_cells_by),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      # scale_color_gradientn(colors = full_spectrum_timepoint)+
      viridis::scale_color_viridis(option = "C")

  }
  else {

    if (is.null(color_values)) {
      num_colors = unique(plot_df[[color_cells_by]]) %>% sort() %>% length()
      color_values = get_colors(num_colors, "vibrant")
      names(color_values) = unique(plot_df[[color_cells_by]]) %>% sort()
    }

    gp = gp + geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2,
          color = color_cells_by),
      size = cell_size,
      alpha = alpha,
      stroke = 0
    ) +
      scale_color_manual(values = color_values)
  }

  # add labels

  if (plot_labels) {

    label_df = centroids(ccs)

    # if the ccs is grouped on something else, by default label it like this
    # if it is currently not defined
    if (color_cells_by != ccs@info$cell_group & is.null(label_cells_by))
      label_cells_by = color_cells_by

    if (is.null(label_cells_by) == FALSE) {

      label_df = convert_to_col(ccs, label_df, label_cells_by)

    }


    if (repel_labels) {
      gp = gp + ggrepel::geom_label_repel(data = label_df,
                                          mapping = aes(get(paste0("umap_", x)),
                                                        get(paste0("umap_", y)),
                                                        label=cell_group),
                                          size=I(group_label_size),
                                          fill = "white")
    } else {
      gp = gp + geom_text(data = label_df,
                          mapping = aes(get(paste0("umap_", x)),
                                        get(paste0("umap_", y)),
                                        label=cell_group),
                          size=I(group_label_size))
    }

  }

  if (is.null(lab_title)) {
    lab_title = color_cells_by
  }

  gp = gp + labs(color = lab_title) + xlab(paste0("UMAP",x)) + ylab(paste0("UMAP",y))

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
                      state_graph,
                      edge_size = 1,
                      path_color = "black",
                      color_path_by = NULL,
                      cond_b_vs_a_tbl = NULL,
                      directed = TRUE,
                      x = 1,
                      y = 2,
                      switch_label = NULL,
                      ...) {

  if (class(data) == "cell_count_set") {
    ccs = data
  } else if (class(data) == "cell_count_model") {
    ccs = data@ccs
  }

  # if (is.null(switch_label) & !is.null(unique(igraph::V(state_graph)$cell_group)) & is(state_graph, "igraph")) {
  #   switch_label = unique(igraph::V(state_graph)$cell_group)
  # } else {
  #   print("define a switch label grouping that matches the ccm")
  # }

  if (is(state_graph, "igraph")){
    path_df = state_graph %>% igraph::as_data_frame()
  }else{
    path_df = state_graph
  }

  if (is.null(switch_label)) {
    node_intersection = intersect(rownames(ccs), unique(path_df$from, path_df$to))
  } else {
    node_intersection = intersect(ccs@colData[[switch_label]], unique(path_df$from, path_df$to))
  }

  assertthat::assert_that(

    tryCatch(length(node_intersection) > 0,
             error = function(e) FALSE),
    msg = "specify a label switch, path nodes and ccm cell group do not match"

  )


  gp = my_plot_cells(data,  x=x, y=y, color_cells_by = switch_label, ...)

  if (class(data) == "cell_count_set") {
    umap_centers = centroids(data, switch_group = switch_label)
  } else if (class(data) == "cell_count_model") {
    umap_centers = centroids(data@ccs, switch_group = switch_label)
  }



  path_df = add_umap_coords(path_df, umap_centers)

  if (directed) {
    if (is.null(color_path_by) == FALSE) {

      gp = gp +
        ggnewscale::new_scale_color() +
        geom_segment(data = path_df,
                     aes(x = get(paste0("umap_to_", x)),
                         y = get(paste0("umap_to_", y)),
                         xend=get(paste0("umap_from_", x)),
                         yend = get(paste0("umap_from_", y)),
                         color = get(color_path_by)),
                     size=edge_size ) +
        geom_segment(data = path_df,
                     aes(x = umap_from_1,
                         y = umap_from_2,
                         xend = (get(paste0("umap_to_", x))+get(paste0("umap_from_", x)))/2,
                         yend = (get(paste0("umap_to_", y))+get(paste0("umap_from_", y)))/2,
                         color = get(color_path_by)),
                     size=edge_size,
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
                             aes(x = get(paste0("umap_to_", x)),
                                 y = get(paste0("umap_to_", y)),
                                 xend=get(paste0("umap_from_", x)),
                                 yend = get(paste0("umap_from_", y))),
                             color = path_color,
                             size = edge_size) +
        geom_segment(data = path_df,
                     aes(x = get(paste0("umap_from_", x)),
                         y = get(paste0("umap_from_", y)),
                         xend = (get(paste0("umap_to_", x)) + get(paste0("umap_from_", x)))/2,
                         yend = (get(paste0("umap_to_", y)) + get(paste0("umap_from_", y)))/2),
                     size = edge_size,
                     color = path_color,
                     linejoin='mitre',
                     arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))
    }
  } else {

    if (is.null(color_path_by) == FALSE) {
      gp = gp +
        ggnewscale::new_scale_color() +
        geom_segment(data = path_df,
                     aes(x = get(paste0("umap_to_", x)),
                         y = get(paste0("umap_to_", y)),
                         xend=get(paste0("umap_from_", x)),
                         yend = get(paste0("umap_from_", y)),
                         color = get(color_path_by)),
                     size=edge_size) +
        scale_color_gradient2(
          low = "#122985",
          mid = "white",
          high = "red4",
          na.value = "white"
        )

    } else {
      gp = gp + geom_segment(data = path_df,
                             aes(x = get(paste0("umap_to_", x)),
                                 y = get(paste0("umap_to_", y)),
                                 xend = get(paste0("umap_from_", x)),
                                 yend = get(paste0("umap_from_", y))),
                             color = path_color,
                             size = edge_size)
    }


  }


  return(gp)

}


add_path_edge <- function(gp,
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
  n = network::network(edges %>% dplyr::select(from,to), directed = T, loops = T)
  nodes = nodes[match(network::network.vertex.names(n), nodes$id),]
  n %v% "id" = network::network.vertex.names(n)

  merged_coords = ggnetwork(n) %>% select(id, vertex.names) %>% unique() %>%
    left_join(umap_centers, by = c("id"="cell_group"))
  rownames(merged_coords) = merged_coords$id
  coords = merged_coords[network::network.vertex.names(n),] %>% select(umap_1,umap_2)
  geo = as.matrix(sapply(coords, as.numeric))

  g <- ggnetwork(x = n, layout = geo, arrow.gap=arrow.gap, scale = scale)

  show(ggplot(g, aes(x, y, xend = xend, yend = yend)) +
         geom_edges(arrow = arrow(length = unit(6, "pt"), type="closed")) +
         geom_nodes(size = 7,colour="black",shape=21) +
         geom_nodetext_repel(aes(label = id), size=3) +
         theme_blank())

}

is_concordant <- function(cond_b_v_a_tbl, state_transition_graph) {

  parent_child_foldchanges(state_transition_graph, cond_b_v_a_tbl) %>%
    mutate(concordant = case_when( fold_changes == "parent increase, descendents increase" ~ TRUE,
                                   fold_changes == "parent decrease, descendents decrease" ~ TRUE,
                                   TRUE~ FALSE)) %>% ungroup #%>%
  # group_by(consitent) %>% tally() %>%
  # mutate(pct = n/sum(n))

}

is_discordant <- function(cond_b_v_a_tbl, state_transition_graph) {

  parent_child_foldchanges(state_transition_graph, cond_b_v_a_tbl) %>%
    mutate(discordant = case_when( fold_changes == "parent decrease, descendents no change" ~ TRUE,
                                   fold_changes == "parent decrease, descendents increase" ~ TRUE,
                                   fold_changes == "parent increase, descendents no change" ~ TRUE,
                                   fold_changes == "parent increase, descendents decrease" ~ TRUE,
                                   TRUE~ FALSE)) %>% ungroup

}

collect_psg_node_metadata <- function(ccm,
                                      color_nodes_by,
                                      label_nodes_by,
                                      group_nodes_by)
{
  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs@cds)[,c(color_nodes_by,
                                                label_nodes_by,
                                                group_nodes_by), drop=F] %>%
    as.data.frame
  cell_group_metadata$cell_group = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  if (is.null(color_nodes_by) == FALSE){
    color_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(color_by_metadata) = c("cell_group", "color_nodes_by")
    node_metadata = left_join(node_metadata, color_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE){
    group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
    node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE){
    label_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(label_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    node_metadata = left_join(node_metadata, label_by_metadata, by=c("id"="cell_group"))
  }else{
    node_metadata$label_nodes_by = node_metadata$id
  }
  node_metadata = node_metadata %>% distinct() %>% as.data.frame(stringsAsFactor=FALSE)
  row.names(node_metadata) = node_metadata$id

  return(node_metadata)
}

layout_state_graph <- function(G, node_metadata, edge_labels, weighted=FALSE)
{

  if (weighted){
    G_nel = graph::graphAM(igraph::get.adjacency(G, attr="weight") %>% as.matrix(),
                           edgemode = 'directed', values=list(weight=1)) %>% as("graphNEL")
  } else{
    G_nel = graph::graphAM(igraph::get.adjacency(G) %>% as.matrix(),
                           edgemode = 'directed') %>%as("graphNEL")
  }

  edge_weights = unlist(graph::edgeWeights(G_nel))
  names(edge_weights) = stringr::str_replace_all(names(edge_weights), "\\.", "~")

  make_subgraphs_for_groups <- function(grouping_set, G_nel, node_meta_df){
    nodes = node_meta_df %>% filter(group_nodes_by == grouping_set) %>% pull(id) %>% as.character
    sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel))
    return (sg)
  }

  subgraph_df = node_metadata %>% group_by(group_nodes_by) %>%
    summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
                                    .x = group_nodes_by,
                                    G_nel,
                                    node_metadata))

  if (is.null(edge_labels)== FALSE) {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraph_df$subgraph, edgeAttrs=list(label=edge_labels), recipEdges="distinct")
    label_df = data.frame("x" = gvizl@renderInfo@edges$labelX, "y" = gvizl@renderInfo@edges$labelY) %>%
      rownames_to_column("edge_name") %>%
      left_join(tibble("edge_name" = names(edge_labels), label=edge_labels))# %>% rownames_to_column("edge_name"), by = "edge_name")

  } else {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraph_df$subgraph, recipEdges="distinct")
    label_df=NULL
  }

  gvizl_coords = cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)

  beziers = lapply(gvizl@renderInfo@edges$splines, function(bc) {
    bc_segments = lapply(bc, Rgraphviz::bezierPoints)
    bezier_cp_df = do.call(rbind, bc_segments) %>% as.data.frame
    colnames(bezier_cp_df) = c("x", "y")
    #bezier_cp_df$point = "control"
    #bezier_cp_df$point[1] = "end"
    #bezier_cp_df$point[nrow(bezier_cp_df)] = "end"
    bezier_cp_df
    #control_point_coords = lapply(bc, function(cp) Rgraphviz::getPoints)
    #control_point_coords = rbind(control_point_coords)
  })
  bezier_df = do.call(rbind, beziers)
  bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]
  bezier_df$from = stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[,1]
  bezier_df$to = stringr::str_split_fixed(bezier_df$edge_name, "~", 2)[,2]
  #bezier_df = bezier_df %>% mutate(x = ggnetwork:::scale_safely(x),
  #                                 y = ggnetwork:::scale_safely(y))

  bezier_df = left_join(bezier_df, tibble(edge_name=names(gvizl@renderInfo@edges$direction), edge_direction=gvizl@renderInfo@edges$direction))
  bezier_df = bezier_df %>% dplyr::distinct()
  return(list(gvizl_coords=gvizl_coords, bezier_df=bezier_df, label_df=label_df))
}

#' Plot the state graph with annotations
#' @export
plot_state_graph_annotations <- function(ccm,
                                         state_graph,
                                         color_nodes_by=NULL,
                                         label_nodes_by=NULL,
                                         group_nodes_by=NULL,
                                         label_edges_by=NULL,
                                         edge_weights=NULL,
                                         arrow.gap=0.03,
                                         arrow_unit = 2,
                                         bar_unit = .075,
                                         node_size = 2,
                                         num_layers=10,
                                         min_edge_size=0.1,
                                         max_edge_size=2,
                                         fract_expr = 0.0,
                                         mean_expr = 0.0,
                                         unlabeled_groups = c("Unknown"),
                                         label_groups=TRUE,
                                         hide_unlinked_nodes=TRUE,
                                         group_label_font_size=6,
                                         edge_label_font_size=2,
                                         label_conn_linetype="dotted",
                                         legend_position = "none",
                                         con_colour = "darkgrey",
                                         group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = collect_psg_node_metadata(ccm, color_nodes_by, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
  }

  if (is.null(label_edges_by)){
    edge_info = edges %>% select(from, to)
  }else{
    if (is(state_graph, "igraph")){
      edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    }else{
      edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

    edges = edges %>% left_join(edge_info)
    #print(edges)
  }

  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE){
    G_df = igraph::as_data_frame(G)
    edge_names =  stringr::str_c(G_df$from, G_df$to, sep="~")
    edge_labels = igraph::E(G)$label
    names(edge_labels) = edge_names
    #print(edge_labels)
    #edge_labels = NULL
  }else{
    edge_labels=NULL
  }

  layout_info = layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size)
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=group_label_font_size,
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
      if (label_groups){
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            label = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
                                        label.fontsize=group_label_font_size,
                                        label.fontface="plain",
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }else{
        p = p + ggforce::geom_mark_rect(aes(x, y,
                                            fill = group_nodes_by,
                                            filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                        size=0,
                                        label.fontsize=group_label_font_size,
                                        con.linetype=label_conn_linetype,
                                        con.colour=con_colour,
                                        data=g)
      }

    }

  }

  if (is.null(color_nodes_by) == FALSE) {
    if (is.null(label_nodes_by) == FALSE){
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(data = g,
                                    aes(x, y,
                                        fill = !!sym(color_nodes_by),
                                        label = label_nodes_by),
                                    size = node_size) +
          labs(fill = color_nodes_by)
        p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodelabel(data = g,
                                          aes(x, y,
                                              fill = color_nodes_by,
                                              label = label_nodes_by),
                                          size = node_size) +
          labs(fill = color_nodes_by)

      }
    }else{
      # if numerical
      if (is.numeric(g[[color_nodes_by]])) {
        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodes(data = g,
                                    aes(x, y,
                                        color = !!sym(color_nodes_by)),
                                    size = node_size) +
          labs(color = color_nodes_by)
        p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }
      else {
        # if categorical
        p = p + ggnetwork::geom_nodes(data = g,
                                          aes(x, y,
                                              color = color_nodes_by),
                                          size = node_size) +
          labs(color = color_nodes_by)
      }
    }

  } else {
    if (is.null(label_nodes_by) == FALSE){
    p = p + ggnetwork::geom_nodelabel(data = g,
                                      aes(x, y, xend = xend, yend = yend,
                                          label = label_nodes_by),
                                      size = node_size)
    }else{
      p = p + ggnetwork::geom_nodes(data = g,
                                        aes(x, y, xend = xend, yend = yend),
                                        size = node_size)
    }
  }

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    p = p + geom_text(data = label_df,
                      aes(x,y, label = label),
                      size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}

#' This function plots when each cell state is lost following a perturbation
#'
#' FIXME: this is really a prototype and needs lots of cleanup and refactoring.
#' It could plot the earliest time a loss is detected, the latest, etc. But
#' right now it only plots the loss at the time of peak abundance in the control
#' see estimate_loss_timing() for more details.
#'
#' @export
plot_state_graph_losses <- function(perturbation_ccm,
                                    state_graph,
                                    start_time,
                                    stop_time,
                                    interval_step,
                                    interval_col,
                                    log_abund_detection_thresh,
                                    q_val,
                                    label_nodes_by=NULL,
                                    group_nodes_by=NULL,
                                    label_edges_by=NULL,
                                    edge_weights=NULL,
                                    arrow.gap=0.03,
                                    arrow_unit = 2,
                                    bar_unit = .075,
                                    node_size = 2,
                                    num_layers=10,
                                    edge_size=2,
                                    fract_expr = 0.0,
                                    mean_expr = 0.0,
                                    unlabeled_groups = c("Unknown"),
                                    hide_unlinked_nodes=TRUE,
                                    label_font_size=6,
                                    label_conn_linetype="dotted",
                                    legend_position = "none",
                                    group_outline=FALSE,
                                    ...)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = hooke:::collect_psg_node_metadata(perturbation_ccm, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(edge_weights)){
    edges = edges %>% select(from, to)
  }else{
    edges = edges %>% select(from, to, weight=!!sym(edge_weights))
  }

  if (is.null(label_edges_by)){
    edge_info = edges %>% select(from, to)
  }else{
    if (is(state_graph, "igraph")){
      edge_info = state_graph %>% igraph::as_data_frame() %>% select(from, to, label=!!sym(label_edges_by))
    }else{
      edge_info = state_graph %>% select(from, to, label=!!sym(label_edges_by))
    }

    edges = edges %>% left_join(edge_info)
    #print(edges)
  }

  earliest_loss_tbl = hooke:::estimate_loss_timing(perturbation_ccm,
                                                   start_time=start_time,
                                                   stop_time=stop_time,
                                                   interval_step = interval_step,
                                                   interval_col=interval_col,
                                                   log_abund_detection_thresh=log_abund_detection_thresh,
                                                   q_val = q_val,
                                                   ...)
  #print (earliest_loss_tbl)
  node_metadata = node_metadata %>% left_join(earliest_loss_tbl, by=c("id" ="cell_group"))

  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  if (is.null(igraph::E(G)$label) == FALSE){
    G_df = igraph::as_data_frame(G)
    edge_names =  stringr::str_c(G_df$from, G_df$to, sep="~")
    edge_labels = igraph::E(G)$label
    names(edge_labels) = edge_names
    #print(edge_labels)
    edge_labels = NULL
  }else{
    edge_labels=NULL
  }

  layout_info = hooke:::layout_state_graph(G, node_metadata, edge_labels, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = edge_size * edge_score)
  }else{
    bezier_df$edge_thickness = edge_size
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  print (g)
  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness), data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){


    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=label_font_size,
                                      con.linetype=label_conn_linetype,
                                      data=g)
    } else {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          fill = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0,
                                      label.fontsize=label_font_size,
                                      con.linetype=label_conn_linetype,
                                      data=g)
    }

  }


  # if numerical

  p = p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodelabel(data = g,
                              aes(x, y,
                                  fill = peak_loss_time,
                                  label = label_nodes_by),
                              size = node_size)
  p = p + scale_fill_stepsn(n.breaks=5, colours = terrain.colors(5)) #+ scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")





  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$label_df
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}


num_extract <- function(string, as_char = TRUE){

  if (as_char)
    stringr::str_trim(
      format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
             scientific = FALSE,
             trim = TRUE)
    )

  else
    as.numeric(
      stringr::str_trim(
        format(stringr::str_extract(string, "\\-*\\d+\\.*\\d*"),
               scientific = FALSE,
               trim = TRUE)
      )
    )

}

calc_sig_ind <- function(p_value, html = TRUE) {

  #p_value <- suppressWarnings(
  #  num_extract(p_value, as_char = FALSE)
  #)

  if (html) {
    dplyr::case_when(
      p_value <= 0 ~ "",
      p_value <= 0.001 ~ "\\***",
      p_value <= 0.01 ~ "\\**",
      p_value <= 0.05 ~ "\\*",
      p_value <= 0.1 ~ ".",
      p_value <= 1 ~ "",
      TRUE ~ ""
    )
  } else {
    dplyr::case_when(
      p_value <= 0 ~ "",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      p_value <= 0.1 ~ ".",
      p_value <= 1 ~ "",
      TRUE ~ ""
    )
  }

}

#' Plot the expression of a single gene on the state transition graph
#' @export
plot_state_graph_abundance_changes <- function(ccm,
                                               state_graph,
                                               comp_abund_table,
                                               contrast = "contrast",
                                               label_nodes_by=NULL,
                                               group_nodes_by=NULL,
                                               edge_labels=NULL,
                                               fc_limits=c(-3,3),
                                               arrow.gap=0.03,
                                               arrow_unit = 2,
                                               bar_unit = .075,
                                               node_size = 6,
                                               num_layers=10,
                                               edge_size=0.5,
                                               fract_expr = 0.0,
                                               mean_expr = 0.0,
                                               unlabeled_groups = c("Unknown"),
                                               label_subset=NULL,
                                               hide_unlinked_nodes=TRUE,
                                               label_font_size=6,
                                               label_conn_linetype="dotted",
                                               legend_position = "none",
                                               group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  comp_abund_table[["contrast"]] = comp_abund_table[[contrast]]
  comp_abund_table$contrast = as.factor(comp_abund_table$contrast)

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = collect_psg_node_metadata(ccm, color_nodes_by=NULL, label_nodes_by=label_nodes_by, group_nodes_by=group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }


  G = edges %>% select(from, to) %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  layout_info = layout_state_graph(G, node_metadata, edge_labels)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)



  #gene_ids = rowData(ccm@ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()

  #gene_expr = aggregated_expr_data(ccm@ccs@cds[gene_ids,], group_cells_by = ccm@ccs@info$cell_group)
  #sub_gene_expr = gene_expr %>%
  #  filter(gene_short_name %in% genes)

  #method = get(method)
  # sub_gene_expr = sub_gene_expr %>%
  #   group_by(cell_group) %>%
  #   dplyr::summarise(fraction_expressing = method(fraction_expressing),
  #                    mean_expression = method(mean_expression),
  #                    specificity = method(specificity)) %>%
  #   mutate(gene_short_name =  paste(genes, collapse = "-"))


  #node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

  abund_fc_df = comp_abund_table

  #sub_gene_expr = sub_gene_expr %>%
  #  group_by(gene_short_name) %>%
  #  mutate(
  #    fraction_max = mean_expression / max(mean_expression),
  #    gene_expr = case_when(
  #      fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
  #      TRUE ~ FALSE))

  if (is.null(fc_limits)) {
    fc_limits = range(abund_fc_df$delta_log_abund)
  } else {
    min = fc_limits[1]
    max = fc_limits[2]
    abund_fc_df = abund_fc_df %>%
      mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
      mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  }

  abund_fc_df = abund_fc_df %>%
    mutate(
      delta_q_value = pmax(0.0001, delta_q_value),
      q_value_sig_code = calc_sig_ind(delta_q_value, html=FALSE))

  g = left_join(g, abund_fc_df, by=c("name"="cell_group"))

  #g$gene_short_name = factor(g$gene_short_name, levels=genes)
  color_nodes_by = "delta_log_abund"
  # group_outline = TRUE

  p <- ggplot(aes(x,y), data=g)
  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)

  if (is.null(label_subset)) {
    label_subset = unique(node_metadata$group_nodes_by)
  }

  label_subset = label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE){

    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          fill = contrast,
                                          color = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      expand = unit(2.5, "mm"),
                                      radius = unit(1, "mm"),
                                      label.fontsize=label_font_size,
                                      label.buffer = unit(1, "mm"),
                                      con.linetype=label_conn_linetype,
                                      show.legend = F) +
        ggforce::geom_mark_rect(aes(x, y,
                                    fill = contrast,
                                    color = group_nodes_by,
                                    label = group_nodes_by,
                                    filter = group_nodes_by %in% label_subset),
                                size=0.5,
                                expand = unit(2.5, "mm"),
                                radius = unit(1, "mm"),
                                label.fontsize=label_font_size,
                                label.buffer = unit(1, "mm"),
                                con.linetype=label_conn_linetype,
                                show.legend = F)
    } else {
      # p = p + ggnetwork::geom_nodetext_repel(data = g %>% filter(group_nodes_by %in% label_subset),
      #                                       mapping = aes(x, y,
      #                                           label = group_nodes_by),
      #                                           color=I("black"),
      #                                       size=label_font_size/2)
      values = rep("white",length(label_subset))
      names(values) = label_subset

      p = p +
        ggforce::geom_mark_rect(aes(x, y,
                                    fill = contrast,
                                    color = group_nodes_by,
                                    label = group_nodes_by,
                                    filter = group_nodes_by %in% label_subset),
                                size=0.5,
                                expand = unit(2.5, "mm"),
                                radius = unit(1, "mm"),
                                label.fontsize=label_font_size,
                                label.buffer = unit(1, "mm"),
                                con.linetype=label_conn_linetype,
                                show.legend = F) +
        scale_color_manual(values = c(values))

    }

    p = p + scale_fill_manual(values=rep("white", length(unique(g$contrast))))

  }
  p = p + guides(fill = "none")
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              size = -log10(delta_q_value)*1.2),
                          color=I("black")) +
    ggnetwork::geom_nodes(data = g,
                          aes(x, y,
                              size = -log10(delta_q_value),
                              color=delta_log_abund)) +
    ggnetwork::geom_nodetext(data = g,
                             aes(x, y,
                                 label = q_value_sig_code),
                             color=I("black"))
  # ggnetwork::geom_nodetext_repel(data = g,
  #                                aes(x, y,
  #                                    label = q_value_sig_code),
  #                                color=I("black"))
  #labs(fill = color_nodes_by)

  p = p + scale_color_gradient2(low = "royalblue3", mid = "white", high="orangered3")

  if (is.null(edge_labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  p = p + facet_wrap(~contrast)

  p = p + scale_size(range=c(1, node_size)) +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position) + guides(fill = "none", size="none")
  return(p)
}

#' Plot the expression of a single gene on the state transition graph
#' @export
plot_state_graph_gene_expression <- function(ccm,
                                             state_graph,
                                             genes,
                                             method = "min",
                                             gene_expr = NULL,
                                             label_nodes_by=NULL,
                                             group_nodes_by=NULL,
                                             edge_labels=NULL,
                                             arrow.gap=0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             node_size = 6,
                                             num_layers=10,
                                             edge_size=0.5,
                                             fract_expr = 0.0,
                                             mean_expr = 0.0,
                                             unlabeled_groups = c("Unknown"),
                                             scale_to_range = FALSE,
                                             hide_unlinked_nodes=TRUE,
                                             switch_label=NULL,
                                             label_font_size=6,
                                             label_conn_linetype="dotted",
                                             legend_position = "none",
                                             group_outline=FALSE)
{

  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  }else{
    edges = state_graph
  }

  # # check that the
  # if (is.null(switch_label)) {
  #   node_intersection = intersect(rownames(ccm@ccs), unique(edges$from, edges$to))
  # } else {
  #   node_intersection = intersect(ccm@ccs@colData[[switch_label]], unique(edges$from, edges$to))
  # }
  # assertthat::assert_that(
  #
  #   tryCatch(length(node_intersection) > 0,
  #            error = function(e) FALSE),
  #   msg = "specify a label switch, state_graph and ccm cell group do not match"
  #
  # )


  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  node_metadata = collect_psg_node_metadata(ccm, color_nodes_by=NULL, label_nodes_by=label_nodes_by, group_nodes_by=group_nodes_by)

  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }


  G = edges %>% select(from, to) %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  layout_info = layout_state_graph(G, node_metadata, edge_labels)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  gene_ids = rowData(ccm@ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()

  if (is.null(gene_expr)) {
    gene_expr = aggregated_expr_data(ccm@ccs@cds[gene_ids,], group_cells_by = ccm@ccs@info$cell_group)
  }

  sub_gene_expr = gene_expr %>%
    filter(gene_short_name %in% genes)

  method = get(method)
  # sub_gene_expr = sub_gene_expr %>%
  #   group_by(cell_group) %>%
  #   dplyr::summarise(fraction_expressing = method(fraction_expressing),
  #                    mean_expression = method(mean_expression),
  #                    specificity = method(specificity)) %>%
  #   mutate(gene_short_name =  paste(genes, collapse = "-"))


  #node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

  sub_gene_expr = sub_gene_expr %>%
    group_by(gene_short_name) %>%
    mutate(
      fraction_max = mean_expression / max(mean_expression),
      gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE))

  color_nodes_by = "mean_expression"
  expression_legend_label = "mean expression"

  if (scale_to_range) {
    sub_gene_expr = sub_gene_expr %>%
      mutate(value = mean_expression) %>%
      group_by(gene_short_name) %>%
      dplyr::mutate(max_val_for_feature = max(value),
                    min_val_for_feature = min(value)) %>%
      dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
    expression_legend_label = "% Max"
    color_nodes_by = "value"
  }

  g = left_join(g, sub_gene_expr, by=c("name"="cell_group"))

  g$gene_short_name = factor(g$gene_short_name, levels=genes)

  # group_outline = TRUE

  p <- ggplot(aes(x,y), data=g)
  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){
    p = p + ggforce::geom_mark_rect(aes(x, y,
                                        fill = gene_short_name,
                                        color = group_nodes_by,
                                        label = group_nodes_by,
                                        filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                    size=0.5,
                                    expand = unit(2.5, "mm"),
                                    radius = unit(1, "mm"),
                                    label.fontsize=label_font_size,
                                    label.buffer = unit(1, "mm"),
                                    con.linetype=label_conn_linetype)
    p = p + scale_fill_manual(values=rep("white", length(unique(g$gene_short_name))))

  }
  p = p + guides(fill = "none")
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max,
                              color = get(color_nodes_by))) +
    labs(color = expression_legend_label)

  p = p + viridis::scale_color_viridis(option = "viridis")

  if (is.null(edge_labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  p = p + facet_wrap(~gene_short_name)

  p = p + scale_size(labels = scales::percent, range=c(1, node_size)) +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position) + guides(fill = "none")
  return(p)
}

#' Plot a graph that summarized an assembled state transition graph
#'
#' @export
plot_state_transition_graph <- function(ccm,
                                        edges,
                                        # contract_nodes_by=NULL,
                                        color_nodes_by=NULL,
                                        label_nodes_by=NULL,
                                        group_nodes_by=NULL,
                                        layer_nodes_by=NULL,
                                        cond_b_v_a_tbl = NULL,
                                        fc_limits = NULL,
                                        genotype_models_tbl = NULL,
                                        concordant = TRUE,
                                        genes = list(),
                                        labels = NULL,
                                        qval_threshold = 0.05,
                                        arrow.gap=0.03,
                                        arrow_unit = 2,
                                        bar_unit = .075,
                                        node_size = 2,
                                        num_layers=10,
                                        edge_size=0.5,
                                        fract_expr = 0.0,
                                        mean_expr = 0.0,
                                        method = "min",
                                        unlabeled_groups = c("Unknown"),
                                        hide_unlinked_nodes=TRUE,
                                        label_font_size=6,
                                        label_conn_linetype="dotted",
                                        legend_position = "none",
                                        group_outline=FALSE){


  # if (is.null(contract_nodes_by) == FALSE) {
  #   state_graph = edges %>% igraph::graph_from_data_frame()
  #   edges = contract_state_graph(ccm, state_graph, group_nodes_by = contract_nodes_by) %>%
  #     igraph::as_data_frame()
  #   ccm = contract_ccm(ccm, group_nodes_by = contract_nodes_by)
  #
  # }


  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs@cds)[,c(color_nodes_by,
                                                label_nodes_by,
                                                group_nodes_by,
                                                layer_nodes_by), drop=F] %>%
    as.data.frame
  cell_group_metadata$cell_group = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)

  if (is.null(color_nodes_by) == FALSE){
    color_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(color_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(color_by_metadata) = c("cell_group", "color_nodes_by")
    node_metadata = left_join(node_metadata, color_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(group_nodes_by) == FALSE){
    group_by_metadata = cell_group_metadata[,c("cell_group", group_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(group_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(group_by_metadata) = c("cell_group", "group_nodes_by")
    node_metadata = left_join(node_metadata, group_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(layer_nodes_by) == FALSE){
    layer_by_metadata = cell_group_metadata[,c("cell_group", layer_nodes_by)] %>%
      as.data.frame
    if (is.numeric(cell_group_metadata[,c(layer_nodes_by)])){
      layer_by_metadata = layer_by_metadata %>%
        group_by(cell_group) %>%
        summarize(mean_layer_var = mean(!!sym(layer_nodes_by), na.rm=TRUE)) %>%
        mutate(layer = ntile(desc(mean_layer_var),num_layers)) %>% dplyr::select(-mean_layer_var)
    }else{
      layer_by_metadata = layer_by_metadata %>%
        count(cell_group, !!sym(layer_nodes_by)) %>%
        group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    }
    colnames(layer_by_metadata) = c("cell_group", "layer_nodes_by")
    node_metadata = left_join(node_metadata, layer_by_metadata, by=c("id"="cell_group"))
  }
  if (is.null(label_nodes_by) == FALSE){
    label_by_metadata = cell_group_metadata[,c("cell_group", color_nodes_by)] %>%
      as.data.frame %>%
      count(cell_group, !!sym(label_nodes_by)) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    node_metadata = left_join(node_metadata, label_by_metadata, by=c("id"="cell_group"))
  }else{
    node_metadata$label_nodes_by = node_metadata$id
  }
  node_metadata = node_metadata %>% distinct() %>% as.data.frame(stringsAsFactor=FALSE)
  row.names(node_metadata) = node_metadata$id
  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  if (is.null(cond_b_v_a_tbl) == FALSE) {
    cond_b_v_a_tbl = cond_b_v_a_tbl %>%
      mutate(delta_log_abund = ifelse(delta_q_value < qval_threshold, delta_log_abund, 0))

    if (is.null(fc_limits) == FALSE) {
      min = fc_limits[1]
      max = fc_limits[2]
      cond_b_v_a_tbl = cond_b_v_a_tbl %>%
        mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
        mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
    }

    node_metadata = node_metadata %>% left_join(cond_b_v_a_tbl, by = c("id" = "cell_group"))

    # might be better way to do this
    color_nodes_by = "delta_log_abund"
    #idk what to do about this currently
    group_outline = TRUE
  }


  if (length(genes) > 0) {

    # ensembl_ids = rowData(ccm@ccs@cds) %>% as.data.frame %>% filter(gene_short_name %in% genes) %>% rownames()

    gene_expr = aggregated_expr_data(ccm@ccs@cds, group_cells_by = ccm@ccs@info$cell_group)
    sub_gene_expr = gene_expr %>%
      filter(gene_short_name %in% genes)

    method = get(method)
    sub_gene_expr = sub_gene_expr %>%
      group_by(cell_group) %>%
      dplyr::summarise(fraction_expressing = method(fraction_expressing),
                       mean_expression = method(mean_expression),
                       specificity = method(specificity)) %>%
      mutate(gene_short_name =  paste(genes, collapse = "-"))

    node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

    node_metadata = node_metadata %>%
      mutate(gene_expr = case_when(
        fraction_expressing >= fract_expr & mean_expression >= mean_expr ~ TRUE,
        TRUE ~ FALSE))

    color_nodes_by = "mean_expression"
    group_outline = TRUE
  }

  G = edges %>% select(from, to) %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)

  # run sugiyama layout
  layers = NULL
  if (is.null(layer_nodes_by) == FALSE) {
    layers=igraph::V(G)$layer_nodes_by
  }


  #lay1 <- igraph::layout_with_sugiyama(G, layers=layers, maxiter=1000)
  #sg1 = test_graph %>% graph::subGraph(snodes=c("13", "9", "10", "1", "33"))
  #subGList <- vector(mode="list", length=1)
  #subGList[[1]] <- list(graph=sg1)
  #subGList[[2]] <- list(graph=sg2, cluster=FALSE)
  #subGList[[3]] <- list(graph=sg3)


  G_nel = graph::graphAM(igraph::get.adjacency(G) %>% as.matrix(),
                         edgemode = 'directed') %>%
    as("graphNEL")

  make_subgraphs_for_groups <- function(grouping_set, G_nel, node_meta_df){
    nodes = node_meta_df %>% filter(group_nodes_by == grouping_set) %>% pull(id) %>% as.character
    sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel))
    return (sg)
  }

  subgraph_df = node_metadata %>% group_by(group_nodes_by) %>%
    summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
                                    .x = group_nodes_by,
                                    G_nel,
                                    node_metadata))

  if (is.null(labels)== FALSE) {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraph_df$subgraph, edgeAttrs=list(label=labels))
    label_df = data.frame("x" = gvizl@renderInfo@edges$labelX, "y" = gvizl@renderInfo@edges$labelY) %>%
      rownames_to_column("edge_name") %>%
      left_join(data.frame("label" = labels) %>% rownames_to_column("edge_name"), by = "edge_name")

  } else {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraph_df$subgraph)
  }

  gvizl_coords = cbind(gvizl@renderInfo@nodes$nodeX, gvizl@renderInfo@nodes$nodeY)


  beziers = lapply(gvizl@renderInfo@edges$splines, function(bc) {
    bc_segments = lapply(bc, Rgraphviz::bezierPoints)
    bezier_cp_df = do.call(rbind, bc_segments) %>% as.data.frame
    colnames(bezier_cp_df) = c("x", "y")
    #bezier_cp_df$point = "control"
    #bezier_cp_df$point[1] = "end"
    #bezier_cp_df$point[nrow(bezier_cp_df)] = "end"
    bezier_cp_df
    #control_point_coords = lapply(bc, function(cp) Rgraphviz::getPoints)
    #control_point_coords = rbind(control_point_coords)
  })
  bezier_df = do.call(rbind, beziers)
  bezier_df$edge_name = stringr::str_split_fixed(row.names(bezier_df), "\\.", 2)[,1]
  #bezier_df = bezier_df %>% mutate(x = ggnetwork:::scale_safely(x),
  #                                 y = ggnetwork:::scale_safely(y))


  if (is.null(genotype_models_tbl) == FALSE) {

    state_transition_graph = edges %>% igraph::graph_from_data_frame()
    if (concordant) {
      edge_count = genotype_models_tbl %>%
        mutate("concordant" = purrr::map(.f = is_concordant,
                                         .x = genotype_eff,
                                         state_transition_graph)) %>%
        select(gene_target, concordant) %>%
        tidyr::unnest(c(concordant)) %>%
        select(gene_target, parent, child, concordant) %>%
        group_by(parent, child, concordant) %>%
        tally() %>% ungroup() %>%
        filter(concordant) %>%
        mutate(edge_name = paste(parent,child, sep="~"))
    } else {
      edge_count = genotype_models_tbl %>%
        mutate("discordant" = purrr::map(.f = is_discordant,
                                         .x = genotype_eff,
                                         state_transition_graph)) %>%
        select(gene_target, discordant) %>%
        tidyr::unnest(c(discordant)) %>%
        select(gene_target, parent, child, discordant) %>%
        group_by(parent, child, discordant) %>%
        tally() %>% ungroup() %>%
        filter(discordant) %>%
        mutate(edge_name = paste(parent,child, sep="~"))
    }



    bezier_df = left_join(bezier_df, edge_count, by = "edge_name")
    bezier_df$n = tidyr::replace_na(bezier_df$n, -1)
  }

  # if (is.null(edge_count) == FALSE) {
  #   bezier_df = left_join()
  #   edge_weights = state_transition_graph %>%
  #     igraph::as_data_frame() %>%
  #     left_join(edge_count, by = c("from" = "parent", "to" = "child")) %>%
  #     filter(concordant) %>%
  #     mutate(edge_name = paste(from, to, sep="~"))
  #
  #   bezier_df = left_join(bezier_df, edge_weights, by = "edge_name")
  #   bezier_df$n = replace_na(bezier_df$n, -1)
  #
  #   bezier_df %>% mutate("" = ifelse())
  #
  # }


  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)


  # add level information
  #g = g %>% left_join(level_df %>% rownames_to_column("id"), by = c("vertex.names"="id"))
  #g = g %>% left_join(regulator_score_df, by = c("vertex.names" = "gene_id") )


  # p <- ggplot(mapping = aes(x, y, xend = xend, yend = yend)) +
  #   # draw activator edges
  #   ggnetwork::geom_edges(data = g,
  #                         arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))


  if (is.null(genotype_models_tbl) == FALSE) {

    if (concordant) {
      edge_colors = c("black", "orangered3")
      names(edge_colors) = c(-1, 1) %>% as.factor()
    } else {
      edge_colors = c("black", "royalblue3")
      names(edge_colors) = c(-1, 1) %>% as.factor()
    }


    p <- ggplot() +
      ggplot2::geom_path(aes(x, y, group=edge_name,
                             size = abs(n)/5, color = as.factor(sign(n))), data=bezier_df,
                         arrow = arrow(length = unit(arrow_unit, "pt"), type="closed")) +
      scale_color_manual(values = edge_colors)
  }
  else {
    p <- ggplot() +
      ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  }




  # draw activator edges
  #ggforce::geom_bezier(aes(x = x, y = y, group=edge_name, linetype = "cubic"),
  #                     data = bezier_df)
  if (is.null(group_nodes_by) == FALSE){

    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      label.fontsize=label_font_size,
                                      con.linetype=label_conn_linetype,
                                      data=g)
    } else {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          fill = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0,
                                      label.fontsize=label_font_size,
                                      con.linetype=label_conn_linetype,
                                      data=g)
    }

  }

  if (is.null(color_nodes_by) == FALSE) {

    # if numerical
    if (is.numeric(g[[color_nodes_by]])) {

      if (color_nodes_by == "mean_expression") {

        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(data = g %>% filter(gene_expr),
                                    aes(x, y,
                                        fill = !!sym(color_nodes_by),
                                        label = label_nodes_by),
                                    size = node_size) +
          ggnetwork::geom_nodelabel(data = g %>% filter(!gene_expr),
                                    aes(x, y,
                                        label = label_nodes_by),
                                    fill="white",
                                    size = node_size) +
          labs(fill = color_nodes_by)

        p = p + viridis::scale_fill_viridis(option = "viridis")
      } else {

        p = p + ggnewscale::new_scale_fill() +
          ggnetwork::geom_nodelabel(data = g,
                                    aes(x, y,
                                        fill = !!sym(color_nodes_by),
                                        label = label_nodes_by),
                                    size = node_size) +
          labs(fill = color_nodes_by)
        p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
      }

    }
    else {
      # if categorical
      p = p + ggnetwork::geom_nodelabel(data = g,
                                        aes(x, y,
                                            fill = color_nodes_by,
                                            label = label_nodes_by),
                                        size = node_size) +
        labs(fill = color_nodes_by)

    }

  } else {
    p = p + ggnetwork::geom_nodelabel(data = g,
                                      aes(x, y, xend = xend, yend = yend,
                                          label = label_nodes_by),
                                      size = node_size)
  }

  if (is.null(labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p + scale_size_identity() +
    monocle3:::monocle_theme_opts() +
    ggnetwork::theme_blank() +
    theme(legend.position=legend_position)
  return(p)
}

#' Simple function for inspecting origins for each cell state
#'
#' @export
plot_origins <- function(data,
                         path_df = path_df,
                         edge_size=2,
                         path_color = "black",
                         color_cells_by = "cluster",
                         color_path_by = NULL,
                         cond_b_vs_a_tbl = NULL) {



  if (class(data) == "cell_count_set") {
    umap_centers = centroids(data)
  } else if (class(data) == "cell_count_model") {
    umap_centers = centroids(data@ccs)
  }

  path_df = path_df %>% dplyr::select(destination=cell_group, origin, to, from, origin_score, weight, distance_from_root) %>%
    hooke:::add_umap_coords(umap_centers)

  print (path_df)
  if (is.null(color_path_by) == FALSE) {

    gp = ggplot(data=path_df) +
      ggnewscale::new_scale_color() +
      geom_segment(aes(x = umap_to_1,
                       y = umap_to_2,
                       xend=umap_from_1,
                       yend = umap_from_2,
                       color = get(color_path_by)),
                   size=edge_size) +
      geom_segment(aes(x = umap_from_1,
                       y = umap_from_2,
                       xend = (umap_to_1+umap_from_1)/2,
                       yend = (umap_to_2+umap_from_2)/2,
                       color = get(color_path_by)),
                   size=edge_size,
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white"
      )

  } else {
    gp = ggplot(data = path_df) +
      geom_segment(aes(x = umap_to_1,
                       y = umap_to_2,
                       xend=umap_from_1,
                       yend = umap_from_2),
                   color = path_color,
                   size = edge_size ) +
      geom_segment(aes(x = umap_from_1,
                       y = umap_from_2,
                       xend = (umap_to_1+umap_from_1)/2,
                       yend = (umap_to_2+umap_from_2)/2),
                   size= edge_size,
                   color=path_color,
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
      geom_text(aes(umap_from_1, umap_from_2, label=from), color="blue", data=path_df %>% filter (from == origin)) +
      geom_text(aes(umap_to_1, umap_to_2, label=to), color="red", data=path_df %>% filter (to == destination))
  }

  return(gp)

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

# '
#' @param ccs
#' @param x_col
#' @param normalize
#' @param plot_zeroes
#' @param plot_points
#' @param log_scale
#' @param nrow
#' @param legend_position
plot_cells_per_sample = function(ccs,
                                 x_col,
                                 normalize = F,
                                 plot_zeroes = F,
                                 plot_points = F,
                                 log_scale = F,
                                 nrow = 1,
                                 legend_position = "none") {

  ccs_coldata = colData(ccs) %>% as.data.frame
  ccs_coldata$SF = size_factors(ccs)
  count_df = counts(ccs) %>%
    as.matrix %>%
    as.data.frame() %>%
    rownames_to_column("cell_group") %>%
    tidyr::pivot_longer(-cell_group, names_to = "sample", values_to = "count") %>%
    left_join(ccs_coldata, by = "sample")

  if (normalize) {
    count_df = count_df %>% mutate(count = round(count/SF))
  }

  if (plot_zeroes) {
    count_df = count_df %>% mutate(count = count + 0.001)
  }

  count_df[[x_col]] = as.factor(count_df[[x_col]])

  p = count_df %>%
    ggplot(aes( x = !!sym(x_col), y = count, fill=cell_group)) +
    geom_boxplot() +
    facet_wrap(~cell_group, nrow=nrow) +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()

  if (plot_points) {
    p = p + geom_jitter(aes( x = !!sym(x_col), y = count), size=1)
  }

  if (log_scale) {
    p = p + scale_y_log10()
  }

  return(p)

}


plot_cells_highlight = function(ccs, group_to_highlight, colname) {

  plot_df = as.data.frame(colData(cds))
  plot_df$cell_group = colData(cds)[[colname]]

  plot_df$cell = row.names(plot_df)
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[plot_df$cell,y]

  gp = ggplot() +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2),
      color = "black",
      size = 1.5 * cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>% filter(!cell_group %in% c(group_to_highlight)),
      aes(umap2D_1, umap2D_2),
      color = "white",
      size = cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>% filter(cell_group %in% c(group_to_highlight)),
      aes(umap2D_1, umap2D_2),
      color = "red",
      size = cell_size,
      stroke = 0
    ) +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()

  return(gp)

}


plot_contrast_3d <- function(ccm,
                             cond_b_vs_a_tbl,
                             log_abundance_thresh = -5,
                             scale_shifts_by=c("receiver", "sender", "none"),
                             edge_size=2,
                             cell_size=25,
                             q_value_thresh = 1.0,
                             group_label_size=2,
                             plot_labels = c("significant", "all", "none"),
                             fc_limits=c(-3,3),
                             sender_cell_groups=NULL,
                             receiver_cell_groups=NULL,
                             plot_edges = c("all", "directed", "undirected", "none"),
                             label_cell_groups = list(),
                             repel_labels = TRUE,
                             model_for_pcors="reduced",
                             switch_label = NULL,
                             sub_cds = NULL,
                             alpha = 1.0,
                             x=1,
                             y=2) {


  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)

  plot_df$umap3D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap3D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]
  plot_df$umap3D_3 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,3]

  cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))


  plot_df = dplyr::left_join(plot_df,
                             cond_b_vs_a_tbl,
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

  color_palette = c('#3A5FCD', '#FAFAFA','#CD3700')
  # color_palette = c('#3A5FCD', '#FFFFFF','#CD3700')
  # color_palette = c('#FFFFFF','#CD3700')

  p = plotly::plot_ly(plot_df,
                      x = ~umap3D_1,
                      y = ~umap3D_2,
                      z = ~umap3D_3,
                      type = 'scatter3d',
                      mode="markers",
                      alpha = I(alpha),
                      color = ~delta_log_abund,
                      colors = color_palette,
                      size = I(cell_size),
                      range_color = fc_limits)

  return(p)



}


