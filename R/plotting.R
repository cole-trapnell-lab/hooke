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
                      path_df = path_df,
                      edge_size = 1,
                      path_color = "black",
                      color_path_by = NULL,
                      cond_b_vs_a_tbl = NULL,
                      directed = TRUE,
                      x = 1,
                      y = 2,
                      ...) {

  gp = my_plot_cells(data,  x=x, y=y, ...)

  if (class(data) == "cell_count_set") {
    umap_centers = centroids(data)
  } else if (class(data) == "cell_count_model") {
    umap_centers = centroids(data@ccs)
  }

  path_df = add_umap_coords(path_df, umap_centers)

  if (directed) {
    if (is.null(color_path_by) == FALSE) {

      gp = gp +
        ggnewscale::new_scale_color() +
        geom_segment(data = path_df,
                     aes(x = get(paste0("umap_to_", x)),
                         y = get(paste0("umap_to_", y)),
                         xend=get(paste0("from_to_", x)),
                         yend = get(paste0("from_to_", y)),
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
                         xend=get(paste0("from_to_", x)),
                         yend = get(paste0("from_to_", y)),
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

#' Plot a graph that summarized an assembled state transition graph
#'
#' @export
plot_state_transition_graph <- function(ccm,
                                        edges,
                                        color_nodes_by=NULL,
                                        label_nodes_by=NULL,
                                        group_nodes_by=NULL,
                                        layer_nodes_by=NULL,
                                        cond_b_v_a_tbl = NULL,
                                        fc_limits = NULL,
                                        genes = list(),
                                        qval_threshold = 0.05,
                                        arrow.gap=0.03,
                                        arrow_unit = 2,
                                        bar_unit = .075,
                                        node_size = 2,
                                        num_layers=10,
                                        edge_size=0.5,
                                        fract_expr = 0.1,
                                        unlabeled_groups = c("Unknown"),
                                        hide_unlinked_nodes=TRUE,
                                        label_font_size=6,
                                        label_conn_linetype="dotted",
                                        legend_position = "none",
                                        group_outline=FALSE){

  #edges = hooke:::distance_to_root(edges)
  edges = edges %>% dplyr::ungroup()

  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = colData(ccm@ccs@cds)[,c(color_nodes_by,
                                                label_nodes_by,
                                                group_nodes_by,
                                                layer_nodes_by)] %>%
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
    node_metadata = node_metadata %>% left_join(sub_gene_expr, by = c("id" = "cell_group"))

    node_metadata = node_metadata %>%
      mutate(gene_expr = case_when(
        fraction_expressing >= fract_expr ~ TRUE,
        TRUE ~ FALSE))

    color_nodes_by = "mean_expression"
    group_outline = TRUE

  }
  # else {
  #   # if not specificed, ignore this information
  #   node_metadata = node_metadata %>% mutate(fraction_expressing = 1.0,
  #                                            mean_expression = 1.0,
  #                                            gene_expr = TRUE)
  # }


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

  gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraph_df$subgraph)
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
  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)


  # add level information
  #g = g %>% left_join(level_df %>% rownames_to_column("id"), by = c("vertex.names"="id"))
  #g = g %>% left_join(regulator_score_df, by = c("vertex.names" = "gene_id") )


  # p <- ggplot(mapping = aes(x, y, xend = xend, yend = yend)) +
  #   # draw activator edges
  #   ggnetwork::geom_edges(data = g,
  #                         arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))
  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))
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
      p = p + ggnewscale::new_scale_fill() +
        ggnetwork::geom_nodelabel(data = g,
                                  aes(x, y,
                                      fill = !!sym(color_nodes_by),
                                      label = label_nodes_by),
                                      size = node_size) +
        labs(fill = color_nodes_by)

      if (color_nodes_by == "mean_expression") {
        p = p + viridis::scale_fill_viridis(option = "viridis")
      } else {
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
