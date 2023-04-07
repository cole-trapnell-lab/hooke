#' Plot a UMAP colored by how cells shift in a given contrast
#'
#' @param ccm A cell_count_model object.
#' @param cond_b_vs_a_tbl data.frame A data frame from compare_abundances.
#' @param log_abundance_thresh numeric Select cell groups by log abundance.
#' @param scale_shifts_by string A scale directed graph edges by "sender",
#'    "receiver", or "none".
#' @param edge_size numeric The size of edges in the plot.
#' @param cell_size numeric The size of cells in the plot.
#' @param q_value_thresh numeric Remove contrasts whose change in
#'    q-value exceeds q_value_thresh.
#' @param group_label_size numeric The size of group labels in the plot.
#' @param plot_labels string Choose cell groups to label.
#' @param fc_limits vector The range of cell abundance changes to
#'    include in the plot.
#' @param sender_cell_groups list Sender cell groups of directed graph.
#' @param receiver_cell_groups list Receiver cell groups of directed graph.
#' @param plot_edges string Type of edges to plot.
#' @param label_cell_groups list The cell_group labels to include in the plot.
#' @param repel_labels logical Repel overlapping plot labels.
#' @param model_for_pcors string The model to use for orienting graph
#'    edges from the PLNnetwork.
#' @param switch_label string The name of the cell_data_set column with cell_group identifiers.
#' @param sub_cds string A cell_data_set.
#' @param alpha numeric A the ggplot opacity. A value between 0 and 1.
#' @param x numeric The column number for the UMAP x coordinate.
#' @param y numeric The column number for the UMAP y coordinate.
#' @return A ggplot2 plot object.
#' @import ggplot2
#' @import dplyr
#' @export
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
                          model_for_pcors=c("reduced", "full"),
                          switch_label = NULL,
                          sub_cds = NULL,
                          alpha = 1.0,
                          x=1,
                          y=2){

  assertthat::assert_that(is(ccm, 'cell_count_model'))
  assertthat::assert_that(is.data.frame(cond_b_vs_a_tbl))
  assertthat::assert_that(is.numeric(log_abundance_thresh))
  assertthat::assert_that(is.numeric(edge_size) && edge_size > 0.0)
  assertthat::assert_that(is.numeric(cell_size) && cell_size > 0.0)
  assertthat::assert_that(is.numeric(q_value_thresh) && q_value_thresh> 0.0)
  assertthat::assert_that(is.numeric(group_label_size) && group_label_size > 0.0)
  assertthat::assert_that(is.vector(fc_limits) && length(fc_limits) == 2 && fc_limits[[1]] < fc_limits[[2]])
  assertthat::assert_that(is.list(label_cell_groups))
  assertthat::assert_that(is.logical(repel_labels))
  assertthat::assert_that(is.numeric(alpha) && alpha >= 0.0 && alpha <= 1.0)
  assertthat::assert_that(is.numeric(x) && x >= 1)
  assertthat::assert_that(is.numeric(y) && y >= 1)

  # Check that sender_cell_group and receiver_cell_group names are in cell_group names.
  assertthat::assert_that(is.null(sender_cell_groups) || is.vector(sender_cell_groups))
  assertthat::assert_that(is.null(receiver_cell_groups) || is.vector(receiver_cell_groups))

  assertthat::assert_that(is.null(sub_cds) || is(sub_cds, 'cell_data_set'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(model_for_pcors) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument model_for_pcors must be one of "reduced",',
                ' or "full".'))
  model_for_pcors = match.arg(model_for_pcors)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(scale_shifts_by) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument scale_shifts_by must be one of "receiver",',
                '"sender", or "none".'))
  scale_shifts_by = match.arg(scale_shifts_by)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(plot_labels) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument plot_labels must be one of "significant",',
                '"all", or "none".'))
  plot_labels = match.arg(plot_labels)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(plot_edges) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument plot_edges must be one of "all",',
                '"directed", "undirected", or "none".'))
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
    guides(color=guide_colourbar(title="log(\u0394 Abundance)")) + # unicode for captital delta
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
#' @noRd
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


#' @noRd
#' plot a cds in the same style as plot_contras
my_plot_cells <- function(cds,
                          color_cells_by = NULL,
                          cell_size=1,
                          legend_position="none",
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

  plot_df = as.data.frame(colData(cds))

  plot_df$cell = row.names(plot_df)
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[plot_df$cell,x]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[plot_df$cell,y]

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

    coord_matrix = as.data.frame(reducedDim(cds, "UMAP"))
    grp_assign = cds@colData %>% as.data.frame %>% dplyr::select(all_of(color_cells_by))
    coord_matrix = cbind(grp_assign, coord_matrix[row.names(grp_assign), ])
    colnames(coord_matrix)[1] = "cell_group"
    centroid_coords = aggregate(. ~ cell_group, data = coord_matrix, FUN = mean)

    if (dim(coord_matrix)[2] ==3) {
      colnames(centroid_coords) = c(color_cells_by, "umap2D_1", "umap2D_2")
    } else {
      colnames(centroid_coords) = c(color_cells_by, "umap2D_1", "umap2D_2", "umap2D_3")
    }

    if (repel_labels) {
      gp = gp + ggrepel::geom_label_repel(data = centroid_coords,
                                          mapping = aes(get(paste0("umap2D_", x)),
                                                        get(paste0("umap2D_", y)),
                                                        label=get(color_cells_by)),
                                          size=I(group_label_size),
                                          fill = "white")
    } else {
      gp = gp + geom_text(data = centroid_coords,
                          mapping = aes(get(paste0("umap2D_", x)),
                                        get(paste0("umap2D_", y)),
                                        label=get(color_cells_by)),
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
#' @noRd
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


#' @noRd
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
#' @noRd
plot_map <- function(data, edges, color_nodes_by = "", arrow.gap = 0.02, scale = F) {

  if (class(data) == "cell_count_set") {
    ccs = data
    plot_df = as.data.frame(ccs@cds_coldata)
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

#' @noRd
is_concordant <- function(cond_b_v_a_tbl, state_transition_graph) {

  parent_child_foldchanges(state_transition_graph, cond_b_v_a_tbl) %>%
    mutate(concordant = case_when( fold_changes == "parent increase, descendents increase" ~ TRUE,
                                   fold_changes == "parent decrease, descendents decrease" ~ TRUE,
                                   TRUE~ FALSE)) %>% ungroup #%>%
  # group_by(consitent) %>% tally() %>%
  # mutate(pct = n/sum(n))

}

#' @noRd
is_discordant <- function(cond_b_v_a_tbl, state_transition_graph) {

  parent_child_foldchanges(state_transition_graph, cond_b_v_a_tbl) %>%
    mutate(discordant = case_when( fold_changes == "parent decrease, descendents no change" ~ TRUE,
                                   fold_changes == "parent decrease, descendents increase" ~ TRUE,
                                   fold_changes == "parent increase, descendents no change" ~ TRUE,
                                   fold_changes == "parent increase, descendents decrease" ~ TRUE,
                                   TRUE~ FALSE)) %>% ungroup

}

#' @noRd
collect_psg_node_metadata <- function(ccm,
                                      color_nodes_by,
                                      label_nodes_by,
                                      group_nodes_by)
{
  cell_groups = ccm@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group) %>% unique()
  node_metadata = tibble(id=cell_groups)

  metadata_cols = c(color_nodes_by,
                    group_nodes_by)
  if (is.null(label_nodes_by) == FALSE && label_nodes_by != "cell_group")
    metadata_cols = c(metadata_cols, label_nodes_by)

  #G = edges %>% select(from, to, n, scaled_weight, distance_from_root)  %>% igraph::graph_from_data_frame(directed = T)
  cell_group_metadata = ccm@ccs@cds_coldata[,metadata_cols, drop=F] %>%
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
    label_by_metadata = cell_group_metadata[,c("cell_group", label_nodes_by), drop=F]
    colnames(label_by_metadata) = c("cell_group", "label_nodes_by")
    label_by_metadata = label_by_metadata %>%
      as.data.frame %>%
      count(cell_group, label_nodes_by) %>%
      group_by(cell_group) %>% slice_max(n, with_ties=FALSE) %>% dplyr::select(-n)
    node_metadata = left_join(node_metadata, label_by_metadata, by=c("id"="cell_group"))
  }else{
    node_metadata$label_nodes_by = node_metadata$id
  }
  node_metadata = node_metadata %>% distinct() %>% as.data.frame(stringsAsFactor=FALSE)
  row.names(node_metadata) = node_metadata$id

  return(node_metadata)
}

#' @noRd
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

  # make_subgraphs_for_groups <- function(grouping_set, G_nel, node_meta_df){
  #   nodes = node_meta_df %>% filter(group_nodes_by == grouping_set) %>% pull(id) %>% as.character %>% unique()
  #   #sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel))
  #   sg = graph::subGraph(snodes=nodes, graph=G_nel)
  #   return (sg)
  # }

  make_subgraphs_for_groups <- function(subgraph_ids, G_nel){
    nodes = subgraph_ids %>% pull(id) %>% as.character %>% unique()
    sg = list(graph=graph::subGraph(snodes=nodes, graph=G_nel), cluster=TRUE)
    #sg = graph::subGraph(snodes=nodes, graph=G_nel)
    return (sg)
  }

  subgraph_df = node_metadata %>%
    select(group_nodes_by, id) %>%
    group_by(group_nodes_by) %>%
    tidyr::nest(subgraph_ids = id) %>%
    summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
                                    .x = subgraph_ids,
                                    G_nel))
  subgraphs = subgraph_df$subgraph
  names(subgraphs) = subgraph_df$group_nodes_by

    #summarize(subgraph = graph::subGraph(id,graph=G_nel)
    # summarize(subgraph = purrr::map(.f = purrr::possibly(make_subgraphs_for_groups, NULL),
    #                                 .x = group_nodes_by,
    #                                 G_nel,
    #                                 node_metadata))

  if (is.null(edge_labels)== FALSE) {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraphs, edgeAttrs=list(label=edge_labels), recipEdges="distinct")
    label_df = data.frame("x" = gvizl@renderInfo@edges$labelX, "y" = gvizl@renderInfo@edges$labelY) %>%
      tibble::rownames_to_column("edge_name") %>%
      left_join(tibble("edge_name" = names(edge_labels), label=edge_labels))

  } else {
    gvizl = Rgraphviz::layoutGraph(G_nel, layoutType="dot", subGList=subgraphs, recipEdges="distinct")
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

#' Plot the state graph with annotations.
#'
#' @param ccm cell_count_model A cell_count_model.
#' @param state_graph
#' @param color_nodes_by
#' @param label_nodes_by
#' @param group_nodes_by
#' @param label_edges_by
#' @param edge_weights
#' @param arrow.gap
#' @param arrow_unit
#' @param bar_unit
#' @param node_size
#' @param min_edge_size
#' @param max_edge_size
#' @param unlabeled_groups
#' @param label_groups
#' @param hide_unlinked_nodes
#' @param group_label_font_size
#' @param edge_label_font_size
#' @param label_conn_linetype
#' @param legend_position
#' @param con_colour
#' @param group_outline
#' @return A ggplot2 plot object.
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
                                         min_edge_size=0.1,
                                         max_edge_size=2,
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
    edges$weight = 1
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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
                                     unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
                                     edge_thickness = replace_na(edge_thickness, min_edge_size))
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge = FALSE
  }

  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)

  p <- ggplot() +
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness, linetype=unsupported_edge), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

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
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
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
                                        expand = unit(2, "mm"),
                                        label.buffer=unit(1, "mm"),
                                        radius = unit(1.5, "mm"),
                                        label.margin = margin(1, 1, 1, 1, "mm"),
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
    label_df = layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y=mean(y))
    label_df$label = edge_labels[label_df$edge_name]
    #label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    #p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p = p + ggrepel::geom_text_repel(data = label_df,
                                mapping = aes(x, y, label=label),
                                size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
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
                                    loss_time = c("largest_loss", "largest_loss_time", "earliest_time", "latest_time", "peak_loss_time", "delta_log_abund_at_peak"),
                                    label_nodes_by=NULL,
                                    group_nodes_by=NULL,
                                    label_edges_by=NULL,
                                    edge_weights=NULL,
                                    arrow.gap=0.03,
                                    arrow_unit = 2,
                                    bar_unit = .075,
                                    node_size = 2,
                                    min_edge_size=0.1,
                                    max_edge_size=2,
                                    unlabeled_groups = c("Unknown"),
                                    label_subset = NULL,
                                    label_groups=TRUE,
                                    hide_unlinked_nodes=TRUE,
                                    group_label_font_size=6,
                                    edge_label_font_size=2,
                                    label_conn_linetype="dotted",
                                    legend_position = "none",
                                    con_colour = "darkgrey",
                                    group_outline=FALSE,
                                    control_ccm=perturbation_ccm,
                                    control_start_time=NULL,
                                    control_stop_time=NULL,
                                    ...)
{
  loss_time = match.arg(loss_time)

  if (is.null(control_start_time))
    control_start_time = start_time
  if (is.null(control_stop_time))
    control_stop_time = stop_time

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
                                                   control_ccm=control_ccm,
                                                   control_start_time=control_start_time,
                                                   control_stop_time=control_stop_time,
                                                   ...)

  #earliest_loss_tbl = earliest_loss_tbl %>% mutate(fill_alpha = ifelse(peak_time_in_ctrl_within_perturb_time_range, 1.0, 0.3))
  #print (earliest_loss_tbl)
  node_metadata = node_metadata %>% left_join(earliest_loss_tbl, by=c("id" ="cell_group"))

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

  layout_info = hooke:::layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
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
    ggplot2::geom_path(aes(x, y, group=edge_name, size=edge_thickness), data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

  if (is.null(label_subset)) {
    label_subset = unique(node_metadata$group_nodes_by)
  }

  label_subset = label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE){


    if (group_outline) {
      p = p + ggforce::geom_mark_rect(aes(x, y,
                                          col = group_nodes_by,
                                          label = group_nodes_by,
                                          filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                      size=0.5,
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      data=g)
    } else {
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
    }

  }


  # if numerical

  min_delta_log_abund = -3
  max_delta_log_abund = 3
  g = g %>% mutate(fill_val = !!sym(loss_time))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")){
    g = g %>% mutate(fill_val = ifelse(fill_val < min_delta_log_abund, min_delta_log_abund, fill_val),
                     fill_val = ifelse(fill_val > max_delta_log_abund, max_delta_log_abund, fill_val),
                     fill_val = ifelse(peak_time_in_ctrl_within_perturb_time_range, fill_val, NA),
                     fill_border_color = ifelse(peak_time_in_ctrl_within_perturb_time_range == FALSE,
                                                "out-of-range",
                                                ifelse(is_lost_at_peak,
                                                  "significant",
                                                  "notsignificant")
                                                ))
  }

  p = p + ggnewscale::new_scale_fill() +
    ggnetwork::geom_nodelabel(data = g,
                              aes(x, y,
                                  fill = fill_val,
                                  color = fill_border_color,
                                  label = label_nodes_by),
                              size = node_size)
  p = p + scale_color_manual(values=c("out-of-range"="black", "significant"="black", "notsignificant"="lightgrey"))

  if (loss_time %in% c("largest_loss", "delta_log_abund_at_peak")){
    p = p + scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3",
                                 limits=c(min_delta_log_abund - abs(min_delta_log_abund)*0.05,
                                          max_delta_log_abund + abs(max_delta_log_abund)*0.05))
  }else{
    p = p + scale_fill_stepsn(n.breaks=5, colours = terrain.colors(5)) #+ scale_fill_gradient2(low = "royalblue3", mid = "white", high="orangered3")
  }

  p = p + guides(color = "none")

  #p = p + guides(fill = guide_legend(position=legend_position))

  if (is.null(edge_labels) == FALSE) {
    label_df = layout_info$bezier_df %>%
      group_by(edge_name) %>%
      summarize(x = mean(x), y=mean(y))
    label_df$label = edge_labels[label_df$edge_name]
    #label_df = layout_info$label_df
    #p = p +  ggnetwork::geom_nodetext(data = label_df,
    #                                  aes(x,y, label = label))
    #p = p + geom_text(data = label_df,
    #                  aes(x,y, label = label),
    #                  size=edge_label_font_size)
    p = p + ggrepel::geom_text_repel(data = label_df,
                                     mapping = aes(x, y, label=label),
                                     size=edge_label_font_size)
  }

  p = p + scale_size_identity() +
    ggnetwork::theme_blank() +
    hooke_theme_opts() #+
    #theme(legend.position=legend_position)
  return(p)
}


#' @noRd
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

#' @noRd
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
                                               label_edges_by=NULL,
                                               edge_weights=NULL,
                                               fc_limits=c(-3,3),
                                               arrow.gap=0.03,
                                               arrow_unit = 2,
                                               bar_unit = .075,
                                               node_size = 6,
                                               min_edge_size=0.1,
                                               max_edge_size=2,
                                               unlabeled_groups = c("Unknown"),
                                               label_subset = NULL,
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

  node_metadata = collect_psg_node_metadata(ccm, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)

  comp_abund_table[["contrast"]] = comp_abund_table[[contrast]]
  comp_abund_table$contrast = as.factor(comp_abund_table$contrast)

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

  abund_fc_df = comp_abund_table

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
  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')

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
                                      expand = unit(2, "mm"),
                                      label.buffer=unit(1, "mm"),
                                      radius = unit(1.5, "mm"),
                                      label.margin = margin(1, 1, 1, 1, "mm"),
                                      label.fontsize=group_label_font_size,
                                      label.fontface="plain",
                                      con.linetype=label_conn_linetype,
                                      con.colour=con_colour,
                                      show.legend = F) +
        ggforce::geom_mark_rect(aes(x, y,
                                    fill = contrast,
                                    color = group_nodes_by,
                                    label = group_nodes_by,
                                    filter = group_nodes_by %in% label_subset),
                                size=0.5,
                                expand = unit(2, "mm"),
                                label.buffer=unit(1, "mm"),
                                radius = unit(1.5, "mm"),
                                label.margin = margin(1, 1, 1, 1, "mm"),
                                label.fontsize=group_label_font_size,
                                label.fontface="plain",
                                con.linetype=label_conn_linetype,
                                con.colour=con_colour,
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
                                expand = unit(2, "mm"),
                                label.buffer=unit(1, "mm"),
                                radius = unit(1.5, "mm"),
                                label.margin = margin(1, 1, 1, 1, "mm"),
                                label.fontsize=group_label_font_size,
                                label.fontface="plain",
                                con.linetype=label_conn_linetype,
                                con.colour=con_colour,
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

  #p = p +
  #  ggplot2::geom_path(aes(x, y, group=edge_name), data=bezier_df, arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"))

  p = p + facet_wrap(~contrast)

  p = p + scale_size(range=c(1, node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
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
                                             fract_expr=0.0,
                                             mean_expr=0.0,
                                             scale_to_range = FALSE,
                                             color_nodes_by=NULL,
                                             label_nodes_by=NULL,
                                             group_nodes_by=NULL,
                                             label_edges_by=NULL,
                                             edge_weights=NULL,
                                             arrow.gap=0.03,
                                             arrow_unit = 2,
                                             bar_unit = .075,
                                             min_node_size = 0.25,
                                             max_node_size = 2,
                                             min_edge_size=0.1,
                                             max_edge_size=2,
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

  node_metadata = collect_psg_node_metadata(ccm, color_nodes_by=NULL, label_nodes_by, group_nodes_by)

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

  layout_info = layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  if (is.null(edge_weights) == FALSE){
    bezier_df = left_join(bezier_df, edges)
    bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                     edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
                                     unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
                                     edge_thickness = replace_na(edge_thickness, min_edge_size))
  }else{
    bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge = FALSE
  }

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
      max_expr = max(mean_expression),
      fraction_max = ifelse (max_expr > 0, mean_expression / max_expr, 0),
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
  p = p +
    ggplot2::geom_path(aes(x, y, group=edge_name), colour=con_colour, data=bezier_df %>% distinct(), arrow = arrow(angle=30, length = unit(arrow_unit, "pt"), type="closed"), linejoin='mitre')


  # if (is.null(label_subset)) {
  #   label_subset = unique(node_metadata$group_nodes_by)
  # }
  #
  # label_subset = label_subset[label_subset != unlabeled_groups]

  if (is.null(group_nodes_by) == FALSE){
    p = p + ggforce::geom_mark_rect(aes(x, y,
                                        fill = gene_short_name,
                                        color = group_nodes_by,
                                        filter = group_nodes_by %in% unlabeled_groups == FALSE),
                                    size=0.5,
                                    expand = unit(2, "mm"),
                                    label.buffer=unit(1, "mm"),
                                    radius = unit(1.5, "mm"),
                                    label.margin = margin(1, 1, 1, 1, "mm"),
                                    label.fontsize=group_label_font_size,
                                    label.fontface="plain",
                                    con.linetype=label_conn_linetype,
                                    con.colour=con_colour,
                                    show.legend = F) #+
      # ggforce::geom_mark_rect(aes(x, y,
      #                             fill = gene_short_name,
      #                             color = group_nodes_by,
      #                             label = group_nodes_by,
      #                             filter = group_nodes_by %in% label_subset),
      #                         size=0.5,
      #                         expand = unit(2, "mm"),
      #                         label.buffer=unit(1, "mm"),
      #                         radius = unit(1.5, "mm"),
      #                         label.margin = margin(1, 1, 1, 1, "mm"),
      #                         label.fontsize=group_label_font_size,
      #                         label.fontface="plain",
      #                         con.linetype=label_conn_linetype,
      #                         con.colour=con_colour,
      #                         show.legend = F)

    p = p + scale_fill_manual(values=rep("white", length(unique(g$gene_short_name))))

  }
  p = p + guides(fill = "none")
  p = p + ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr),
                          aes(x, y,
                              size = fraction_max,
                              color = I(con_colour))) +
    ggnewscale::new_scale_color() +
    ggnetwork::geom_nodes(data = g %>% filter(gene_expr & fraction_max > 0),
                          aes(x, y,
                              size = fraction_max,
                              color = get(color_nodes_by))) +
    labs(color = expression_legend_label) +
    viridis::scale_color_viridis(option = "viridis")

  if (is.null(edge_labels) == FALSE) {
    p = p +  ggnetwork::geom_nodetext(data = label_df,
                                      aes(x,y, label = label), size=3)
  }

  p = p + facet_wrap(~gene_short_name)

  p = p + scale_size_continuous(labels = scales::percent, range=c(min_node_size, max_node_size)) +
    ggnetwork::theme_blank() +
    hooke_theme_opts() +
    theme(legend.position=legend_position) + guides(fill = "none")
  return(p)
}


# '
#' @param ccs
#' @param x_col
#' @param normalize
#' @param plot_zeroes
#' @param plot_points
#' @param log_scale
#' @param nrow
#' @param legend_position
#' @noRd
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
    tibble::rownames_to_column("cell_group") %>%
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


#' @noRd
plot_cells_highlight = function(ccs, group_to_highlight, colname) {

  plot_df = as.data.frame(ccs@cds_coldata)
  plot_df$cell_group = ccs@cds_coldata[[colname]]

  plot_df$cell = row.names(plot_df)
  plot_df$umap2D_1 <- ccs@cds_reduced_dims[["UMAP"]][plot_df$cell,x]
  plot_df$umap2D_2 <- ccs@cds_reduced_dims[["UMAP"]][plot_df$cell,y]

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


#' @noRd
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

  plot_df$umap3D_1 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell,1]
  plot_df$umap3D_2 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell,2]
  plot_df$umap3D_3 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell,3]

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

#' Default plotting options for ggplot2
#'
#' return A ggplot2 theme object.
#' @export
hooke_theme_opts <- function(){
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    #theme(axis.line.x = element_line(size=0.25, color="black")) +
    #theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}





