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
                          plot_edges = c("none", "all", "directed", "undirected"),
                          edge_significance = c("both", "one-sided"),
                          keep_colors = FALSE,
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

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(edge_significance) == "", TRUE, TRUE),
             error = function(e) FALSE),
    msg = paste('Argument plot_edges must be one of "both" or
                "one-sided".'))
  edge_significance = match.arg(edge_significance)

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
  umap_centers_delta_abund = dplyr::left_join(umap_centers_delta_abund, cond_b_vs_a_tbl, by=c("cell_group"="cell_group"))
  umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  umap_centers_delta_abund = umap_centers_delta_abund %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))


  corr_edge_coords_umap_delta_abund = hooke:::collect_pln_graph_edges(ccm,
                                                              umap_centers_delta_abund,
                                                              log_abundance_thresh,
                                                              model_for_pcors=model_for_pcors)

  if (edge_significance == "one-sided") {
    corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund %>%
      mutate(edge_type = case_when(
        (from_delta_q_value < q_value_thresh | to_delta_q_value < q_value_thresh ) ~ edge_type,
        TRUE ~ "hidden"
      ))
  } else {
    corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund %>%
      mutate(edge_type = case_when(
        (from_delta_q_value < q_value_thresh & to_delta_q_value < q_value_thresh ) ~ edge_type,
        TRUE ~ "hidden"
      ))
  }



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


  directed_cells_to_keep = unique(union(directed_edge_df$to, directed_edge_df$from))
  undirected_cells_to_keep = unique(union(undirected_edge_df$to, undirected_edge_df$from))
  others = cond_b_vs_a_tbl %>% filter(delta_q_value <= q_value_thresh) %>% pull(cell_group)
  cell_groups_to_keep = unique(union(directed_cells_to_keep, undirected_cells_to_keep))
  cell_groups_to_keep = unique(union(cell_groups_to_keep, others))


  if (keep_colors == FALSE) {

    cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund =
                                                        ifelse(cell_group %in% cell_groups_to_keep, delta_log_abund, 0))

    # cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund =
    #                                                       ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

    # corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund %>%
    #   dplyr::mutate(from_delta_log_abund = ifelse(from_delta_q_value <= q_value_thresh, from_delta_log_abund, 0)) %>%
    #   dplyr::mutate(to_delta_log_abund = ifelse(to_delta_q_value <= q_value_thresh, to_delta_log_abund, 0))
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
      label_df = label_df %>% filter(delta_q_value < q_value_thresh)

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




