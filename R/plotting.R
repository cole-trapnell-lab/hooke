#' Plot a UMAP colored by how cells shift in a given contrast
#'
#' @param ccs A cell_count_set object.
#' @param cond_b_vs_a_tbl data.frame A data frame from compare_abundances.
#' @param mask a list of cell types to gray out in the plots
#' @param cell_size numeric The size of cells in the plot.
#' @param q_value_thresh numeric Remove contrasts whose change in
#'    q-value exceeds q_value_thresh.
#' @param group_label_size numeric The size of group labels in the plot.
#' @param plot_labels string Choose cell groups to label.
#' @param fc_limits vector The range of cell abundance changes to
#'    include in the plot.
#' @param downsample how much to downsample the plots
#' @param x numeric The column number for the UMAP x coordinate.
#' @param y numeric The column number for the UMAP y coordinate.
#' @export
plot_abundance <- function(ccs,
                           cond_b_vs_a_tbl,
                           mask = list(),
                           cell_size = 1,
                           q_value_thresh = 1.0,
                           fc_limits = c(-3, 3),
                           plot_labels = c("significant", "all", "none"),
                           alpha = 1.0,
                           group_label_size = 2,
                           x = 1,
                           y = 2,
                           downsample = NULL) {
  plot_labels <- match.arg(plot_labels)

  cond_b_vs_a_tbl <- cond_b_vs_a_tbl %>% dplyr::mutate(
    delta_log_abund =
      ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0)
  )




  plot_df <- ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell <- row.names(plot_df)

  plot_df <- plot_df[rownames(reducedDim(ccs@cds, type = "UMAP")), ]

  plot_df$umap2D_1 <- reducedDim(ccs@cds, type = "UMAP")[plot_df$cell, x]
  plot_df$umap2D_2 <- reducedDim(ccs@cds, type = "UMAP")[plot_df$cell, y]

  plot_df <- dplyr::left_join(plot_df,
    cond_b_vs_a_tbl,
    by = c("cell_group" = "cell_group")
  )

  if (is.null(fc_limits)) {
    fc_limits <- range(plot_df$delta_log_abund)
  } else {
    min <- fc_limits[1]
    max <- fc_limits[2]
    plot_df <- plot_df %>%
      mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
      mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  }

  # option to downsample for easier plotting
  if (is.null(downsample) == FALSE) {
    n <- min(nrow(plot_df), downsample)
    plot_df <- plot_df[sample(nrow(plot_df), n), ]
  }

  if (length(mask) > 0) {
    plot_df$mask <- ifelse(plot_df$cell_group %in% mask, T, F)
  } else {
    plot_df$mask <- F
  }


  gp <- ggplot() +
    geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2),
      color = "black",
      size = 1.5 * cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>% filter(mask == T),
      aes(umap2D_1, umap2D_2),
      color = "light gray",
      size = cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>% filter(mask == F),
      aes(umap2D_1, umap2D_2),
      color = "white",
      size = cell_size,
      stroke = 0
    ) +
    geom_point(
      data = plot_df %>%
        filter(mask == F) %>%
        arrange(
          !is.na(abs(delta_log_abund)),
          abs(delta_log_abund)
        ),
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
    ) +
    guides(color = guide_colourbar(title = "log(\u0394 Abundance)")) + #
    monocle3:::monocle_theme_opts()






  if (plot_labels != "none") {
    label_df <- centroids(ccs)
    label_df <- left_join(label_df, cond_b_vs_a_tbl, by = "cell_group")

    if (plot_labels == "significant") {
      # sig_cell_groups = cond_b_vs_a_tbl %>% filter(delta_q_value < q_value_thresh) %>% pull(cell_group)
      # label_df = label_df %>% filter(cell_group %in% sig_cell_groups)
      label_df <- label_df %>% filter(delta_q_value < q_value_thresh)
    }

    gp <- gp + ggrepel::geom_label_repel(
      data = label_df,
      mapping = aes(get(paste0("umap_", x)),
        get(paste0("umap_", y)),
        label = cell_group
      ),
      size = I(group_label_size),
      fill = "white"
    )
  }

  gp <- gp + xlab(paste0("UMAP", x)) + ylab(paste0("UMAP", y))

  return(gp)
}




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
#' @param downsample how much to downsample the plots
#' @param x numeric The column number for the UMAP x coordinate.
#' @param y numeric The column number for the UMAP y coordinate.
#' @return A ggplot2 plot object.
#' @import ggplot2
#' @import dplyr
#' @export
plot_contrast <- function(ccm,
                          cond_b_vs_a_tbl,
                          log_abundance_thresh = -5,
                          scale_shifts_by = c("receiver", "sender", "none"),
                          edge_size = 2,
                          cell_size = 1,
                          q_value_thresh = 1.0,
                          group_label_size = 2,
                          plot_labels = c("significant", "all", "none"),
                          fc_limits = c(-3, 3),
                          sender_cell_groups = NULL,
                          receiver_cell_groups = NULL,
                          plot_edges = c("none", "all", "directed", "undirected"),
                          edge_significance = c("both", "one-sided"),
                          keep_colors = FALSE,
                          label_cell_groups = list(),
                          repel_labels = TRUE,
                          model_for_pcors = c("reduced", "full"),
                          switch_label = NULL,
                          sub_cds = NULL,
                          alpha = 1.0,
                          mask = list(),
                          downsample = NULL,
                          x = 1,
                          y = 2) {
  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(is.data.frame(cond_b_vs_a_tbl))
  assertthat::assert_that(is.numeric(log_abundance_thresh))
  assertthat::assert_that(is.numeric(edge_size) && edge_size > 0.0)
  assertthat::assert_that(is.numeric(cell_size) && cell_size > 0.0)
  assertthat::assert_that(is.numeric(q_value_thresh) && q_value_thresh > 0.0)
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

  assertthat::assert_that(is.null(sub_cds) || is(sub_cds, "cell_data_set"))

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(model_for_pcors) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste(
      'Argument model_for_pcors must be one of "reduced",',
      ' or "full".'
    )
  )
  model_for_pcors <- match.arg(model_for_pcors)

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(scale_shifts_by) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste(
      'Argument scale_shifts_by must be one of "receiver",',
      '"sender", or "none".'
    )
  )
  scale_shifts_by <- match.arg(scale_shifts_by)

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(plot_labels) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste(
      'Argument plot_labels must be one of "significant",',
      '"all", or "none".'
    )
  )
  plot_labels <- match.arg(plot_labels)

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(plot_edges) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste(
      'Argument plot_edges must be one of "all",',
      '"directed", "undirected", or "none".'
    )
  )
  plot_edges <- match.arg(plot_edges)

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(edge_significance) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste('Argument plot_edges must be one of "both" or
                "one-sided".')
  )
  edge_significance <- match.arg(edge_significance)

  if (!is.null(sub_cds)) {
    colData(ccm@ccs@cds)[["cell_group"]] <- colData(ccm@ccs@cds)[[ccm@ccs@info$cell_group]]
    # get cell groups in sub
    sub_cell_group <- colData(sub_cds)[[ccm@ccs@info$cell_group]] %>% unique()

    # subset to select group
    ccm@ccs@cds <- ccm@ccs@cds[, colnames(sub_cds)]
    # ccm@ccs@cds = ccm@ccs@cds[,colData(ccm@ccs@cds)[["cell_group"]] %in% sub_cell_group]
    ccm@ccs@metadata[["cell_group_assignments"]] <- ccm@ccs@metadata[["cell_group_assignments"]][colnames(ccm@ccs@cds), ]

    cond_b_vs_a_tbl <- cond_b_vs_a_tbl %>% filter(cell_group %in% sub_cell_group)

    # switch coords to the new ones
    sub_umap <- reducedDims(sub_cds)[["UMAP"]]
    reducedDims(ccm@ccs@cds)[["UMAP"]] <- sub_umap[colnames(ccm@ccs@cds), ]
  }

  umap_centers <- centroids(ccm@ccs)

  umap_centers_delta_abund <- umap_centers
  umap_centers_delta_abund <- dplyr::left_join(umap_centers_delta_abund, cond_b_vs_a_tbl, by = c("cell_group" = "cell_group"))
  umap_centers_delta_abund <- umap_centers_delta_abund %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  umap_centers_delta_abund <- umap_centers_delta_abund %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))


  corr_edge_coords_umap_delta_abund <- hooke:::collect_pln_graph_edges(ccm,
    umap_centers_delta_abund,
    log_abundance_thresh,
    model_for_pcors = model_for_pcors
  )

  if (edge_significance == "one-sided") {
    corr_edge_coords_umap_delta_abund <- corr_edge_coords_umap_delta_abund %>%
      mutate(edge_type = case_when(
        (from_delta_q_value < q_value_thresh | to_delta_q_value < q_value_thresh) ~ edge_type,
        TRUE ~ "hidden"
      ))
  } else {
    corr_edge_coords_umap_delta_abund <- corr_edge_coords_umap_delta_abund %>%
      mutate(edge_type = case_when(
        (from_delta_q_value < q_value_thresh & to_delta_q_value < q_value_thresh) ~ edge_type,
        TRUE ~ "hidden"
      ))
  }



  directed_edge_df <- corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("directed_to_from", "directed_from_to"))
  undirected_edge_df <- corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("undirected"))

  if (is.null(sender_cell_groups) == FALSE) {
    directed_edge_df <- directed_edge_df %>% dplyr::filter(from %in% sender_cell_groups)
    undirected_edge_df <- undirected_edge_df %>% dplyr::filter(from %in% sender_cell_groups | to %in% sender_cell_groups)
  }

  if (is.null(receiver_cell_groups) == FALSE) {
    directed_edge_df <- directed_edge_df %>% dplyr::filter(to %in% receiver_cell_groups)
    undirected_edge_df <- undirected_edge_df %>% dplyr::filter(from %in% receiver_cell_groups | to %in% receiver_cell_groups)
  }

  # corr_edge_coords_umap_delta_abund = left_join(corr_edge_coords_umap_delta_abund,
  #                                               umap_centers,
  #                                               by=c("from"="cell_group"))

  if (scale_shifts_by == "sender") {
    directed_edge_df <- directed_edge_df %>%
      dplyr::group_by(to) %>%
      dplyr::mutate(
        flow_factor = -pmin(0, pcor),
        total_weight = sum(flow_factor),
        scaled_weight = flow_factor / total_weight
      )
    undirected_edge_df <- undirected_edge_df %>%
      dplyr::mutate(scaled_weight = abs(pcor) / max(abs(pcor)))
  } else if (scale_shifts_by == "receiver") {
    directed_edge_df <- directed_edge_df %>%
      dplyr::group_by(from) %>%
      dplyr::mutate(
        flow_factor = -pmin(0, pcor),
        total_weight = sum(flow_factor),
        scaled_weight = flow_factor / total_weight
      )
    undirected_edge_df <- undirected_edge_df %>%
      dplyr::mutate(scaled_weight = abs(pcor) / max(abs(pcor)))
  } else {
    directed_edge_df <- directed_edge_df %>%
      dplyr::mutate(scaled_weight = abs(pcor) / max(abs(pcor)))
    undirected_edge_df <- undirected_edge_df %>%
      dplyr::mutate(scaled_weight = abs(pcor) / max(abs(pcor)))
  }


  directed_cells_to_keep <- unique(union(directed_edge_df$to, directed_edge_df$from))
  undirected_cells_to_keep <- unique(union(undirected_edge_df$to, undirected_edge_df$from))
  others <- cond_b_vs_a_tbl %>%
    filter(delta_q_value <= q_value_thresh) %>%
    pull(cell_group)
  cell_groups_to_keep <- unique(union(directed_cells_to_keep, undirected_cells_to_keep))
  cell_groups_to_keep <- unique(union(cell_groups_to_keep, others))


  if (keep_colors == FALSE) {
    cond_b_vs_a_tbl <- cond_b_vs_a_tbl %>% dplyr::mutate(
      delta_log_abund =
        ifelse(cell_group %in% cell_groups_to_keep, delta_log_abund, 0)
    )

    # cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund =
    #                                                       ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))

    # corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund %>%
    #   dplyr::mutate(from_delta_log_abund = ifelse(from_delta_q_value <= q_value_thresh, from_delta_log_abund, 0)) %>%
    #   dplyr::mutate(to_delta_log_abund = ifelse(to_delta_q_value <= q_value_thresh, to_delta_log_abund, 0))
  }

  gp <- plot_abundance(ccm@ccs, cond_b_vs_a_tbl,
    cell_size = cell_size, mask = mask, q_value_thresh = q_value_thresh,
    fc_limits = fc_limits, alpha = alpha, x = x, y = y, downsample = downsample, group_label_size = group_label_size,
    plot_labels = "none"
  )

  # plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  # plot_df$cell = row.names(plot_df)
  #
  # plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,x]
  # plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,y]
  #
  # plot_df = dplyr::left_join(plot_df,
  #                            cond_b_vs_a_tbl,
  #                            # cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
  #                            by=c("cell_group"="cell_group"))
  #
  #
  # if (is.null(fc_limits)) {
  #   fc_limits = range(plot_df$delta_log_abund)
  # } else {
  #   min = fc_limits[1]
  #   max = fc_limits[2]
  #   plot_df = plot_df %>%
  #     mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
  #     mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  # }
  #
  # # option to downsample for easier plotting
  # if (is.null(downsample) == FALSE) {
  #   n = min(nrow(plot_df), downsample)
  #   plot_df = plot_df[sample(nrow(plot_df), n),]
  # }


  # plot
  # gp = ggplot() +
  #   geom_point(
  #     data = plot_df,
  #     aes(umap2D_1, umap2D_2),
  #     color = "black",
  #     size = 1.5 * cell_size,
  #     stroke = 0
  #   ) +
  #   geom_point(
  #     data = plot_df,
  #     aes(umap2D_1, umap2D_2),
  #     color = "white",
  #     size = cell_size,
  #     stroke = 0
  #   ) +
  #   geom_point(
  #     data = plot_df %>%
  #       arrange(!is.na(abs(delta_log_abund)),
  #               abs(delta_log_abund)),
  #     aes(umap2D_1, umap2D_2, color = delta_log_abund),
  #     size = cell_size,
  #     alpha = alpha,
  #     stroke = 0
  #   ) +
  #   scale_color_gradient2(
  #     low = "royalblue3",
  #     mid = "white",
  #     high = "orangered3",
  #     na.value = "white",
  #     limits = fc_limits
  #   )  +
  #   guides(color=guide_colourbar(title="log(\u0394 Abundance)")) + # unicode for captital delta
  #   #theme_void() +
  #   #theme(legend.position = "none") +
  #   monocle3:::monocle_theme_opts()

  # force to be empty
  if (plot_edges == "directed") {
    undirected_edge_df <- undirected_edge_df %>% dplyr::filter(!edge_type %in% c("undirected"))
  } else if (plot_edges == "undirected") {
    directed_edge_df <- directed_edge_df %>% dplyr::filter(edge_type %in% c("undirected"))
  }

  if (plot_edges != "none") {
    gp <- gp +
      geom_segment(
        data = undirected_edge_df,
        aes(
          x = get(paste0("to_umap_", x)),
          y = get(paste0("to_umap_", y)),
          xend = get(paste0("from_umap_", x)),
          yend = get(paste0("from_umap_", y)),
          size = edge_size * scaled_weight
        ),
        color = "lightgray"
      ) +
      geom_segment(
        data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
        aes(
          x = get(paste0("to_umap_", x)),
          y = get(paste0("to_umap_", y)),
          xend = get(paste0("from_umap_", x)),
          yend = get(paste0("from_umap_", y)),
          size = edge_size * scaled_weight
        ),
        color = "black"
      ) +
      geom_segment(
        data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
        aes(
          x = get(paste0("to_umap_", x)),
          y = get(paste0("to_umap_", y)),
          xend = (get(paste0("to_umap_", x)) + get(paste0("from_umap_", x))) / 2,
          yend = (get(paste0("to_umap_", y)) + get(paste0("from_umap_", x))) / 2,
          size = edge_size * scaled_weight
        ),
        color = "black",
        linejoin = "mitre",
        arrow = arrow(type = "closed", angle = 30, length = unit(1, "mm"))
      ) +
      geom_segment(
        data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
        aes(
          x = get(paste0("from_umap_", x)),
          y = get(paste0("from_umap_", y)),
          xend = get(paste0("to_umap_", x)),
          yend = get(paste0("to_umap_", y)),
          size = edge_size * scaled_weight
        ),
        color = "black"
      ) +
      geom_segment(
        data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
        aes(
          x = get(paste0("from_umap_", x)),
          y = get(paste0("from_umap_", y)),
          xend = (get(paste0("from_umap_", x)) + get(paste0("to_umap_", x))) / 2,
          yend = (get(paste0("from_umap_", y)) + get(paste0("to_umap_", y))) / 2,
          size = edge_size * scaled_weight
        ),
        color = "black",
        linejoin = "mitre",
        arrow = arrow(type = "closed", angle = 30, length = unit(1, "mm"))
      ) +
      scale_size_identity()
  }

  if (is.null(switch_label) == FALSE) {
    label_df <- centroids(ccs, switch_group = switch_label)

    if (length(label_cell_groups) > 0) {
      label_df <- label_df %>% filter(cell_group %in% label_cell_groups)
    }

    # can't have both
    plot_labels <- "none"

    if (repel_labels) {
      gp <- gp + ggrepel::geom_label_repel(
        data = label_df,
        mapping = aes(get(paste0("umap_", x)),
          get(paste0("umap_", y)),
          label = cell_group
        ),
        size = I(group_label_size),
        fill = "white"
      )
    } else {
      gp <- gp + geom_text(
        data = label_df,
        mapping = aes(get(paste0("umap_", x)),
          get(paste0("umap_", y)),
          label = cell_group
        ),
        size = I(group_label_size)
      )
    }
  }


  if (plot_labels != "none") {
    label_df <- umap_centers_delta_abund

    if (plot_labels == "significant") {
      label_df <- label_df %>% filter(delta_q_value < q_value_thresh)
    }

    if (length(label_cell_groups) > 0) {
      label_df <- label_df %>% filter(cell_group %in% label_cell_groups)
    }

    if (repel_labels) {
      gp <- gp + ggrepel::geom_label_repel(
        data = label_df,
        mapping = aes(get(paste0("umap_", x)),
          get(paste0("umap_", y)),
          label = cell_group
        ),
        size = I(group_label_size),
        fill = "white"
      )
    } else {
      gp <- gp + geom_text(
        data = label_df,
        mapping = aes(get(paste0("umap_", x)),
          get(paste0("umap_", y)),
          label = cell_group
        ),
        size = I(group_label_size)
      )
    }
  }

  gp <- gp + xlab(paste0("UMAP", x)) + ylab(paste0("UMAP", y))
  return(gp)
}

#' returns different color palettes
#' @param num_colors the number of colors needed
#' @param
#' @noRd
get_colors <- function(num_colors, type = "rainbow") {
  if (type == "rainbow") {
    colors <-
      c(
        "18h" = "#DF4828",
        "24h" = "#E78C35",
        "36h" = "#F6C141",
        "48h" = "#4EB265",
        "72h" = "#1965B0"
      )
  } else if (type == "vibrant") {
    colors <-
      c(
        "#EE7733", "#0077BB", "#228833", "#33BBEE", "#EE3377", "#CC3311",
        "#AA3377", "#009988", "#004488", "#DDAA33", "#99CC66", "#D590DD"
      )
  } else {
    colors <-
      c(
        "#4477AA",
        "#EE6677",
        "#228833",
        "#CCBB44",
        "#66CCEE",
        "#AA3377",
        "#BBBBBB"
      )
  }

  full_spectrum <- colorRampPalette(colors)(num_colors)

  return(full_spectrum)
}


#' @noRd
#' plot a cds in the same style as plot_contras
my_plot_cells <- function(cds,
                          color_cells_by = NULL,
                          cell_size = 1,
                          legend_position = "none",
                          q_value_thresh = 1.0,
                          fc_limits = c(-3, 3),
                          plot_labels = TRUE,
                          label_cells_by = NULL,
                          group_label_size = 2,
                          repel_labels = TRUE,
                          lab_title = NULL,
                          color_values = NULL,
                          alpha = 1.0,
                          x = 1,
                          y = 2) {
  plot_df <- as.data.frame(colData(cds))

  plot_df$cell <- row.names(plot_df)
  plot_df$umap2D_1 <- reducedDim(cds, type = "UMAP")[plot_df$cell, x]
  plot_df$umap2D_2 <- reducedDim(cds, type = "UMAP")[plot_df$cell, y]

  plot_df$color_cells_by <- plot_df[[color_cells_by]]

  gp <- ggplot() +
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
    num_colors <- unique(plot_df$timepoint) %>%
      sort() %>%
      length()
    full_spectrum_timepoint <- get_colors(num_colors, type = "rainbow")
    names(full_spectrum_timepoint) <- unique(plot_df$timepoint) %>% sort()

    gp <- gp +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
          color = as.character(timepoint)
        ),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      scale_color_manual(values = full_spectrum_timepoint)
  } else if (color_cells_by %in% c("coefficients", "delta_log_needed")) {
    gp <- gp +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
          color = color_cells_by
        ),
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
    gp <- gp +
      # geom_point(
      #   data = plot_df %>% filter(is.na(color_cells_by)),
      #   aes(umap2D_1, umap2D_2),
      #   color = "gray",
      #   size = cell_size,
      #   alpha = 0.2,
      #   stroke = 0) +
      geom_point(
        data = plot_df %>% filter(
          !is.na(color_cells_by),
          color_cells_by > 0
        ),
        aes(umap2D_1, umap2D_2,
          color = log_dose
        ),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      # viridis::scale_color_viridis(option = "C")
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = fc_limits
      )
  } else if (color_cells_by %in% c("viridis", "inferno", "C")) {
    gp <- gp +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
          color = coefficients
        ),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      viridis::scale_color_viridis(option = color_cells_by)
  } else if (color_cells_by == "delta_log_abund") {
    gp <- gp +
      geom_point(
        data = plot_df %>%
          arrange(
            !is.na(abs(delta_log_abund)),
            abs(delta_log_abund)
          ),
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
  } else if (is.numeric(plot_df[[color_cells_by]])) {
    # num_colors = unique(plot_df[[color_cells_by]]) %>% sort() %>% length()
    # full_spectrum_timepoint = get_colors(num_colors, type="rainbow")
    # names(full_spectrum_timepoint) = unique(plot_df$timepoint) %>% sort()

    gp <- gp +
      geom_point(
        data = plot_df,
        aes(umap2D_1, umap2D_2,
          color = color_cells_by
        ),
        size = cell_size,
        alpha = alpha,
        stroke = 0
      ) +
      # scale_color_gradientn(colors = full_spectrum_timepoint)+
      viridis::scale_color_viridis(option = "C")
  } else {
    if (is.null(color_values)) {
      num_colors <- unique(plot_df[[color_cells_by]]) %>%
        sort() %>%
        length()
      color_values <- get_colors(num_colors, "vibrant")
      names(color_values) <- unique(plot_df[[color_cells_by]]) %>% sort()
    }

    gp <- gp + geom_point(
      data = plot_df,
      aes(umap2D_1, umap2D_2,
        color = color_cells_by
      ),
      size = cell_size,
      alpha = alpha,
      stroke = 0
    ) +
      scale_color_manual(values = color_values)
  }

  # add labels

  if (plot_labels) {
    coord_matrix <- as.data.frame(reducedDim(cds, "UMAP"))
    grp_assign <- cds@colData %>%
      as.data.frame() %>%
      dplyr::select(all_of(color_cells_by))
    coord_matrix <- cbind(grp_assign, coord_matrix[row.names(grp_assign), ])
    colnames(coord_matrix)[1] <- "cell_group"
    centroid_coords <- aggregate(. ~ cell_group, data = coord_matrix, FUN = mean)

    if (dim(coord_matrix)[2] == 3) {
      colnames(centroid_coords) <- c(color_cells_by, "umap2D_1", "umap2D_2")
    } else {
      colnames(centroid_coords) <- c(color_cells_by, "umap2D_1", "umap2D_2", "umap2D_3")
    }

    if (repel_labels) {
      gp <- gp + ggrepel::geom_label_repel(
        data = centroid_coords,
        mapping = aes(get(paste0("umap2D_", x)),
          get(paste0("umap2D_", y)),
          label = get(color_cells_by)
        ),
        size = I(group_label_size),
        fill = "white"
      )
    } else {
      gp <- gp + geom_text(
        data = centroid_coords,
        mapping = aes(get(paste0("umap2D_", x)),
          get(paste0("umap2D_", y)),
          label = get(color_cells_by)
        ),
        size = I(group_label_size)
      )
    }
  }

  if (is.null(lab_title)) {
    lab_title <- color_cells_by
  }

  gp <- gp + labs(color = lab_title) + xlab(paste0("UMAP", x)) + ylab(paste0("UMAP", y))

  return(gp)
}

# plot a subset of labels on the plot
# returns a label df, use as follows:
# ggrepel::geom_text_repel(data = label_df, aes(label=cell_type, x=x, y=y), size=3)
my_plot_labels <- function(p, cds, x = 1, y = 2, relevant_cell_types = NULL) {
  colData(cds)$umap3d_1 <- reducedDims(cds)[["UMAP"]][, x]
  colData(cds)$umap3d_2 <- reducedDims(cds)[["UMAP"]][, y]
  label_df <- colData(cds) %>%
    as.data.frame() %>%
    group_by(projection_group, cell_type) %>%
    summarise(x = mean(umap3d_1), y = mean(umap3d_2))

  if (is.null(relevant_cell_types) == FALSE) {
    label_df <- label_df %>%
      filter(cell_type %in% relevant_cell_types)
  }

  return(label_df)
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
    ccs <- data
  } else if (class(data) == "cell_count_model") {
    ccs <- data@ccs
  }

  # if (is.null(switch_label) & !is.null(unique(igraph::V(state_graph)$cell_group)) & is(state_graph, "igraph")) {
  #   switch_label = unique(igraph::V(state_graph)$cell_group)
  # } else {
  #   print("define a switch label grouping that matches the ccm")
  # }

  if (is(state_graph, "igraph")) {
    path_df <- state_graph %>% igraph::as_data_frame()
  } else {
    path_df <- state_graph
  }

  if (is.null(switch_label)) {
    node_intersection <- intersect(rownames(ccs), unique(path_df$from, path_df$to))
  } else {
    node_intersection <- intersect(ccs@colData[[switch_label]], unique(path_df$from, path_df$to))
  }

  assertthat::assert_that(
    tryCatch(length(node_intersection) > 0,
      error = function(e) FALSE
    ),
    msg = "specify a label switch, path nodes and ccm cell group do not match"
  )


  gp <- my_plot_cells(data, x = x, y = y, color_cells_by = switch_label, ...)

  if (class(data) == "cell_count_set") {
    umap_centers <- centroids(data, switch_group = switch_label)
  } else if (class(data) == "cell_count_model") {
    umap_centers <- centroids(data@ccs, switch_group = switch_label)
  }



  path_df <- add_umap_coords(path_df, umap_centers)

  if (directed) {
    if (is.null(color_path_by) == FALSE) {
      gp <- gp +
        ggnewscale::new_scale_color() +
        geom_segment(
          data = path_df,
          aes(
            x = get(paste0("umap_to_", x)),
            y = get(paste0("umap_to_", y)),
            xend = get(paste0("umap_from_", x)),
            yend = get(paste0("umap_from_", y)),
            color = get(color_path_by)
          ),
          size = edge_size
        ) +
        geom_segment(
          data = path_df,
          aes(
            x = umap_from_1,
            y = umap_from_2,
            xend = (get(paste0("umap_to_", x)) + get(paste0("umap_from_", x))) / 2,
            yend = (get(paste0("umap_to_", y)) + get(paste0("umap_from_", y))) / 2,
            color = get(color_path_by)
          ),
          size = edge_size,
          linejoin = "mitre",
          arrow = arrow(type = "closed", angle = 30, length = unit(1, "mm"))
        ) +
        scale_color_gradient2(
          low = "#122985",
          mid = "white",
          high = "red4",
          na.value = "white"
        )
    } else {
      gp <- gp + geom_segment(
        data = path_df,
        aes(
          x = get(paste0("umap_to_", x)),
          y = get(paste0("umap_to_", y)),
          xend = get(paste0("umap_from_", x)),
          yend = get(paste0("umap_from_", y))
        ),
        color = path_color,
        size = edge_size
      ) +
        geom_segment(
          data = path_df,
          aes(
            x = get(paste0("umap_from_", x)),
            y = get(paste0("umap_from_", y)),
            xend = (get(paste0("umap_to_", x)) + get(paste0("umap_from_", x))) / 2,
            yend = (get(paste0("umap_to_", y)) + get(paste0("umap_from_", y))) / 2
          ),
          size = edge_size,
          color = path_color,
          linejoin = "mitre",
          arrow = arrow(type = "closed", angle = 30, length = unit(1, "mm"))
        )
    }
  } else {
    if (is.null(color_path_by) == FALSE) {
      gp <- gp +
        ggnewscale::new_scale_color() +
        geom_segment(
          data = path_df,
          aes(
            x = get(paste0("umap_to_", x)),
            y = get(paste0("umap_to_", y)),
            xend = get(paste0("umap_from_", x)),
            yend = get(paste0("umap_from_", y)),
            color = get(color_path_by)
          ),
          size = edge_size
        ) +
        scale_color_gradient2(
          low = "#122985",
          mid = "white",
          high = "red4",
          na.value = "white"
        )
    } else {
      gp <- gp + geom_segment(
        data = path_df,
        aes(
          x = get(paste0("umap_to_", x)),
          y = get(paste0("umap_to_", y)),
          xend = get(paste0("umap_from_", x)),
          yend = get(paste0("umap_from_", y))
        ),
        color = path_color,
        size = edge_size
      )
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
  path_df <- add_umap_coords(path_df, umap_centers)

  gp <- gp +
    geom_segment(
      data = path_df,
      aes(
        x = umap_to_1,
        y = umap_to_2,
        xend = umap_from_1,
        yend = umap_from_2
      ),
      size = size,
      color = color
    ) +
    geom_segment(
      data = path_df,
      aes(
        x = umap_from_1,
        y = umap_from_2,
        xend = (umap_to_1 + umap_from_1) / 2,
        yend = (umap_to_2 + umap_from_2) / 2
      ),
      size = size,
      color = color,
      linejoin = "mitre",
      arrow = arrow(type = "closed", angle = 30, length = unit(1, "mm"))
    )

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
    ccs <- data
    plot_df <- as.data.frame(ccs@cds_coldata)
  } else if (class(data) == "cell_count_model") {
    ccs <- data@ccs
  } else {
    print("some error message")
  }

  umap_centers <- centroids(ccs)
  cell_groups <- ccs@metadata[["cell_group_assignments"]] %>%
    pull(cell_group) %>%
    unique()
  nodes <- data.frame(id = cell_groups)
  if ("source" %in% colnames(edges)) {
    edges <- edges %>% rename("from" = source, "to" = target)
  }
  edges <- edges %>% select(from, to)

  no_edge <- setdiff(nodes$id, union(edges$from, edges$to))
  edges <- rbind(edges, data.frame("from" = no_edge, "to" = no_edge))
  n <- network::network(edges %>% dplyr::select(from, to), directed = T, loops = T)
  nodes <- nodes[match(network::network.vertex.names(n), nodes$id), ]
  n %v% "id" <- network::network.vertex.names(n)

  merged_coords <- ggnetwork(n) %>%
    select(id, vertex.names) %>%
    unique() %>%
    left_join(umap_centers, by = c("id" = "cell_group"))
  rownames(merged_coords) <- merged_coords$id
  coords <- merged_coords[network::network.vertex.names(n), ] %>% select(umap_1, umap_2)
  geo <- as.matrix(sapply(coords, as.numeric))

  g <- ggnetwork(x = n, layout = geo, arrow.gap = arrow.gap, scale = scale)

  show(ggplot(g, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    geom_nodes(size = 7, colour = "black", shape = 21) +
    geom_nodetext_repel(aes(label = id), size = 3) +
    theme_blank())
}

#' Plot Cells Per Sample
#'
#' This function generates a boxplot of cell counts per sample, with options to customize the plot.
#'
#' @param ccs A data frame or tibble containing cell count data.
#' @param x_col A string specifying the column name to be used for the x-axis.
#' @param y_col A string specifying the column name to be used for the y-axis. Default is "count".
#' @param cell_groups A vector of cell group names to filter the data. Default is an empty vector.
#' @param color_by A string specifying the column name to be used for coloring the plot. Default is "cell_group".
#' @param plot_zeroes A logical value indicating whether to plot zero counts. Default is FALSE.
#' @param plot_points A logical value indicating whether to plot individual points. Default is FALSE.
#' @param log_scale A logical value indicating whether to use a logarithmic scale for the y-axis. Default is FALSE.
#' @param nrow An integer specifying the number of rows for faceting. Default is 1.
#' @param legend_position A string specifying the position of the legend. Default is "none".
#'
#' @return A ggplot object representing the boxplot of cell counts per sample.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang sym
#' @importFrom monocle3 monocle_theme_opts
#' @export
plot_cells_per_sample <- function(ccs,
                                  x_col,
                                  y_col = c("count", "count_per_1000"),
                                  cell_groups = c(),
                                  batch_col = NULL, 
                                  color_by = "cell_group",
                                  plot_zeroes = T,
                                  plot_points = F,
                                  facet = T, 
                                  log_scale = F,
                                  nrow = 1,
                                  legend_position = "none") {
  
  y_col <- match.arg(y_col)
  
  
  if (is.null(batch_col) == FALSE) {
    batches_to_keep = colData(ccs) %>% as.data.frame %>% 
      group_by(!!sym(batch_col), !!sym(x_col)) %>% 
      tally() %>% 
      group_by(!!sym(batch_col)) %>% 
      filter(n()>1) %>% 
      pull(!!sym(batch_col))
    
    ccs = subset_ccs(ccs, !!sym(batch_col) %in% batches_to_keep)
    
  }
  
  count_df <- get_norm_df(ccs)

  if (length(cell_groups) != 0) {
    count_df <- count_df %>% filter(cell_group %in% cell_groups)
  }

  if (plot_zeroes) {
    count_df <- count_df %>% mutate(count = count + 0.001)
  }

  count_df[[x_col]] <- as.factor(count_df[[x_col]])

  p <- count_df %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col), fill = !!sym(color_by))) +
    geom_boxplot() +
    facet_wrap(~cell_group, nrow = nrow) +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()

  if (plot_points) {
    p <- p + geom_jitter(aes(x = !!sym(x_col), y = !!sym(y_col)), size = 1)
  }

  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  if (facet) {
    p <- p + facet_wrap(~cell_group, scale = "free")
  }

  return(p)
}

plot_conditional_cells_per_sample <- function(ccm,
                                              x_col,
                                              y_col = c("count", "count_per_1000"),
                                              cell_groups = c(),
                                              batch_col = NULL, 
                                              color_by = "cell_group",
                                              plot_zeroes = T,
                                              plot_points = F,
                                              facet = T, 
                                              log_scale = F,
                                              nrow = 1,
                                              legend_position = "none", 
                                              newdata = tibble()) {
  
  y_col <- match.arg(y_col)
  
  sample_metadata <- colData(ccm@ccs) %>% tidyr::as_tibble()
  sel_ccs_counts <- monocle3::normalized_counts(ccm@ccs, 
                                                norm_method = "size_only", 
                                                pseudocount = 0)
  
  # override the columns with the newdata columns
  for (c in colnames(newdata)) {
    sample_metadata[[c]] <- newdata[[c]]
  }
  
  conditional_counts <- hooke::estimate_abundances_cond(ccm,
                                                        newdata = sample_metadata,
                                                        cond_responses = sel_ccs_counts,
                                                        pln_model = "reduced"
  )
  
  count_df <- conditional_counts %>%
    dplyr::select(sample, cell_group, log_abund) %>%
    mutate(num_cells = exp(log_abund)) %>%
    select(-log_abund) %>% 
    mutate(count = as.integer(num_cells))
  
  count_df <- left_join(count_df, 
                        colData(ccm@ccs) %>% as.data.frame, by = "sample")
  
  cell_type_total <- Matrix::colSums(counts(ccm@ccs))
  geometric_mean <- exp(mean(log(cell_type_total)))
  
  count_df <- count_df %>% mutate(count_per_1000 = count * 1000 / geometric_mean)
  
  if (length(cell_groups) != 0) {
    count_df <- count_df %>% filter(cell_group %in% cell_groups)
  }
  
  if (plot_zeroes) {
    count_df <- count_df %>% mutate(count = count + 0.001)
  }
  
  count_df[[x_col]] <- as.factor(count_df[[x_col]])
  
  p <- count_df %>%
    ggplot(aes(x = !!sym(x_col), y = !!sym(y_col), fill = !!sym(color_by))) +
    geom_boxplot() +
    facet_wrap(~cell_group, nrow = nrow) +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts()
  
  if (plot_points) {
    p <- p + geom_jitter(aes(x = !!sym(x_col), y = !!sym(y_col)), size = 1)
  }
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  if (facet) {
    p <- p + facet_wrap(~cell_group, scale = "free")
  }
  
  return(p)
}

#' @noRd
plot_cells_highlight <- function(ccs, group_to_highlight, colname) {
  plot_df <- as.data.frame(ccs@cds_coldata)
  plot_df$cell_group <- ccs@cds_coldata[[colname]]

  plot_df$cell <- row.names(plot_df)
  plot_df$umap2D_1 <- ccs@cds_reduced_dims[["UMAP"]][plot_df$cell, x]
  plot_df$umap2D_2 <- ccs@cds_reduced_dims[["UMAP"]][plot_df$cell, y]

  gp <- ggplot() +
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
                             scale_shifts_by = c("receiver", "sender", "none"),
                             edge_size = 2,
                             cell_size = 25,
                             q_value_thresh = 1.0,
                             group_label_size = 2,
                             plot_labels = c("significant", "all", "none"),
                             fc_limits = c(-3, 3),
                             sender_cell_groups = NULL,
                             receiver_cell_groups = NULL,
                             plot_edges = c("all", "directed", "undirected", "none"),
                             label_cell_groups = list(),
                             repel_labels = TRUE,
                             model_for_pcors = "reduced",
                             switch_label = NULL,
                             sub_cds = NULL,
                             alpha = 1.0,
                             x = 1,
                             y = 2) {
  plot_df <- ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell <- row.names(plot_df)

  plot_df$umap3D_1 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell, 1]
  plot_df$umap3D_2 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell, 2]
  plot_df$umap3D_3 <- ccm@ccs@cds_reduced_dims[["UMAP"]][plot_df$cell, 3]

  cond_b_vs_a_tbl <- cond_b_vs_a_tbl %>% dplyr::mutate(delta_log_abund = ifelse(delta_q_value <= q_value_thresh, delta_log_abund, 0))


  plot_df <- dplyr::left_join(plot_df,
    cond_b_vs_a_tbl,
    by = c("cell_group" = "cell_group")
  )

  fc_limits <- range(plot_df$delta_log_abund)
  range_min <- min(abs(fc_limits))

  plot_df <- plot_df %>%
    mutate(delta_log_abund = ifelse(delta_log_abund > range_min, max, delta_log_abund)) %>%
    mutate(delta_log_abund = ifelse(delta_log_abund < -range_min, -range_min, delta_log_abund))

  # if (is.null(fc_limits)) {
  #   fc_limits = range(plot_df$delta_log_abund)
  # } else {
  #   min = fc_limits[1]
  #   max = fc_limits[2]
  #   plot_df = plot_df %>%
  #     mutate(delta_log_abund = ifelse(delta_log_abund > max, max, delta_log_abund)) %>%
  #     mutate(delta_log_abund = ifelse(delta_log_abund < min, min, delta_log_abund))
  # }

  color_palette <- c("#3A5FCD", "#FAFAFA", "#CD3700")
  # color_palette = c('#3A5FCD', '#FFFFFF','#CD3700')
  # color_palette = c('#FFFFFF','#CD3700')

  p <- plotly::plot_ly(plot_df,
    x = ~umap3D_1,
    y = ~umap3D_2,
    z = ~umap3D_3,
    type = "scatter3d",
    mode = "markers",
    alpha = I(alpha),
    size = I(cell_size),
    color = ~delta_log_abund,
    colors = color_palette
    # marker = list(size = I(cell_size),
    #               color = ~delta_log_abund,
    #               colorscale = colorscale)
  )

  return(p)
}

#' Default plotting options for ggplot2
#'
#' return A ggplot2 theme object.
#' @export
hooke_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(panel.border = element_blank()) +
    # theme(axis.line.x = element_line(size=0.25, color="black")) +
    # theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank())
}
