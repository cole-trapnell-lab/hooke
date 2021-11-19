#' Plot a UMAP colored by how cells shift in a given contrast
#'
#' @param ccm A cell_count_model object.
#' @param criterion a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".
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
                          cell_size=1){

  umap_centers = centroids(ccm@ccs)

  corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                              umap_centers,
                                                              cond_b_vs_a_tbl,
                                                              log_abundance_thresh)
  directed_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("directed_to_from", "directed_from_to"))
  undirected_edge_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(edge_type %in% c("undirected"))
  
  umap_centers_delta_abund = umap_centers
  umap_centers_delta_abund = dplyr::left_join(umap_centers_delta_abund, cond_b_vs_a_tbl, by=c("cell_group"="cell_group"))
  umap_centers_delta_abund = umap_centers_delta_abund %>% dplyr::mutate(max_log_abund = pmax(log_abund_x, log_abund_y))
  umap_centers_delta_abund = umap_centers_delta_abund %>%
    dplyr::mutate(max_log_abund = ifelse(max_log_abund < log_abundance_thresh, log_abundance_thresh, max_log_abund))


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

  #cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% mutate(cluster = stringr::str_split_fixed(cell_group, "\\.", 3)[,3])
  plot_df = dplyr::left_join(plot_df,
                      cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                      by=c("cell_group"="cell_group"))

  #directed_edge_df = directed_edge_df %>% filter(edge_type %in% c("directed_to_from", "directed_from_to"))

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
      na.value = "white"
    )  + 
    theme_void() +
    theme(legend.position = "none") +

    #plot_oriented_pln_network(directed_edge_df) +

    #geom_point(data = umap_centers_delta_abund, aes(umap2D_1, umap2D_2, color=delta_log_abund, size=max_log_abund)) +
    # #geom_text(data = umap_centers_delta_abund, aes(umap2D_1, umap2D_2, color=delta_log_abund, label=label)) +
    #scale_color_gradient2(low = 'steelblue', mid = 'grey', high = 'darkred') +
    # #theme(legend.position = "none") +
    monocle3:::monocle_theme_opts() #+ ggtitle(paste(cond_b, "vs",cond_a, "pos partial corr w/ umap penalty"))

  gp = gp  +
    geom_segment(data = undirected_edge_df,
                 aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=umap_from_1,
                     yend = umap_from_2,
                     size=edge_size * scaled_weight),
                 #size=edge_size / 4,
                 color="lightgray") +
    geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                 aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=umap_from_1,
                     yend = umap_from_2,
                     size=edge_size * scaled_weight),
                 color="black") +
    geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                 aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=(umap_to_1+umap_from_1)/2,
                     yend = (umap_to_2+umap_from_2)/2,
                     size=edge_size * scaled_weight),
                 color="black",
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
    geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                 aes(x = umap_from_1,
                     y = umap_from_2,
                     xend=umap_to_1,
                     yend = umap_to_2,
                     size=edge_size * scaled_weight),
                 color="black") +
    geom_segment(data = directed_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                 aes(x = umap_from_1,
                     y = umap_from_2,
                     xend=(umap_from_1+umap_to_1)/2,
                     yend = (umap_from_2+umap_to_2)/2,
                     size=edge_size * scaled_weight),
                 color="black",
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
    scale_size_identity()
  
  return(gp)
}

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
                      cond_b_vs_a_tbl = NULL) {
  
  if (class(data) == "cell_count_set") {
    cds = data@cds
    plot_df = as.data.frame(colData(cds))
  } else if (class(data) == "cell_data_set") {
    cds = data
    plot_df = as.data.frame(colData(cds))
    plot_df = data@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else if (class(data) == "cell_count_model") {
    cds = data@ccs@cds
    plot_df = as.data.frame(colData(cds))
    plot_df$cell_group = data@ccs@metadata[["cell_group_assignments"]] %>% pull(cell_group)
  } else {
    print("some error message")
  }

  
  plot_df$cell = row.names(plot_df)
  plot_df$color_cells_by = plot_df[[color_cells_by]]
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[plot_df$cell,2]
  
  
  
  
  # could be an vector that corresponds to a label
  if (!is.null(residuals)) {
    res_df = as.data.frame(residuals) %>% rownames_to_column("cell_group")
    colnames(res_df)[2] = value
    plot_df = plot_df %>% left_join(res_df, by="cell_group")
    color_cells_by = "residuals"

  }

  if (!is.null(cond_b_vs_a_tbl)) {

    plot_df = dplyr::left_join(plot_df,
                               cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                               by=c("cell_group"="cell_group"))
    
    color_cells_by = "delta_log_abund"

  }

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

  } else if (color_cells_by == "residuals") {
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
        # limits = c(-3,3)
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
        # limits=c(-3,3)
      )
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


plot_path <- function(ccm,
                      path_df = path_df,
                      edge_size=2, 
                      path_color = "black", 
                      color_cells_by = "cluster", 
                      residuals = NULL, 
                      cond_b_vs_a_tbl = NULL) {


  gp = my_plot_cells(ccm, color_cells_by = color_cells_by, residuals = residuals, cond_b_vs_a_tbl = cond_b_vs_a_tbl)
  
  umap_centers = centroids(ccm@ccs)
  
  path_df = add_umap_coords(path_df, umap_centers)
  
  gp = gp + geom_segment(data = path_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=umap_from_1,
                   yend = umap_from_2),
               color=path_color,
               size=2 ) + 
    geom_segment(data = path_df,
                 aes(x = umap_from_1,
                     y = umap_from_2,
                     xend=(umap_to_1+umap_from_1)/2,
                     yend = (umap_to_2+umap_from_2)/2),
                 size=2,
                 color=path_color,
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))
  
  return(gp)

}


