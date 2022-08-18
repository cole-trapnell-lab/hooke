library(monocle3)
library(hooke)


# get WT time series data -----------------------------------------------------

# get a WT time series cds
# ref_cds <- readRDS(file = "R_objects/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")
# meso_cell_types = c("mature fast muscle",
#                     "mature slow muscle",
#                     "slow-committed myocyte",
#                     "fast-committed myocyte",
#                     "myoblast",
#                     "mesodermal progenitor cells (contains PSM)",
#                     "satellite cells")
#
# meso_time_cds = ref_cds[,colData(ref_cds)$cell_type_broad %in% meso_cell_types]
# meso_time_cds = meso_time_cds[, !is.na(meso_time_cds$timepoint)]
#
# meso_time_cds = preprocess_cds(meso_time_cds, num_dim = 50) %>%
#             align_cds(residual_model_formula_str = "~log.n.umi") %>%
#             reduce_dimension(max_components = 3, preprocess_method = "Aligned")
#
# meso_time_cds = cluster_cells(meso_time_cds, random_seed = 42, resolution = 1e-4)
#
# colData(meso_time_cds)$cluster = monocle3::clusters(meso_time_cds)
# colData(meso_time_cds)$partition = partitions(meso_time_cds)
# colData(meso_time_cds)$Size_Factor = size_factors(meso_time_cds)
#
# saveRDS(meso_time_cds, "my_R_objects/ref-meso-all-timepoint.rds")

#meso_time_cds <- readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/my_R_objects/ref-meso-all-timepoint.rds")


# I'm just going to use timepoints between 18 and 24 hpf here for now ---------

#meso_18_24_cds <- meso_time_cds[,colData(meso_time_cds)$timepoint <= 24]

meso_18_24_cds = meso_cds
#
# meso_18_24_cds = preprocess_cds(meso_18_24_cds, num_dim = 50) %>%
#   align_cds(residual_model_formula_str = "~log.n.umi") %>%
#   reduce_dimension(max_components = 3, preprocess_method = "Aligned")
#
# meso_18_24_cds = cluster_cells(meso_18_24_cds, random_seed = 42, resolution = 1e-4)
#
# colData(meso_18_24_cds)$cluster = monocle3::clusters(meso_18_24_cds)

plot_cells(meso_18_24_cds)


# run hooke -------------------------------------------------------------------

meso_ccs = new_cell_count_set(meso_18_24_cds,
                              sample_group = "Oligo",
                              cell_group = "cluster")

#meso_ccm = new_cell_count_model(meso_ccs,
#                                model_formula_str = "~ timepoint + gene_target")


# Now try it with a white and black list:
paga_graph = hooke:::get_paga_graph(meso_ccs@cds)
edge_whitelist = igraph::as_data_frame(paga_graph)
edge_blacklist = igraph::as_data_frame(igraph::complementer(paga_graph))

meso_ccm  = new_cell_count_model(meso_ccs,
                            model_formula_str = "~timepoint + gene_target",
                            whitelist=edge_whitelist,
                            base_penalty=1)

meso_ccm = select_model(meso_ccm, sparsity_factor = 1)

#plot(meso_ccm@best_model, output = "corrplot")

#coef(meso_ccm@best_model)

# do all combinations of contrasts --------------------------------------------

timepoints = unique(colData(meso_18_24_cds)$timepoint) %>% sort()
timepoint_combos = combn(timepoints, 2) %>%
  t() %>%
  as.data.frame() %>%
  filter(abs(V1-V2) > 2) # get rid of really close combos
timepoint_combos

# first the farthest timepoints ------------------------------------------------

time_18 = estimate_abundances(meso_ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
#time_20 = estimate_abundances(meso_ccm, data.frame(timepoint="20"))
#time_22 = estimate_abundances(meso_ccm, data.frame(timepoint="22"))
time_24 = estimate_abundances(meso_ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))


# large neg combo -------------------------------------------------------------

cond_18_vs_24_tbl = compare_abundances(meso_ccm, time_18, time_24)
plot_contrast(meso_ccm, cond_18_vs_24_tbl)

# all combos in between ---------------------------------------------------------
#
# cond_18_vs_20_tbl = compare_abundances(meso_ccm, time_18, time_20)
# plot_contrast(meso_ccm, cond_18_vs_20_tbl)
#
# cond_18_vs_22_tbl = compare_abundances(meso_ccm, time_18, time_22)
# plot_contrast(meso_ccm, cond_18_vs_22_tbl)
#
# cond_20_vs_22_tbl = compare_abundances(meso_ccm, time_20, time_22)
# plot_contrast(meso_ccm, cond_20_vs_22_tbl)
#
# cond_20_vs_24_tbl = compare_abundances(meso_ccm, time_20, time_24)
# plot_contrast(meso_ccm, cond_20_vs_24_tbl)
#
# cond_22_vs_24_tbl = compare_abundances(meso_ccm, time_22, time_24)
# plot_contrast(meso_ccm, cond_22_vs_24_tbl)


# -----------------------------------------------------------------------------

umap_centers = centroids(meso_ccm@ccs)

# -----------------------------------------------------------------------------

corr_edge_coords_umap_delta_abund_18_24 = hooke:::collect_pln_graph_edges(meso_ccm,
                                                                  cond_18_vs_24_tbl)
corr_edge_coords_umap_delta_abund = corr_edge_coords_umap_delta_abund_18_24

# corr_edge_coords_umap_delta_abund_18_20 = hooke::collect_pln_graph_edges(meso_ccm,
#                                                                   cond_18_vs_20_tbl,
#                                                                   log_abundance_thresh)



# -----------------------------------------------------------------------------

# max flow

calc_max_flow <- function(edges, source, target) {

  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  mf = igraph::max_flow(G, source = source, target = target, capacity = igraph::E(G)$pcor)
  igraph::E(G)$flow <- mf$flow
  igraph::V(G)$pass <- igraph::strength(G,mode="in",weights=mf$flow)

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(G)

  # if flow values are neg, reverse them
  switch_dir_df = new_pos_edge_coords_df %>% dplyr::filter(flow < 0) %>%
    dplyr::rename("from"="to",
                  "umap_from_1"="umap_to_1",
                  "umap_from_2"="umap_to_2",
                  "to"="from",
                  "umap_to_1"="umap_from_1",
                  "umap_to_2"="umap_from_2") %>%
    dplyr::mutate(flow = -flow)
  max_flow_edge_df = rbind(dplyr::filter(new_pos_edge_coords_df, flow >=0),
                           switch_dir_df)
  return(max_flow_edge_df)

}


# find positive correlations
pos_edge_coords_df = corr_edge_coords_umap_delta_abund %>% dplyr::filter(pcor > 0 )

max_flow_edge_df = calc_max_flow(pos_edge_coords_df, source = "1", target = "4")



calc_shortest_path <- function(edges, source, target) {
  edges = edges %>% dplyr::select(from, to, from_delta_log_abund, to_delta_log_abund, pcor)
  edges = edges %>% dplyr::mutate(weight = abs(from_delta_log_abund - to_delta_log_abund) / pcor)
  G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))

  mf = igraph::shortest_paths(G, from = source, to = target, weights = igraph::E(G)$weight, output="epath")
  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])

  #igraph::E(G)$flow <- mf$flow
  #igraph::V(G)$pass <- igraph::strength(G,mode="in",weights=mf$flow)

  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph) %>% dplyr::select(from, to, weight)

  # if flow values are neg, reverse them
  # switch_dir_df = new_pos_edge_coords_df %>% dplyr::filter(flow < 0) %>%
  #   dplyr::rename("from"="to",
  #                 "umap_from_1"="umap_to_1",
  #                 "umap_from_2"="umap_to_2",
  #                 "to"="from",
  #                 "umap_to_1"="umap_from_1",
  #                 "umap_to_2"="umap_from_2")
  #max_flow_edge_df = rbind(dplyr::filter(new_pos_edge_coords_df, flow >=0),
  #                         switch_dir_df)
  return(new_pos_edge_coords_df)
}
debug(calc_shortest_path)
sp_edge_df = calc_shortest_path(pos_edge_coords_df, source = "1", target = "4")

#corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
sp_edge_df = dplyr::left_join(sp_edge_df, cond_b_vs_a_tbl %>% setNames(paste0('to_', names(.))), by=c("to"="to_cell_group")) #%>%
#dplyr::rename(log_abund_x,
#              to_delta_log_abund = delta_log_abund)
sp_edge_df = dplyr::left_join(sp_edge_df, cond_b_vs_a_tbl %>% setNames(paste0('from_', names(.))), by=c("from"="from_cell_group"))# %>%
#  dplyr::rename(from_delta_log_abund = delta_log_abund)

sp_edge_df_umap_delta_abund = dplyr::left_join(sp_edge_df,
                                                     umap_centers %>% dplyr::select("from" = cell_group,
                                                                                    "umap_from_1" = V1,
                                                                                    "umap_from_2" = V2))
sp_edge_df_umap_delta_abund = dplyr::left_join(sp_edge_df_umap_delta_abund ,
                                                     umap_centers %>% dplyr::select("to" = cell_group,
                                                                                    "umap_to_1" = V1,
                                                                                    "umap_to_2" = V2))
max_flow_edge_df = sp_edge_df_umap_delta_abund
max_flow_edge_df$flow = 1/max_flow_edge_df$weight

plot_cells(meso_18_24_cds, color_cells_by = "cluster")

return_baseplot <- function(ccm,
                            cond_b_vs_a_tbl,
                            cell_size=1,
                            edge_size=2) {


  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)

  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]

  #cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% mutate(cluster = stringr::str_split_fixed(cell_group, "\\.", 3)[,3])
  plot_df = dplyr::left_join(plot_df,
                             cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                             by=c("cell_group"="cell_group"))

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
        dplyr::arrange(!is.na(abs(delta_log_abund)),
                abs(delta_log_abund)),
      aes(umap2D_1, umap2D_2),
      size = cell_size,
      stroke = 0
    ) +
    scale_color_gradient2(
      low = "#122985",
      mid = "white",
      high = "red4",
      na.value = "white"
    ) +
    theme_void() +
    theme(legend.position = "none") +
    monocle3:::monocle_theme_opts()
  return(gp)
}

bp = return_baseplot(meso_ccm, cond_18_vs_24_tbl)
bp + geom_segment(data = pos_edge_coords_df,
                  aes(x = umap_to_1,
                      y = umap_to_2,
                      xend=umap_from_1,
                      yend = umap_from_2),
                  color="black")


ggplot(data = max_flow_edge_df %>% dplyr::arrange(-flow) %>% head(10)) +
  geom_segment(
    aes(x = umap_from_1,
        y = umap_from_2,
        xend=(umap_to_1+umap_from_1)/2,
        yend = (umap_to_2+umap_from_2)/2,
        color=flow),
    # color="black",
    linejoin='mitre',
    arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
  geom_segment(
    aes(x = umap_from_1,
        y = umap_from_2,
        xend=umap_to_1,
        yend = umap_to_2,
        color=flow),
    # color="black"
  ) +
  monocle3:::monocle_theme_opts()



bp +
  geom_segment(
    data = max_flow_edge_df %>% dplyr::arrange(-flow) %>% dplyr::filter(flow > 0),
    aes(x = umap_from_1,
        y = umap_from_2,
        xend=(umap_to_1+umap_from_1)/2,
        yend = (umap_to_2+umap_from_2)/2,
        color=flow),
    # color="black",
    linejoin='mitre',
    arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
  geom_segment(
    data = max_flow_edge_df %>% dplyr::arrange(-flow) %>% dplyr::filter(flow > 0),
    aes(x = umap_from_1,
        y = umap_from_2,
        xend=umap_to_1,
        yend = umap_to_2,
        color=flow),
    # color="black"
  ) +
  monocle3:::monocle_theme_opts()


# cMST ------------------------------------------------------------------------


calc_mst <- function(edges) {
  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  G_mst <- mst(G, weights=igraph::E(G)$pcor)
  mst_df <- igraph::as_data_frame(G_mst)
  return(mst_df)
}


mst_df <- calc_mst(pos_edge_coords_df)

bp +
  geom_segment(data = mst_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=(umap_to_1+umap_from_1)/2,
                   yend = (umap_to_2+umap_from_2)/2),
               color="black",
               linejoin='mitre',
               arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
  geom_segment(data = mst_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=umap_from_1,
                   yend = umap_from_2,
                   color=flow),
               color="black") +
  monocle3:::monocle_theme_opts()


