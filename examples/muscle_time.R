library(monocle3)
library(hooke)
setwd("~/OneDrive/UW/Trapnell/hooke/")
# devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")
# actually i should switch the branch of this
# devtools::load_all("~/OneDrive/UW/Trapnell/bin/monocle3-dev/")
# devtools::load_all("~/OneDrive/UW/Trapnell/bin/monocle3-dev-test/")

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
meso_time_cds = preprocess_cds(meso_time_cds, num_dim = 50) %>%
            align_cds(residual_model_formula_str = "~log.n.umi") %>%
            reduce_dimension(max_components = 3, 
                             preprocess_method = "Aligned", 
                             build_nn_index = T)
 
meso_time_cds = cluster_cells(meso_time_cds, random_seed = 42, resolution = 1e-4)
# meso_time_cds@preprocess_aux$gene_loadings = NULL
# meso_time_cds@preprocess_aux$prop_var_expl = NULL
# meso_time_cds@preprocess_aux$beta = NULL

# save_transform_models(meso_time_cds, "../my_R_objects/ref-meso-all-timepoint-models-new")
 
# colData(meso_time_cds)$cluster = monocle3::clusters(meso_time_cds)
# colData(meso_time_cds)$partition = partitions(meso_time_cds)
# colData(meso_time_cds)$Size_Factor = size_factors(meso_time_cds)
# 
# saveRDS(meso_time_cds, "my_R_objects/ref-meso-all-timepoint.rds")

meso_time_cds <- readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/my_R_objects/ref-meso-all-timepoint.rds")


# I'm just going to use timepoints between 18 and 24 hpf here for now ---------

# meso_18_24_cds <- meso_time_cds[,colData(meso_time_cds)$timepoint <= 24]
# 
# meso_18_24_cds = preprocess_cds(meso_18_24_cds, num_dim = 50) %>%
#               align_cds(residual_model_formula_str = "~log.n.umi") %>%
#               reduce_dimension(max_components = 3, 
#                                preprocess_method = "Aligned",
#                                build_nn_index = T)
# 
# save_transform_models(meso_18_24_cds, "../my_R_objects/ref-meso-18-24-models")
# 
# meso_18_24_cds = cluster_cells(meso_18_24_cds, random_seed = 42, resolution = 1e-4)
# 
# colData(meso_18_24_cds)$cluster = monocle3::clusters(meso_18_24_cds)
# 
# plot_cells(meso_18_24_cds) %>% ggsave(filename = "examples/meso_18_24_by_cluster.png")

# saveRDS(meso_18_24_cds, "../gap-notebook-ct-1/my_R_objects/meso_18_24_cds.rds")

meso_18_24_cds <- readRDS("../gap-notebook-ct-1/my_R_objects/meso_18_24_cds.rds")


# project the mutant data into wt space and transfer cluster labels -----------




# run hooke -------------------------------------------------------------------

meso_ccs = new_cell_count_set(meso_18_24_cds,
                              sample_group = "embryo",
                              cell_group = "cluster")

meso_ccm = new_cell_count_model(meso_ccs,
                                model_formula_str = "~ as.factor(timepoint)")


plot(meso_ccm@best_model, output = "corrplot")



# do all combinations of contrasts --------------------------------------------

timepoints = unique(colData(meso_18_24_cds)$timepoint) %>% sort()
timepoint_combos = combn(timepoints, 2) %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    filter(abs(V1-V2) > 2) # get rid of really close combos
timepoint_combos

# first the farthest timepoints ------------------------------------------------

time_18 = estimate_abundances(meso_ccm, data.frame(timepoint="18"))
time_20 = estimate_abundances(meso_ccm, data.frame(timepoint="20"))
time_22 = estimate_abundances(meso_ccm, data.frame(timepoint="22"))
time_24 = estimate_abundances(meso_ccm, data.frame(timepoint="24"))


# large neg combo -------------------------------------------------------------

cond_18_vs_24_tbl = compare_abundances(meso_ccm, time_18, time_24)
plot_contrast(meso_ccm, cond_18_vs_24_tbl)

# all combos in between ---------------------------------------------------------

cond_18_vs_20_tbl = compare_abundances(meso_ccm, time_18, time_20)
plot_contrast(meso_ccm, cond_18_vs_20_tbl)

cond_18_vs_22_tbl = compare_abundances(meso_ccm, time_18, time_22)
plot_contrast(meso_ccm, cond_18_vs_22_tbl)

cond_20_vs_22_tbl = compare_abundances(meso_ccm, time_20, time_22)
plot_contrast(meso_ccm, cond_20_vs_22_tbl)

cond_20_vs_24_tbl = compare_abundances(meso_ccm, time_20, time_24)
plot_contrast(meso_ccm, cond_20_vs_24_tbl)

cond_22_vs_24_tbl = compare_abundances(meso_ccm, time_22, time_24)
plot_contrast(meso_ccm, cond_22_vs_24_tbl)


# -----------------------------------------------------------------------------

umap_centers = centroids(meso_ccm@ccs)



# create a matrix Z
# zero out the values in the WT pcor matrix 
# that aren't included in the green edges

# zero out bottom half 
zero_out_pcor_matrix <- function(cmm, green_edges=NULL) {
  
  pcor_matrix = as.matrix(cmm@best_model$latent_network(type="partial_cor"))
  
  if (!is.null(green_edges)) {
    pcor_matrix = pcor_matrix %>%
      as.data.frame() %>%
      rownames_to_column("from") %>%
      pivot_longer(-"from", names_to = "to") %>%
      mutate(value = case_when(
        from %in% green_edges$from & to %in% green_edges$to ~ value,
        TRUE ~ 0
      )) %>%
      pivot_wider(names_from = "to", values_from = value) %>%
      column_to_rownames("from")
  } else {
    pcor_matrix[lower.tri(pcor_matrix)] = 0
  }
  
  
  return(pcor_matrix)
}

return_baseplot <- function(ccm, 
                            cond_b_vs_a_tbl, 
                            cell_size=1,
                            edge_size=2, 
                            legend_position="none") {
  
  
  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)
  
  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]
  
  #cond_b_vs_a_tbl = cond_b_vs_a_tbl %>% mutate(cluster = stringr::str_split_fixed(cell_group, "\\.", 3)[,3])
  plot_df = dplyr::left_join(plot_df,
                             cond_b_vs_a_tbl %>% dplyr::select(cell_group, delta_log_abund),
                             by=c("cell_group"="cell_group"))
  
  # geom_segment(data = pos_edge_coords_df,
  #              aes(x = umap_to_1,
  #                  y = umap_to_2,
  #                  xend=umap_from_1,
  #                  yend = umap_from_2),
  #              color="black")
  
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
      # limits=c(-3,3)
    ) +
    theme_void() +
    theme(legend.position = legend_position) +
    monocle3:::monocle_theme_opts() 
  return(gp)
}

# -----------------------------------------------------------------------------

corr_edge_coords_umap_delta_abund_18_24 = collect_pln_graph_edges(meso_ccm,
                                                            umap_centers,
                                                            cond_18_vs_24_tbl,
                                                            log_abundance_thresh)


corr_edge_coords_umap_delta_abund_18_20 = collect_pln_graph_edges(meso_ccm,
                                                                umap_centers,
                                                                cond_18_vs_20_tbl,
                                                                log_abundance_thresh)



# -----------------------------------------------------------------------------

# max flow 

calc_max_flow <- function(edges, source, target) {
  
  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  mf = max_flow(G, source = source, target = target, capacity = igraph::E(G)$pcor)
  igraph::E(G)$flow <- mf$flow
  igraph::V(G)$pass <- igraph::strength(G,mode="in",weights=mf$flow)
  
  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(G)
  
  # if flow values are neg, reverse them
  switch_dir_df = new_pos_edge_coords_df %>% filter(flow < 0) %>% 
    dplyr::rename("from"="to", 
                  "umap_from_1"="umap_to_1",
                  "umap_from_2"="umap_to_2",
                  "to"="from", 
                  "umap_to_1"="umap_from_1", 
                  "umap_to_2"="umap_from_2") %>% 
    mutate(flow = -flow)
  max_flow_edge_df = rbind(filter(new_pos_edge_coords_df, flow >=0), 
                           switch_dir_df) %>%
                     filter(flow > 0)
  return(max_flow_edge_df)
  
}

max_flow_edge_df = calc_max_flow(pos_edge_coords_df, source = "4", target = "1")


# find positive correlations
pos_edge_coords_df = corr_edge_coords_umap_delta_abund_18_24 %>% filter(pcor > 0 )

plot_cells(meso_18_24_cds, color_cells_by = "cluster")

bp = return_baseplot(meso_ccm, cond_18_vs_24_tbl)
bp + geom_segment(data = pos_edge_coords_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=umap_from_1,
                   yend = umap_from_2),
               color="black")  


ggplot(data = max_flow_edge_df %>% arrange(-flow) %>% head(10)) + 
  geom_segment(
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=(umap_to_1+umap_from_1)/2,
                   yend = (umap_to_2+umap_from_2)/2, 
                   color=flow),
               # color="black",
               linejoin='mitre',
               arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
  geom_segment(
                      aes(x = umap_to_1,
                          y = umap_to_2,
                          xend=umap_from_1,
                          yend = umap_from_2, 
                          color=flow), 
                      # color="black"
               ) + 
  monocle3:::monocle_theme_opts() 


# MST ------------------------------------------------------------------------


calc_mst <- function(edges) {
  G <- igraph::graph_from_data_frame(edges, directed=FALSE)
  G_mst <- mst(G, weights=igraph::E(G)$pcor)
  mst_df <- igraph::as_data_frame(G_mst)
  return(mst_df)
}


mst_df <- calc_mst(pos_edge_coords_df)

plot_edges <- function(df, bp) {
  
  show(
    # bp +
    ggplot(data = df)+
    geom_segment(aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=(umap_to_1+umap_from_1)/2,
                     yend = (umap_to_2+umap_from_2)/2),
                 color="black",
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
    geom_segment( aes(x = umap_to_1,
                     y = umap_to_2,
                     xend=umap_from_1,
                     yend = umap_from_2, 
                     color=flow), 
                 color="black") + 
    monocle3:::monocle_theme_opts() )
}



# try with paga also ----------------------------------------------------------

paga_graph = hooke:::get_paga_graph(meso_ccs@cds)
edge_whitelist = igraph::as_data_frame(paga_graph)
edge_blacklist = igraph::as_data_frame(igraph::complementer(paga_graph))


# plot paga graph on umap -----------------------------------------------------


plot_paga_edges <- function(ccm, edge_whitelist) {
  
  
  plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
  plot_df$cell = row.names(plot_df)
  
  plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
  plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]
  
  umap_centers = centroids(ccm@ccs)
  
  edge_whitelist = edge_whitelist %>% 
    left_join(umap_centers, by=c("from" = "cell_group")) %>% 
    left_join(umap_centers, by=c("to" = "cell_group"))
  
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
    theme(legend.position = "none") +
    monocle3:::monocle_theme_opts() + 
    geom_segment(data = edge_whitelist, 
                    aes(x = V1.x, 
                        xend= V1.y, 
                        y= V2.x, 
                        yend = V2.y))
  
  return(gp)
  
}


get_umap_coords <- function(df,umap_centers) {
  df %>% 
    left_join(umap_centers, by = c("from"="cell_group")) %>% rename("umap_from_1"=V1, "umap_from_2"=V2, "umap_from_3"=V3) %>% 
    left_join(umap_centers, by = c("to"="cell_group")) %>% rename("umap_to_1"=V1, "umap_to_2"=V2, "umap_to_3"=V3)
  
}

# distance penalties 

dist_matrix = as.matrix(dist(umap_centers, method = "euclidean", upper=T, diag = T))
dist_df = dist_matrix %>% 
  as.data.frame %>% 
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to") %>% 
  get_umap_coords(umap_centers)

init_pm <- init_penalty_matrix(meso_ccs)

init_penalty_df = init_pm %>% 
  as.data.frame %>% 
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to") %>% 
  get_umap_coords(umap_centers)

paga_pm <- init_penalty_matrix(meso_ccs, whitelist = edge_whitelist)

paga_penalty_df = paga_pm %>% 
  as.data.frame %>% 
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to") %>% 
  get_umap_coords(umap_centers)

paga_base = plot_paga_edges(meso_ccm, edge_whitelist)

paga_base + 
  geom_segment(data = paga_penalty_df, aes(x = umap_to_1,
    y = umap_to_2,
    xend = umap_from_1,
    yend = umap_from_2,
    color = value)) + scale_color_viridis_c()


paga_pm_bl<- init_penalty_matrix(meso_ccs, whitelist = edge_whitelist, blacklist = edge_blacklist)

paga_penalty_bl_df = paga_pm_bl %>% 
  as.data.frame %>% 
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to") %>% 
  get_umap_coords(umap_centers)



dist_df %>% rename("dist"=value) %>% 
  left_join(init_penalty_df %>% select(from, to, "init"=value), by = c("from", "to")) %>% 
  left_join(paga_penalty_df %>% select(from, to, "paga"=value), by = c("from", "to")) %>%
  # left_join(paga_penalty_bl_df%>% select(from, to, "paga_bl"=value),  by = c("from", "to")) %>%
  pivot_longer(-c(from, to, dist, umap_from_1, umap_from_2, umap_from_3, umap_to_1, umap_to_2, umap_to_3)) %>% 
  ggplot(aes(dist,value)) + geom_point() + facet_wrap(~name) + monocle3:::monocle_theme_opts() 

# compare the distribution of distance penalties to actual distances

dist_df %>% rename("dist"=value) %>% 
  left_join(init_penalty_df %>% select(from, to, "init"=value), by = c("from", "to")) %>% 
  left_join(paga_penalty_df %>% select(from, to, "paga"=value), by = c("from", "to")) %>% 
  pivot_longer(-c(from, to, umap_from_1, umap_from_2, umap_from_3, umap_to_1, umap_to_2, umap_to_3)) %>% 
  ggplot(aes(value)) + geom_density() + facet_wrap(~name, scales = "free") + monocle3:::monocle_theme_opts() 


# -------------------------------------------------------------------------------


meso_ccm_paga  = new_cell_count_model(meso_ccs,
                                      model_formula_str = "~ as.factor(timepoint)",
                                      whitelist=edge_whitelist,
                                      base_penalty=25)

plot(meso_ccm_paga@best_model, output = "corrplot")

time_18_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="18"))
time_20_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="20"))
time_22_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="22"))
time_24_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="24"))


# large neg combo -------------------------------------------------------------

cond_18_vs_24_tbl_paga = compare_abundances(meso_ccm_paga, time_18_paga, time_24_paga)
plot_contrast(meso_ccm_paga, cond_18_vs_24_tbl_paga)

# all combos in between ---------------------------------------------------------

cond_18_vs_20_tbl_paga = compare_abundances(meso_ccm_paga, time_18_paga, time_20_paga)
plot_contrast(meso_ccm_paga, cond_18_vs_20_tbl_paga)

cond_18_vs_22_tbl_paga = compare_abundances(meso_ccm_paga, time_18_paga, time_22_paga)
plot_contrast(meso_ccm_paga, cond_18_vs_22_tbl_paga)

cond_20_vs_22_tbl_paga = compare_abundances(meso_ccm_paga, time_20_paga, time_22_paga)
plot_contrast(meso_ccm_paga, cond_20_vs_22_tbl_paga)

cond_20_vs_24_tbl_paga = compare_abundances(meso_ccm_paga, time_20_paga, time_24_paga)
plot_contrast(meso_ccm_paga, cond_20_vs_24_tbl_paga)

cond_22_vs_24_tbl_paga = compare_abundances(meso_ccm_paga, time_22_paga, time_24_paga)
plot_contrast(meso_ccm_paga, cond_22_vs_24_tbl_paga)

log_abundance_thresh = -5
corr_edge_coords_umap_delta_abund_18_24_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                  umap_centers,
                                                                  cond_18_vs_24_tbl_paga,
                                                                  log_abundance_thresh)


corr_edge_coords_umap_delta_abund_18_20_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                  umap_centers,
                                                                  cond_18_vs_20_tbl_paga,
                                                                  log_abundance_thresh)

corr_edge_coords_umap_delta_abund_18_22_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                       umap_centers,
                                                                       cond_18_vs_22_tbl_paga,
                                                                       log_abundance_thresh)


corr_edge_coords_umap_delta_abund_20_22_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                       umap_centers,
                                                                       cond_20_vs_22_tbl_paga,
                                                                       log_abundance_thresh)

corr_edge_coords_umap_delta_abund_20_24_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                       umap_centers,
                                                                       cond_20_vs_24_tbl_paga,
                                                                       log_abundance_thresh)

# pos_edge_df = corr_edge_coords_umap_delta_abund_18_20_paga %>% 
#   filter(pcor > 0) 
# 
# pos_edge_df %>% arrange(-pcor) %>% head(5)
# calc_max_flow(pos_edge_df, source = "8", target = "9") %>% 
#   # filter(flow > 0.1) %>% 
#   plot_edges() 
  


# try new weighting scheme -----------------------------------------------------

# weighting scheme
# three sources of information: 
# 1. whether the node change the same way in the contrast
# 2. wehether they are correlated across indiviuals
# 3. how similar their transcriptomes are



time_18_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="18"))
time_24_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="24"))

cond_18_vs_24_tbl_paga = compare_abundances(meso_ccm_paga, time_18_paga, time_24_paga)
pc_plot = plot_contrast(meso_ccm_paga, cond_18_vs_24_tbl_paga)
ggsave(pc_plot, filename = "plot_contrast_18vs24_paga_base25.png")

# let's focus on path from cluster 4 to cluster 9 for now ---------------------

corr_edge_coords_umap_delta_abund_18_24_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                       umap_centers,
                                                                       cond_18_vs_24_tbl_paga,
                                                                       log_abundance_thresh=-5)

neg_edges = corr_edge_coords_umap_delta_abund_18_24_paga %>% 
  filter(edge_type!="undirected")

edge_4_to_9_df =filter(corr_edge_coords_umap_delta_abund_18_24_paga,from == "4", to == "")

# plot this only 
bp_18_24_paga = return_baseplot(meso_ccm_paga,cond_18_vs_24_tbl_paga)

# edge_plot = bp_18_24_paga + 
#   geom_segment(data = edge_4_to_9_df,
#   aes(x = umap_from_1,
#       y = umap_from_2,
#       xend=(umap_to_1+umap_from_1)/2,
#       yend = (umap_to_2+umap_from_2)/2),
#   color="black",
#   linejoin='mitre',size=2,
#   arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
#   geom_segment(data=edge_4_to_9_df,
#     aes(x = umap_to_1,
#         y = umap_to_2,
#         xend = umap_from_1,
#         yend = umap_from_2), 
#     color="black", size=2) 
# ggsave(edge_plot, filename = "examples/neg_corr_cluster4_to_cluster9.png")


# plot the pos edges only --------------------------------------------------

pos_edge_df = corr_edge_coords_umap_delta_abund_18_24_paga %>% filter(pcor > 0) 

pos_edge_plot = bp_18_24_paga + 
  geom_segment(data=pos_edge_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend = umap_from_1,
                   yend = umap_from_2), 
               color="gray")

# ggsave(pos_edge_plot, filename = "examples/pos_corr_cluster4_to_cluster9.png")
ggsave(pos_edge_plot, filename = "examples/pos_corr_edges_basepenalty_25.png")


# now let's weight the edges by FC
# delta_log_abund_df = cond_18_vs_24_tbl_paga %>% select(cell_group, delta_log_abund)
# get adjacent clusters from edges

# dist_plot = ggplot()+
#   geom_segment(data = edge_4_to_9_df,
#                aes(x = umap_from_1,
#                    y = umap_from_2,
#                    xend=(umap_to_1+umap_from_1)/2,
#                    yend = (umap_to_2+umap_from_2)/2),
#                color="black",
#                linejoin='mitre', size=2,
#                arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))+
#   geom_segment(data=edge_4_to_9_df,
#                aes(x = umap_to_1,
#                    y = umap_to_2,
#                    xend = umap_from_1,
#                    yend = umap_from_2), 
#                color="black", size=2) +
#   geom_segment(data=pos_edge_df,
#                aes(x = umap_to_1,
#                    y = umap_to_2,
#                    xend = umap_from_1,
#                    yend = umap_from_2, 
#                    color=z_z),size=2) + monocle3:::monocle_theme_opts() + 
#   scale_color_gradient2(
#     low = "#122985",
#     mid = "white",
#     high = "green4",
#     na.value = "white"
#   ) 
# 
# ggsave(total_weight, filename = "examples/total_weight_18v24.png")
# ggsave(fc_plot, filename = "examples/foldchange_weight_18v24.png")
# ggsave(pcor_plot, filename = "examples/pcor_weight_18v24.png")
# ggsave(dist_plot, filename = "examples/dist_weight_18v24.png")


calc_shortest_path <- function(edges, source, target) {
  edges = edges %>% dplyr::select(from, to, x, y, z, x_z, y_z, z_z, weight)
  G <- igraph::as.directed(igraph::graph_from_data_frame(edges, directed=FALSE))
  
  mf = igraph::shortest_paths(G, from = source, to = target, weights = igraph::E(G)$weight, output="epath")
  directed_subgraph = igraph::subgraph.edges(G, mf$epath[[1]])
  
  # turn back into a dataframe
  new_pos_edge_coords_df = igraph::as_data_frame(directed_subgraph)
  
  return(new_pos_edge_coords_df)
}


get_pos_pcor_edges <- function(umap_centers, 
                               corr_edge_coords_umap_delta_abund, 
                               alpha= 1, 
                               beta = 1, 
                               gamma = 1, 
                               sum_weights = F) {
  
  row.names(umap_centers) <- umap_centers$cell_group
  
  dist_df = dist(umap_centers[,-1], method = "euclidean", upper=T, diag = T) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("from") %>% 
    pivot_longer(-from, names_to = "to", values_to = "z")
  
  
  
  pos_edge_df = corr_edge_coords_umap_delta_abund %>% 
    filter(pcor > 0)  %>% 
    left_join(dist_df, by=c("from", "to")) %>% 
    mutate(x = abs(from_delta_log_abund-to_delta_log_abund), 
           y = 1/pcor) %>% 
    mutate(x_z = (x-min(x))/(max(x)-min(x)), 
           y_z = (y-min(y))/(max(y)-min(y)), 
           z_z = (z-min(z))/(max(z)-min(z))) 
  
  if (sum_weights) {
    pos_edge_df = pos_edge_df %>% mutate(weight = alpha*x_z + beta*y_z + gamma*z_z)
  } else {
    pos_edge_df = pos_edge_df %>% mutate(weight = alpha*x_z * beta*y_z * gamma*z_z)
  }
  
  return(pos_edge_df)
  
}


plot_shortest_path <- function(ccm, 
                               cond_a_vs_b_tbl, 
                               source, 
                               target,
                               alpha=1, 
                               beta=1,
                               gamma=1,
                               log_abundance_thresh = -5, 
                               plot_sp = T, 
                               plot_weights = T, 
                               sum_weights=T, 
                               edge_size=2) {
  
  umap_centers = centroids(ccm@ccs)
  corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(ccm,
                                                              umap_centers,
                                                              cond_a_vs_b_tbl,
                                                              log_abundance_thresh)
  
  bp = return_baseplot(ccm, 
                       cond_a_vs_b_tbl)
  

  target_edge_df = filter(corr_edge_coords_umap_delta_abund,
                         from == source, 
                         to == target)
  
  pos_edge_df = get_pos_pcor_edges(umap_centers, 
                                   corr_edge_coords_umap_delta_abund, 
                                   alpha = alpha, 
                                   beta = beta, 
                                   gamma = gamma, 
                                   sum_weights=sum_weights)
  
  path_df = calc_shortest_path(pos_edge_df, source, target) %>% select(from,to,weight) %>% 
    get_umap_coords(umap_centers)

  gp <- bp +
      geom_segment(data = target_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = umap_to_1,
                       y = umap_to_2,
                       xend=umap_from_1,
                       yend = umap_from_2),
                   color="black",size=edge_size ) +
      geom_segment(data = target_edge_df %>% dplyr::filter(edge_type == "directed_to_from"),
                   aes(x = umap_to_1,
                       y = umap_to_2,
                       xend=(umap_to_1+umap_from_1)/2,
                       yend = (umap_to_2+umap_from_2)/2,
                       ),size=edge_size,
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) + 
      geom_segment(data = target_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                      aes(x = umap_from_1,
                          y = umap_from_2,
                          xend=umap_to_1,
                          yend = umap_to_2,
                          ),size=edge_size ,
                      color="black") +
      geom_segment(data = target_edge_df %>% dplyr::filter(edge_type == "directed_from_to"),
                   aes(x = umap_from_1,
                       y = umap_from_2,
                       xend=(umap_from_1+umap_to_1)/2,
                       yend = (umap_from_2+umap_to_2)/2),
                   size=edge_size,
                   color="black",
                   linejoin='mitre',
                   arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))
  

  if (plot_sp) {
    gp = gp + geom_segment(data = path_df,
                           aes(x = umap_to_1,
                               y = umap_to_2,
                               xend = umap_from_1,
                               yend = umap_from_2),
                           color = "green4", 
                           size = 2)
  } 
  if (plot_weights) {
    gp = gp + geom_segment(data=pos_edge_df,
                   aes(x = umap_to_1,
                       y = umap_to_2,
                       xend = umap_from_1,
                       yend = umap_from_2,
                       color = weight))
  }
     
  return(gp)
}



# with base penalty 25
# dir.create("examples/shortest_path_plots_bp_25")

lapply(1:nrow(neg_edges), function(i) {
  source = neg_edges[i,]$from
  target = neg_edges[i,]$to
  gp = plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, source, target, 
                          plot_sp = T, plot_weights = F, sum_weights = F)
  ggsave(gp,
         filename = paste0("examples/shortest_path_plots_bp_25/shortest_path_18_24_s=",source,"_t=",target, ".png"))
  
  gp = plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, source, target, 
                          plot_sp = T, plot_weights = F, sum_weights = T)
  ggsave(gp, 
         filename = paste0("examples/shortest_path_plots_bp_25/shortest_path_18_24_sumweightss=",source,"_t=",target, ".png"))
})



plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, "4", "9", 
                   plot_sp = T, plot_weights = F)
plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, "4", "9", 
                   plot_sp = T, plot_weights = F, sum_weights = T)

plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, "6", "8", 
                   plot_sp = T, plot_weights = T)

plot_shortest_path(meso_ccm_paga, cond_18_vs_24_tbl_paga, "6", "10", 
                   plot_sp = T, plot_weights = T, sum_weights = F)

# collate green edges ---------------------------------------------------------


collate_green_edges <- function(corr_edge_coords_umap_delta_abund, umap_centers) {
  
  neg_edges = corr_edge_coords_umap_delta_abund %>% filter(edge_type!="undirected")
  
  if (nrow(neg_edges) > 0) {
    green_edges = lapply(1:nrow(neg_edges), function(i) {
      source = neg_edges[i,]$from
      target = neg_edges[i,]$to
      get_pos_pcor_edges(umap_centers, 
                         corr_edge_coords_umap_delta_abund, 
                         alpha = 1, 
                         beta = 1, 
                         gamma = 1, 
                         sum_weights=T) %>% 
        calc_shortest_path(source, target) %>% select(from,to,weight)
    })
    
    green_edge_df = do.call(rbind,green_edges) %>% 
      group_by(from,to) %>% 
      summarise(weight=sum(weight), n=n()) %>% 
      get_umap_coords(umap_centers)
    return(green_edge_df)  
  }
  
}

green_edge_df = collate_green_edges(corr_edge_coords_umap_delta_abund_18_24_paga, umap_centers)

# plot just the green lines 
green_line_plot = bp_18_24_paga + 
    geom_segment(data = green_edge_df,
    aes(x = umap_from_1,
        y = umap_from_2,
        xend=umap_to_1,
        yend = umap_to_2),
    color="green4", size=2) 
ggsave(green_line_plot, filename = "examples/green_line_concat.png")

# use this graph as a whitelist


green_edge_blacklist = igraph::graph_from_data_frame(green_edge_df) %>% 
                       igraph::complementer() %>% 
                       igraph::as_data_frame()


meso_ccm_green  = new_cell_count_model(meso_ccs,
                                      model_formula_str = "~ as.factor(timepoint)",
                                      whitelist=green_edge_df,
                                      # blacklist = green_edge_blacklist,
                                      base_penalty=25)


plot(meso_ccm_green@best_model, output="corrplot")

time_18_green = estimate_abundances(meso_ccm_green, data.frame(timepoint="18"))
time_24_green = estimate_abundances(meso_ccm_green, data.frame(timepoint="24"))

cond_18_vs_24_tbl_green = compare_abundances(meso_ccm_green, time_18_green, time_24_green)
plot_contrast(meso_ccm_green, cond_18_vs_24_tbl_green)
# plot_contrast(meso_ccm_paga, cond_18_vs_24_tbl_paga)

# iterate over multiple contrasts ---------------------------------------------

timepoints = unique(colData(meso_18_24_cds)$timepoint) %>% sort()
timepoint_combos = combn(timepoints, 2) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("t1"=V1, "t2"=V2)
i = 1


greenedges = lapply(1:nrow(timepoint_combos), function(i) {
  t1 = timepoint_combos[i,]$t1
  t2 = timepoint_combos[i,]$t2
  
  time_t1 = estimate_abundances(meso_ccm_paga, data.frame(timepoint=t1))
  time_t2 = estimate_abundances(meso_ccm_paga, data.frame(timepoint=t2))
  
  cond_t1_vs_t2_tbl = compare_abundances(meso_ccm_paga, time_t1, time_t2)
  
  corr_edge_coords_umap_delta_abund = collect_pln_graph_edges(meso_ccm_paga,
                                                              umap_centers,
                                                              cond_t1_vs_t2_tbl,
                                                              log_abundance_thresh=-5)
  
  green_edge_df = collate_green_edges(corr_edge_coords_umap_delta_abund, 
                                      umap_centers) %>% mutate(timepoint_x=t1, timepoint_y=t2)
  
  green_edge_df
  
})
greenedges_time_df = do.call(rbind,greenedges)

total_green_edges = greenedges_time_df %>% 
  group_by(from,to) %>% 
  summarise(total_weight = sum(weight), total_n = sum(n)) %>% 
  get_umap_coords(umap_centers)


comb_green_edge_plot = bp_18_24_paga + 
  geom_segment(data = total_green_edges, 
               aes(x = umap_from_1,
                   y = umap_from_2,
                   xend=umap_to_1,
                   yend = umap_to_2, 
                   size = total_weight),
               color="green4") 

# ggsave(comb_green_edge_plot, filename="examples/muscle_time_plots/comb_green_edge_concat.png")

# add comb green in whitelist 

total_green_edge_blacklist = igraph::graph_from_data_frame(total_green_edges) %>% 
  igraph::complementer() %>% 
  igraph::as_data_frame()


meso_ccm_green_total = new_cell_count_model(meso_ccs,
                                       model_formula_str = "~ as.factor(timepoint)",
                                       whitelist = total_green_edges,
                                       # blacklist = total_green_edge_blacklist,
                                       base_penalty=25)

plot(meso_ccm_green_total@best_model, output="corrplot")

time_18_green_total = estimate_abundances(meso_ccm_green_total, data.frame(timepoint="18"))
time_24_green_total = estimate_abundances(meso_ccm_green_total, data.frame(timepoint="24"))

cond_18_vs_24_tbl_green_total = compare_abundances(meso_ccm_green_total, 
                                                   time_18_green_total, time_24_green_total)

plot_contrast(meso_ccm_paga, cond_18_vs_24_tbl_paga)
plot_contrast(meso_ccm_green, cond_18_vs_24_tbl_green)
plot_contrast(meso_ccm_green_total, cond_18_vs_24_tbl_green_total)

# from to       pcor from_delta_log_abund to_delta_log_abund umap_from_1 umap_from_2  umap_to_1 umap_to_2
# 1    1  5 -0.2684833           -2.0100033          2.2373482    4.149991  -3.0304539  0.7074403  6.825218
# 2    2 10 -0.1371792           -2.0153915          0.3293776   -7.027485  -0.3268598 -2.3459778  4.619233
# 3    8 10 -0.0154430           -2.1464688          0.3293776   -6.322590   0.2195906 -2.3459778  4.619233
# 4    9 10 -0.0477801           -0.9478383          0.3293776    2.967358  -3.8445025 -2.3459778  4.619233
# 5    9  5 


