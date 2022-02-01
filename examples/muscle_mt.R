
# run just on the mutant for now

# predict on mutant data -------------------------------------------------------


# gap16 <- readRDS("../gap-notebook-ct-1/R_objects/gap16_all-cells_anno_clean_2.7M_cds.RDS")
#
# meso_mt_cds <- gap16[,colData(gap16)$timepoint %in% c(18,24)]
#
# meso_mt_cds <- meso_mt_cds[,colData(meso_mt_cds)$gene_target %in% c("ctrl-inj", "tbx16", "tbx16-msgn1")]
# meso_mt_cds = meso_mt_cds[,colData(meso_mt_cds)$cell_type_broad %in% meso_cell_types]


# meso_mt_cds <- load_transform_models(meso_mt_cds, "../my_R_objects/ref-meso-18-24-models")
#
# meso_mt_cds <- preprocess_transform(meso_mt_cds, method="PCA")
# meso_mt_cds <- align_beta_transform(meso_mt_cds)
# meso_mt_cds <- reduce_dimension_transform(meso_mt_cds, method="UMAP")
#
# meso_mt_umap_dims = reducedDims(meso_mt_cds)[["UMAP"]]
# meso_res = uwot:::annoy_search(meso_mt_umap_dims,
#                               k = 10,
#                               ann = meso_mt_cds@reduce_dim_aux$UMAP$nn_index$annoy_index)
#
# cluster_labels  = get_nn_labels(meso_umap_dims,
#                                               meso_res,
#                                               as.data.frame(colData(meso_18_24_cds)),
#                                               transfer_type = "cluster")
# # colData(meso_mt_cds)$cluster= NULL
# #colData(meso_mt_cds)$cluster_transfer = cluster_labels
# colData(meso_mt_cds)$cluster = cluster_labels
# colData(meso_mt_cds)$new_cluster = paste0("cluster_",cluster_labels)
#
# plot_cells(meso_mt_cds, color_cells_by = "cell_type_broad")
# plot_cells(meso_mt_cds, color_cells_by="new_cluster")
# plot_cells(meso_mt_cds, color_cells_by="cluster_transfer")
# plot_cells(meso_18_24_cds, color_cells_by="new_cluster")
#
# meso_mt_cds = meso_mt_cds[,!is.na(colData(meso_mt_cds)$cluster)]

#
# mt_clusters = as.factor(colData(meso_mt_cds)$cluster)
# names(mt_clusters) = colData(meso_mt_cds)$cell
# colData(meso_mt_cds)$cluster = mt_clusters
# meso_mt_cds@clusters[["UMAP"]] = list(clusters = mt_clusters)
#
# meso_mt_cds@clusters
# meso_18_24_cds@clusters
#
# clusters <- factor(igraph::membership(cluster_result$optim_res))
#


# colData(meso_mt_cds)$sample = NULL
# saveRDS(meso_mt_cds, "../gap-notebook-ct-1/my_R_objects/meso_tbx16-like_18_24_cds.rds")
# saveRDS(meso_mt_cds, "../gap-notebook-ct-1/my_R_objects/meso_tbx16-like_18_24_proj_cds.rds")

meso_mt_cds <- readRDS("examples/R_objects/meso_tbx16-like_18_24_proj_cds.rds")

#meso_comb_cds <- combine_cds(list(meso_18_24_cds, meso_mt_cds), keep_reduced_dims = T, cell_names_unique = T)
#saveRDS(meso_comb_cds, "../gap-notebook-ct-1/my_R_objects/meso_comb_18_24_proj_cds.rds")

#plot_cells(meso_comb_cds, color_cells_by = "cluster")
#plot_cells(meso_comb_cds, color_cells_by = "gene_target")

# ---------------------------------------------------------------------------

# clusters(meso_mt_cds) <- NULL
meso_mt_cds = meso_mt_cds[,!is.na(colData(meso_mt_cds)$cluster)]


meso_mt_ccs = new_cell_count_set(meso_mt_cds,
                                 sample_group = "embryo",
                                 cell_group = "cluster")

meso_mt_ccm = new_cell_count_model(meso_mt_ccs,
                                   main_model_formula_str = "~ as.factor(timepoint) + gene_target",
                                   #whitelist = total_green_edges,
                                   base_penalty = 5)

plot(model(meso_mt_ccm, "reduced"), output="corrplot")

time_18_ctrlinj = estimate_abundances(meso_mt_ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
time_18_tbx16 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="18", gene_target="tbx16"))
time_24_ctrlinj = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))
time_24_tbx16 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target="tbx16"))


cond_mt_18_ctrlinj_v_tbx16_tbl = compare_abundances(meso_mt_ccm, time_18_ctrlinj, time_18_tbx16)
cond_mt_24_ctrlinj_v_tbx16_tbl = compare_abundances(meso_mt_ccm, time_24_ctrlinj, time_24_tbx16)
cond_mt_18_24_ctrlinj_tbl = compare_abundances(meso_mt_ccm, time_18_ctrlinj, time_24_ctrlinj)


plot_contrast(meso_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, scale_shifts_by = "receiver", q_value_thresh = 0.05)
plot_contrast(meso_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, scale_shifts_by = "sender", q_value_thresh = 0.05)
plot_contrast(meso_mt_ccm, cond_mt_24_ctrlinj_v_tbx16_tbl, q_value_thresh = 0.05)

colData(meso_mt_cds)$umap_1 = reducedDims(meso_mt_cds)[["UMAP"]][,1]
colData(meso_mt_cds)$umap_2 = reducedDims(meso_mt_cds)[["UMAP"]][,2]

colData(meso_mt_cds) %>%
  as.data.frame() %>%
  ggplot(aes(umap_1, umap_2, color=cluster)) + geom_point(size=0.1)

mt_umap_centers = centroids(meso_mt_ccs)
corr_edge_coords_umap_delta_abund = hooke:::collect_pln_graph_edges(meso_mt_ccm,
                                                            cond_mt_18_ctrlinj_v_tbx16_tbl,
                                                            log_abundance_thresh=-5)



sp = plot_shortest_path(meso_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, "10", "2",
                        plot_sp = T, plot_weights = F, sum_weights = T)
sp
sp = plot_shortest_path(meso_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, "10", "8",
                        plot_sp = T, plot_weights = F, sum_weights = T)

sp = plot_shortest_path(meso_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, "1", "4",
                        plot_sp = T, plot_weights = F, sum_weights = T)
sp
# ggsave(sp, filename = "examples/muscle_time_plots/shortest_path_18_ctrlinj_tbx16_cluster2_to_4.png")

# Fast muscle alone

selected_cell_types = c("fast-committed myocyte (fusing)",
                        "mature fast muscle 1",
                        #"mature fast muscle 2",
                        #"mature fast muscle 3",
                        #"mature fast muscle 4",
                        #"mature fast muscle 5",
                        #"mature fast muscle 6",
                        "fast-committed myocyte (pre-fusion)",
                        "mesodermal progenitor cells (contains PSM)",
                        #"satellite cells",
                        "myoblast")
fast_muscle_cds = meso_mt_cds[,colData(meso_mt_cds)$cell_type_sub %in% selected_cell_types]
colData(fast_muscle_cds)$timepoint_num = as.numeric(colData(fast_muscle_cds)$timepoint)

fast_muscle_cds = reduce_dimension(fast_muscle_cds)
#fast_muscle_cds = cluster_cells(fast_muscle_cds)
#plot_cells(fast_muscle_cds, color_cells_by="cluster")

plot_cells(fast_muscle_cds, color_cells_by="cell_type_sub")


fast_mt_ccs = new_cell_count_set(fast_muscle_cds,
                                 sample_group = "embryo",
                                 cell_group = "cell_type_sub")

fast_mt_ccm = new_cell_count_model(fast_mt_ccs,
                                   main_model_formula_str = "~ as.factor(timepoint) + gene_target")

fast_mt_ccm = select_model(fast_mt_ccm, criterion="EBIC", sparsity_factor=1)

plot(model(fast_mt_ccm,"reduced"), output="corrplot")

time_18_ctrlinj = estimate_abundances(fast_mt_ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
time_18_tbx16 = estimate_abundances(fast_mt_ccm, data.frame(timepoint="18", gene_target="tbx16"))
time_24_ctrlinj = estimate_abundances(fast_mt_ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))
time_24_tbx16 = estimate_abundances(fast_mt_ccm, data.frame(timepoint="24", gene_target="tbx16"))


cond_mt_18_ctrlinj_v_tbx16_tbl = compare_abundances(fast_mt_ccm, time_18_ctrlinj, time_18_tbx16)
cond_mt_24_ctrlinj_v_tbx16_tbl = compare_abundances(fast_mt_ccm, time_24_ctrlinj, time_24_tbx16)
cond_mt_18_24_ctrlinj_tbl = compare_abundances(fast_mt_ccm, time_18_ctrlinj, time_24_ctrlinj)


plot_contrast(fast_mt_ccm, cond_mt_18_24_ctrlinj_tbl, q_value_thresh = 1)

plot_contrast(fast_mt_ccm, cond_mt_18_24_ctrlinj_tbl, q_value_thresh = 1, plot_labels="none") +
  ggsave("fast_muscle_wt_18_24.png", width=7, height=6)

plot_contrast(fast_mt_ccm, cond_mt_24_ctrlinj_v_tbx16_tbl, p_value_thresh = 1, plot_labels="none") +
  ggsave("fast_muscle_tbx16_18_24.png", width=7, height=6)

plot_contrast(fast_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, scale_shifts_by = "receiver", q_value_thresh = 1)
plot_contrast(fast_mt_ccm, cond_mt_18_ctrlinj_v_tbx16_tbl, scale_shifts_by = "sender", q_value_thresh = 0.05)
plot_contrast(fast_mt_ccm, cond_mt_24_ctrlinj_v_tbx16_tbl, q_value_thresh = 1)



# mutant alone ----------------------------------------------------------------
# whitelist=total_green_edges,
# blacklist = total_green_edge_blacklist,

# mt_paga_graph = hooke:::get_paga_graph(meso_mt_ccs@cds)
# mt_edge_whitelist = igraph::as_data_frame(mt_paga_graph)
#
# meso_mt_ccm = new_cell_count_model(meso_mt_ccs,
#                                    model_formula_str = "~ as.factor(timepoint) + gene_target",
#                                    # whitelist = mt_edge_whitelist,
#                                    base_penalty = 1)
#
# plot(meso_mt_ccm@best_model, output = "corrplot")
#
#
# time_mt_18 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
# time_mt_24 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))
#
# cond_mt_18_vs_24_tbl = compare_abundances(meso_mt_ccm, time_mt_18, time_mt_24)
#
# plot_contrast(meso_mt_ccm, cond_mt_18_vs_24_tbl)
#
# # let's plot shortest path - source = 1 and target = 2
# plot_shortest_path(meso_mt_ccm, cond_mt_18_vs_24_tbl, "1", "2",
#                    plot_sp = T, plot_weights = F, sum_weights = T)
#
#
# time_mt_18 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
# time_mt_18_tbx16 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="18", gene_target="tbx16"))
#
# cond_mt_18_vs_18_tbx16_tbl = compare_abundances(meso_mt_ccm, time_mt_18, time_mt_18_tbx16)
#
# plot_contrast(meso_mt_ccm, cond_mt_18_vs_18_tbx16_tbl)
#
# sp = plot_shortest_path(meso_mt_ccm, cond_mt_18_vs_18_tbx16_tbl, "2", "4",
#                         plot_sp = T, plot_weights = F, sum_weights = T)
# ggsave(sp, filename = "examples/muscle_time_plots/shortest_path_18_ctrlinj_tbx16_cluster2_to_4.png")
#
# sp = plot_shortest_path(meso_mt_ccm, cond_mt_18_vs_18_tbx16_tbl, "2", "8",
#                         plot_sp = T, plot_weights = F, sum_weights = T)
# ggsave(sp, filename = "examples/muscle_time_plots/shortest_path_18_ctrlinj_tbx16_cluster2_to_8.png")
#
#
# time_mt_24 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))
# time_mt_24_tbx16 = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target="tbx16"))
#
# cond_mt_24_vs_24_tbx16_tbl = compare_abundances(meso_mt_ccm, time_mt_24, time_mt_24_tbx16)
#
# sp = plot_shortest_path(meso_mt_ccm, cond_mt_24_vs_24_tbx16_tbl, "2", "4",
#                         plot_sp = T, plot_weights = F, sum_weights = T)
# ggsave(sp, filename = "examples/muscle_time_plots/shortest_path_24_ctrlinj_tbx16_cluster2_to_4.png")
#
# sp = plot_shortest_path(meso_mt_ccm, cond_mt_24_vs_24_tbx16_tbl, "2", "8",
#                         plot_sp = T, plot_weights = F, sum_weights = T)
# ggsave(sp, filename = "examples/muscle_time_plots/shortest_path_24_ctrlinj_tbx16_cluster2_to_8.png")


