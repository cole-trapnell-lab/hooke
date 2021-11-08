library(monocle3)
library(hooke)

meso_cds = readRDS("/Users/coletrap/google_drive/Projects/SDG/GAPFISH/data/gap16/analysis/R_objects/gap16_mesoderm_tbx16_18-24h_sub-cds.RDS")


# assign best celltype column and reduce dims
colData(meso_cds)$cell_type = colData(meso_cds)$cell_type_broad
colData(meso_cds)$timepoint = stringr::str_remove(colData(meso_cds)$timepoint, "h")

#meso_cds = reduce_dimension(meso_cds, max_components = 2, preprocess_method = 'PCA')
meso_cds = cluster_cells(meso_cds, random_seed = 42, resolution = 1e-4)

# add attributes to colData
colData(meso_cds)$cluster = monocle3::clusters(meso_cds)
colData(meso_cds)$partition = partitions(meso_cds)
colData(meso_cds)$Size_Factor = size_factors(meso_cds)

plot_cells(meso_cds, color_cells_by = "cluster")

colData(meso_cds)$sample = NULL

#plot_cells(cds)

ccs = new_cell_count_set(meso_cds,
                         sample_group = "Oligo",
                         cell_group = "cluster")


ccm  = new_cell_count_model(ccs,
                            model_formula_str = "~timepoint + gene_target")


cond_a = estimate_abundances(ccm, data.frame(timepoint="18", gene_target="ctrl-inj"))
cond_b = estimate_abundances(ccm, data.frame(timepoint="24", gene_target="ctrl-inj"))
#plot_abundance_shift(meso_cds, best_model_umap, umap_centers, time_a, time_b, log_abundance_thresh=-1, scale_shifts_by="receiver", edge_size=2)

cond_b_vs_a_tbl = compare_abundances(ccm, cond_a, cond_b)

plot_contrast(ccm, cond_b_vs_a_tbl)


# Now try it with a white and black list:
paga_graph = hooke:::get_paga_graph(ccs@cds)
edge_whitelist = igraph::as_data_frame(paga_graph)
edge_blacklist = igraph::as_data_frame(igraph::complementer(paga_graph))

ccm  = new_cell_count_model(ccs,
                            model_formula_str = "~timepoint + gene_target",
                            whitelist=edge_whitelist,
                            base_penalty=1)


cond_a = estimate_abundances(ccm, tibble(timepoint="18", gene_target="ctrl-inj"))
cond_b = estimate_abundances(ccm, tibble(timepoint="24", gene_target="ctrl-inj"))
#plot_abundance_shift(meso_cds, best_model_umap, umap_centers, time_a, time_b, log_abundance_thresh=-1, scale_shifts_by="receiver", edge_size=2)

cond_b_vs_a_tbl = compare_abundances(ccm, cond_a, cond_b)

plot_contrast(ccm, cond_b_vs_a_tbl)


cond_a = estimate_abundances(ccm, tibble(timepoint="24", gene_target="ctrl-inj"))
cond_b = estimate_abundances(ccm, tibble(timepoint="24", gene_target="tbx16"))
#plot_abundance_shift(meso_cds, best_model_umap, umap_centers, time_a, time_b, log_abundance_thresh=-1, scale_shifts_by="receiver", edge_size=2)

cond_b_vs_a_tbl = compare_abundances(ccm, cond_a, cond_b)

plot_contrast(ccm, cond_b_vs_a_tbl)

#hooke:::get_paga_graph(ccm)
