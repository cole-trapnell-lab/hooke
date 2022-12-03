library(hooke)

cds <- monocle3::load_a549()
colData(cds)$sample <- NULL

cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-2)
colData(cds)$cluster <- clusters(cds)
#plot_cells(cds)

ccs <- new_cell_count_set(cds,
                         sample_group = "replicate",
                         cell_group = "cluster")

ccm  <- new_cell_count_model(ccs,
                            main_model_formula_str = "~log_dose")


dose_a <- estimate_abundances(ccm, tibble(log_dose=1.000434077))
dose_b <- estimate_abundances(ccm, tibble(log_dose=2.000043427))
#plot_abundance_shift(meso_cds, best_model_umap, umap_centers, time_a, time_b, log_abundance_thresh=-1, scale_shifts_by="receiver", edge_size=2)

cond_b_vs_a_tbl <- compare_abundances(ccm, dose_a, dose_b)

plot_contrast(ccm, cond_b_vs_a_tbl)
