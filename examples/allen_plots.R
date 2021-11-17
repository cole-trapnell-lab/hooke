library(monocle3)
library(hooke)
setwd("~/OneDrive/UW/Trapnell/hooke/")
library(devtools)
load_all("~/OneDrive/UW/Trapnell/bin/monocle3-dev/")

dir.create("examples/allen_plots")

meso_18_24_cds <- readRDS("../gap-notebook-ct-1/my_R_objects/meso_18_24_cds.rds")

plot_cells(meso_18_24_cds, color_cells_by = "cell_type_broad")

colData(meso_18_24_cds)$new_cluster = paste0("cluster_", 
                                             colData(meso_18_24_cds)$cluster)

meso_ccs = new_cell_count_set(meso_18_24_cds,
                              sample_group = "embryo",
                              cell_group = "cluster")


# -----------------------------------------------------------------------------

color_by_colname <- function(cds, 
                             colname = "new_cluster",
                               cell_size=1,
                               legend_position="none") {
  
  plot_df= colData(cds) %>% as.data.frame()
  plot_df$umap2D_1 <- reducedDim(cds, type="UMAP")[,1]
  plot_df$umap2D_2 <- reducedDim(cds, type="UMAP")[,2]
  
  plot_df$cell_group = plot_df[[colname]]
  
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
      data = plot_df, 
      aes(umap2D_1, umap2D_2, 
          color=cell_group),
      size = cell_size,
      stroke = 0
    ) +
    theme_void() +
    theme(legend.position = "bottom") + 
    monocle3:::monocle_theme_opts()
  return(gp)
}

# undebug(color_by_colname)
color_by_colname(meso_18_24_cds, "cell_type_broad")

# color_by_colname(meso_18_24_cds, "timepoint")

# right now only works for by cluster
paga_graph = hooke:::get_paga_graph(meso_ccs@cds)
edge_whitelist = igraph::as_data_frame(paga_graph)
# edge_whitelist$from = paste0("cluster_",edge_whitelist$from )
# edge_whitelist$to = paste0("cluster_",edge_whitelist$to )

meso_ccm_paga  = new_cell_count_model(meso_ccs,
                                      model_formula_str = "~ as.factor(timepoint)",
                                      whitelist=edge_whitelist,
                                      base_penalty=30)

time_18_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="18"))
time_24_paga = estimate_abundances(meso_ccm_paga, data.frame(timepoint="24"))

cond_18_vs_24_tbl_paga = compare_abundances(meso_ccm_paga, time_18_paga, time_24_paga)
plot_contrast(meso_ccm_paga, cond_18_vs_24_tbl_paga)


bp = return_baseplot(meso_ccm_paga, 
         cond_18_vs_24_tbl_paga, legend_position="bottom")
bp

#ggsave(bp, filename="examples/allen_plots/meso_ref_18v24.png", width = 10, height=7)




# get green edges --------------------------------------------------------------


umap_centers = centroids(meso_ccs)

corr_edge_coords_umap_delta_abund_18_24_paga = collect_pln_graph_edges(meso_ccm_paga,
                                                                       umap_centers,
                                                                       cond_18_vs_24_tbl_paga,
                                                                       log_abundance_thresh=-5)


neg_edges = corr_edge_coords_umap_delta_abund_18_24_paga %>% filter(edge_type!="undirected")


green_edges = lapply(1:nrow(neg_edges), function(i) {
  source = neg_edges[i,]$from
  target = neg_edges[i,]$to
  pos_edges = get_pos_pcor_edges(umap_centers, 
                     corr_edge_coords_umap_delta_abund_18_24_paga, 
                     alpha = 1, 
                     beta = 1, 
                     gamma = 1, 
                     sum_weights=T) 
  pos_edges %>% calc_shortest_path(source, target) %>% select(from,to,weight)
})

green_edge_df = collate_green_edges(corr_edge_coords_umap_delta_abund_18_24_paga, umap_centers)

bp_18_24_paga = return_baseplot(meso_ccm_paga,cond_18_vs_24_tbl_paga)

green_line_plot = bp_18_24_paga +
  geom_segment(data = green_edge_df,
               aes(x = umap_to_1,
                   y = umap_to_2,
                   xend=umap_from_1,
                   yend = umap_from_2),
               color="black",size=2 ) + 
  geom_segment(data = green_edge_df,
               aes(x = umap_from_1,
                   y = umap_from_2,
                   xend=(umap_to_1+umap_from_1)/2,
                   yend = (umap_to_2+umap_from_2)/2),
                   size=2,
               color="black",
               linejoin='mitre',
               arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))

green_line_plot
ggsave(green_line_plot, filename="examples/allen_plots/wt_blackine_plot.png",width = 10, height=7)


# ------------------------------------------------------------------------------

# use model on the projected mutant -------------------------------------------
# use previous graph as whitelist prior 

# meso_mt_cds <- readRDS("../gap-notebook-ct-1/my_R_objects/meso_tbx16-like_18_24_proj_cds.rds")


plot_cells(meso_mt_cds, color_cells_by = "new_cluster")

colData(meso_mt_cds)$x = reducedDims(meso_mt_cds)[["UMAP"]][,1]
colData(meso_mt_cds)$y= reducedDims(meso_mt_cds)[["UMAP"]][,2]

colData(meso_mt_cds) %>% as.data.frame() %>% 
  ggplot(aes(x,y, color=new_cluster)) + geom_point()



plot_cells(meso_18_24_cds, color_cells_by = "new_cluster")

meso_mt_ccs = new_cell_count_set(meso_mt_cds,
                                 sample_group = "embryo",
                                 cell_group = "new_cluster")

green_edge_blacklist = igraph::graph_from_data_frame(green_edge_df) %>% 
  igraph::complementer() %>% 
  igraph::as_data_frame()

# whitelist_edges = green_edge_df
whitelist_edges = rbind(green_edge_df %>% select(from,to),
      neg_edges %>% select(from,to))

meso_mt_ccm = new_cell_count_model(meso_mt_ccs,
                                   model_formula_str = "~ as.factor(timepoint) + gene_target",
                                   whitelist = whitelist_edges,
                                   blacklist = green_edge_blacklist,
                                   base_penalty = 30)

plot(meso_mt_ccm@best_model, output="corrplot")


# ------------------------------------------------------------------------------

model_coef = coef(meso_mt_ccm@best_model)
gene_target_coefs = as.data.frame(model_coef)[["gene_targettbx16-msgn1"]]
names(gene_target_coefs) = rownames(as.data.frame(model_coef))


# wt pcor 

wt_pcor_matrix = as.matrix(meso_ccm_paga@best_model$latent_network(type="partial_cor"))
gene_target_coefs_updated = gene_target_coefs

for (i in 1:nrow(green_edge_df)) {
  
  x = paste0("cluster_",green_edge_df[i,]$from)
  y = paste0("cluster_",green_edge_df[i,]$to)
  gene_target_coefs_updated[y] = gene_target_coefs_updated[y] - gene_target_coefs[x]
  # gene_target_coefs_updated[y] = gene_target_coefs_updated[y] - gene_target_coefs[x]*wt_pcor_matrix[x,y]

}

pcor_weights = left_join(whitelist_edges,
          wt_pcor_matrix %>% as.data.frame %>% rownames_to_column() %>% pivot_longer(-rowname),
          by= c("from"="rowname", "to"="name")) %>%
  select(from, to, value)

pcor_weights

mt_umap_centers = centroids(meso_mt_ccs)


# make green edges

  
#-------------------------------------------------------------------------------

color_by_residuals <- function(ccm, 
                               coefficients,
                               cell_size=1,
                               legend_position="none") {
    
    plot_df = ccm@ccs@metadata[["cell_group_assignments"]] %>% dplyr::select(cell_group)
    plot_df$cell = row.names(plot_df)
    
    plot_df$umap2D_1 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,1]
    plot_df$umap2D_2 <- reducedDim(ccm@ccs@cds, type="UMAP")[plot_df$cell,2]
    
    res_df = as.data.frame(coefficients) %>% rownames_to_column("cell_group")
    colnames(res_df)[2] = "value"
    
    plot_df = plot_df %>% left_join(res_df, by="cell_group")

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
        data = plot_df, 
        aes(umap2D_1, umap2D_2, 
            color=value),
        size = cell_size,
        stroke = 0
      ) +
      scale_color_gradient2(
        low = "#122985",
        mid = "white",
        high = "red4",
        na.value = "white",
        limits = c(-3,3)
      ) +
      theme_void() +
      theme(legend.position = legend_position) +
      monocle3:::monocle_theme_opts()
    return(gp)
}


color_by_residuals(meso_mt_ccm, gene_target_coefs) %>% 
  ggsave(filename="examples/allen_plots/mt_before_plot.png", width = 10, height=7)
color_by_residuals(meso_mt_ccm, gene_target_coefs_updated) %>% 
  ggsave(filename="examples/allen_plots/mt_after_plot.png",width = 10, height=7)
# color_by_residuals(meso_mt_ccm, gene_target_coefs_updated-gene_target_coefs)
# color_by_residuals(meso_mt_ccm, (gene_target_coefs-gene_target_coefs_updated)/gene_target_coefs)

# what do these look like compared to beta binomial 

time_24_mt = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target ="tbx16"))
time_24_wt = estimate_abundances(meso_mt_ccm, data.frame(timepoint="24", gene_target ="ctrl-inj"))

cond_24_wt_vs_mt_tbl = compare_abundances(meso_mt_ccm, time_24_wt, time_24_mt)
return_baseplot(meso_mt_ccm, 
                cond_24_wt_vs_mt_tbl)





# ITERATION 3 -----------------------------------------------------------------

# make a combo one

meso_comb_cds <- combine_cds(list(meso_18_24_cds, meso_mt_cds), 
                             keep_reduced_dims = T, cell_names_unique = T)

colData(meso_comb_cds)$sample = NULL
colData(meso_comb_cds)$cluster_transfer = NULL
colData(meso_comb_cds)$cluster2 = NULL
colData(meso_comb_cds)$umap_1 = NULL
colData(meso_comb_cds)$umap_2 = NULL
colData(meso_comb_cds)$group_cluster = NULL

# undebug(new_cell_count_set)
meso_comb_ccs = new_cell_count_set(meso_comb_cds,
                              sample_group = "embryo",
                              cell_group = "new_cluster")

meso_comb_ccm = new_cell_count_model(meso_comb_ccs,
                                   model_formula_str = "~ as.factor(timepoint) + gene_target",
                                   whitelist = whitelist_edges,
                                   base_penalty = 40)

plot(meso_comb_ccm@best_model, output="corrplot")



