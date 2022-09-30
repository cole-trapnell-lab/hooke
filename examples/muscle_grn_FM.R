library(monocle3)
library(msigdbr)
library(igraph)
library(tidygraph)
library(network)
library(ggnetwork)
library(tibble)
# library(hooke)

devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")

setwd("~/OneDrive/UW/Trapnell/hooke/")

if (!file.exists("examples/R_objects/meso_18_24_cds.rd")) {
  drive_download("https://drive.google.com/file/d/1QMeWa2Pzg9coQNaybTMELwFH6MzR58uo/")
}

if (!file.exists("examples/R_objects/meso_tbx16-like_18_24_proj_cds.rds")) {
  drive_download("https://drive.google.com/file/d/1P-W-HXD6f7BhAqD5W-53DcHYbh4pLeuS/")
}



# mutant ----------------------------------------------------------------------

meso_mt_cds <- readRDS("examples/R_objects/meso_tbx16-like_18_24_proj_cds.rds")
meso_mt_cds = estimate_size_factors(meso_mt_cds)
meso_mt_cds = detect_genes(meso_mt_cds)

colData(meso_mt_cds)$cluster = gsub("cluster_", "", colData(meso_mt_cds)$new_cluster) %>% as.numeric()
names(colData(meso_mt_cds)$cluster) = rownames(colData(meso_mt_cds))
meso_mt_ccs = new_cell_count_set(meso_mt_cds,
                                 sample_group = "embryo",
                                 cell_group = "cluster")
meso_mt_pb_cds = pseudobulk(meso_mt_ccs)

# this is from muscle_grn.R
weighted_edges <- readRDS("examples/muscle_grn/weighted_edges.rds")

edges = data.frame("from" = c(4, 4), 
                   "to" = c(8, 9)) %>% 
  mutate(branch = row_number()) %>% 
  dplyr::mutate(shortest_path = purrr::map2(.f = get_shortest_path, 
                                            .x = from, .y = to, 
                                            weighted_edges)) %>% 
  dplyr::rename("root" = from, "leaf" = to) %>%
  select(root, leaf, branch, shortest_path) %>% 
  dplyr::mutate(full_path = purrr::map(.f = get_path_order, 
                                       .x = shortest_path)) %>% 
  tidyr::unnest(shortest_path) %>% 
  arrange(branch, distance_from_root) %>% 
  mutate(sub_states = purrr::map2(.f = select_states,
                                  .x = full_path,
                                  .y = from)) 


fm_path = edges %>% filter(branch == 2) 
fm_states = union(fm_path$from, fm_path$to)

knockouts = c("ENSDARG00000007329", "ENSDARG00000070546")
allowed_reg_ids = get_allowed_regulators(meso_mt_pb_cds, knockouts) 
fm_path_gene_mod = get_gene_modules(meso_mt_pb_cds, unique(unlist(fm_path$degs)), allowed_reg_ids)



fm_path = fm_path %>% 
  dplyr::mutate(states = purrr::map2(.f = purrr::possibly(
    collect_transition_states, NA_real_), 
    .x = from,
    .y = to )) 

fm_path = fm_path %>% 
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_), 
                                  .x = states, 
                                  model_formula_str = "~ as.factor(cell_group)", 
                                  meso_mt_pb_cds)) 

fm_path = fm_path %>% 
  dplyr::mutate(gene_count_model = purrr::map2(.f = purrr::possibly(
    build_pln_model_on_genes, NA_real_),
    .x = degs,
    .y = states,
    cds = meso_mt_pb_cds,
    knockout_genes = knockouts, 
    model_formula_str = "~cell_state",
    regulatory_genes = allowed_reg_ids,
    gene_module_df = fm_path_gene_mod)) 

# score these regulators 
fm_path = fm_path %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(rank_regulators, NA_real_), 
                                         .x = states,
                                         .y = gene_count_model, 
                                         meso_mt_pb_cds))


fm_regulators = fm_path %>% 
  select(regulators) %>% 
  tidyr::unnest(regulators)  


fm_regulators %>% ggplot(aes(regulator_score)) + 
  geom_density() + 
  geom_vline(xintercept = 0.1) 


filtered_fm_regulators = fm_regulators %>% 
  filter(abs(regulator_score) > 0.05) %>% 
  pull(gene_id) %>% 
  unique() # 72, 250



fm_pln_model = build_pln_model_on_genes(genes = filtered_fm_regulators,
                                        states = fm_states,
                                        cds = meso_mt_pb_cds,
                                        knockout_genes = knockouts,
                                        regulatory_genes = allowed_reg_ids,
                                        model_formula_str = "~cell_state",
                                        gene_module_df = fm_path_gene_mod)

saveRDS(fm_pln_model, "examples/muscle_grn/transcription_regulators/fm_pln_model.rds")
fm_pln_model <- readRDS("examples/muscle_grn/transcription_regulators/fm_pln_model.rds")


cov_graph = fm_pln_model@best_model %>% return_igraph()
cov_edges <- igraph::as_data_frame(cov_graph, what="edges") 

cov_edges %>% ggplot(aes(weight)) + geom_density() + geom_vline(xintercept = 0.02)

cov_edges %>% dplyr::filter(abs(weight) > 0.01) %>% nrow()

  

fm_path = fm_path %>% 
  dplyr::mutate(full_regulators = purrr::map(.f = purrr::possibly(rank_regulators, NA_real_), 
                                             .x = states,
                                             gene_model_ccm = fm_pln_model, 
                                             cds = meso_mt_pb_cds)) 
  

fm_path = fm_path %>% 
  dplyr::mutate(filter_regulators = purrr::map(full_regulators, ~filter(., regulator_score!=0)))

# this predict changes
fm_path = fm_path %>% 
  dplyr::mutate(cond_b_vs_a_tbl = 
                  purrr::map2(.f = comp_abund, 
                              .x = from,
                              .y = to,
                              ccm = fm_pln_model)) 

# maybe filter
pval_threshold = 1.0
fm_path = fm_path %>% 
          dplyr::mutate(filter_cond_b_vs_a_tbl = purrr::map(cond_b_vs_a_tbl, ~filter(., delta_p_value < pval_threshold))) 

# get edges 
fm_path = fm_path %>% 
  dplyr::mutate(edges = 
                  purrr::map(.f = collect_pln_gene_edges, 
                             .x = filter_cond_b_vs_a_tbl, 
                             ccm = fm_pln_model)) 

# filter edges for larger pcor values 
fm_path = fm_path %>% 
  dplyr::mutate(filter_edges = purrr::map(edges, ~filter(., abs(pcor) > 0.05)))

# this currently takes in filter_edges
fm_path = fm_path %>%
  dplyr::mutate(dir_edges =
                  purrr::map(.f = purrr::possibly(find_all_edges, NA_real_),
                             .x = sub_states,
                             path = fm_path))


# saveRDS(fm_path, "examples/muscle_grn/transcription_regulators/fm_path.rds")
# fm_path <- readRDS("examples/muscle_grn/transcription_regulators/fm_path.rds")



# get TF timing --------------------------------------------------------------

agg_expr_data <- aggregated_expr_data(meso_mt_cds, group_cells_by = "new_cluster")


# run top_markers 

# top_markers = monocle3::top_markers(meso_mt_cds, group_cells_by = "new_cluster")
gene_id_level = agg_expr_data %>% 
  filter(cell_group %in% paste0("cluster_", fm_states)) %>%
  filter(gene_id %in% unique(fm_regulators$gene_id)) %>%
  group_by(gene_id) %>% 
  top_n(wt = mean_expression, 1)


# plot nodes according to level -----------------------------------------------


# make nodes

all_dir_edges = fm_path %>% 
  select(dir_edges) %>% 
  filter(!is.na(dir_edges)) %>% 
  tidyr::unnest(dir_edges) 

saveRDS(all_dir_edges, "examples/muscle_grn/transcription_regulators/FM_dir_edges.rds")
# all_dir_edges <- readRDS("examples/muscle_grn/transcription_regulators/FM_dir_edges.rds")

ensembl_to_shortname = as.data.frame(rowData(meso_mt_cds)) %>% select(id, gene_short_name)



plot_grn <- function(all_dir_edges, 
                     path,
                     num_levels = 2,
                     color_nodes_by = NULL,
                     arrow.gap=0.03,
                     activator_curvature = 0.0,
                     repressor_curvature = 0.1,
                     arrow_unit = 2,
                     bar_unit = .075,
                     node_size = 2) {


  # remove any edge duplicates
  ade = all_dir_edges %>%
        select(to, from, edge_type, pcor) %>% 
        group_by(from,to) %>%
        arrange(edge_type) %>% # needs to deal w this later, only ggnetwork can't handle it
        slice(1) %>% # currently chooses the activator if multiple
        ungroup() %>% 
        mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  

  G = ade %>% select(from, to, edge_type, pcor, scaled_weight)  %>% igraph::graph_from_data_frame(directed = T)

  state_order = path %>% select(to, distance_from_root) %>%
    rbind(data.frame("to"="4", distance_from_root=0)) %>%
    mutate("cell_group" = paste0("cluster_", to))

  level_df = data.frame("name" = V(G)$name) %>%
    left_join(gene_id_level, by = c("name" = "gene_id")) %>%
    left_join(state_order, by = "cell_group") %>%
    group_by(distance_from_root) %>%
    mutate(rn = row_number()) %>% 
    mutate(group = cut(rn, num_levels, labels=F)) %>% 
    mutate(group_label = as.numeric(distance_from_root) + (group-1)*(0.75/num_levels)) %>% 
    ungroup() %>%  
    tibble::column_to_rownames("name")
  
  
  regulator_score_df =
    path %>%
    select(regulators) %>%
    tidyr::unnest(regulators) %>%
    select(gene_id, regulator_score) %>%
    group_by(gene_id) %>%
    arrange(-abs(regulator_score)) %>% # take the highest available score? 
    slice(1) %>% 
    ungroup()
  

  # run sugiyama layout
  lay1 <- layout_with_sugiyama(G, layers=level_df[V(G)$name,][["group_label"]])

  g = ggnetwork(igraph::as_data_frame(G), layout = lay1$layout, arrow.gap = arrow.gap)
  
  # add level information
  g = g %>% left_join(level_df %>% rownames_to_column("id"), by = c("vertex.names"="id"))
  g = g %>% left_join(regulator_score_df, by = c("vertex.names" = "gene_id") )
  

  p <- ggplot(mapping = aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
      # draw unvdirected edges
      geom_edges(data = g %>% filter(edge_type == "undirected"),
                 alpha = 0.2) +
      # draw activator edges
      geom_edges(data = g %>% filter(edge_type == "activator"),
                 arrow = arrow(length = unit(arrow_unit, "pt"), type="closed"),
                 curvature = activator_curvature) +
      # draw repressor edges
      geom_edges(data = g %>% filter(edge_type == "repressor"),
                 arrow = arrow(angle = 90, length = unit(bar_unit, "cm")),
                 curvature = repressor_curvature)

  if (is.null(color_nodes_by) == FALSE) {
    
    # if numerical 
    if (color_nodes_by == "cell_group" & is.numeric(g[["cell_group"]])) {
      p = p + geom_nodelabel(data = g,
                             aes(fill = as.factor(get(color_nodes_by)),
                                 label = gene_short_name),
                             size = node_size) + 
        labs(fill = color_nodes_by) + 
        scale_fill_gradient2(low = "darkblue", mid = "white", high="red4")
    } else if (color_nodes_by == "regulator_score"){
      p = p + geom_nodelabel(data = g,
                         aes(fill = regulator_score,
                             label = gene_short_name, 
                             size = 3*abs(regulator_score))) + 
        scale_fill_gradient2(low = "darkblue", mid = "white", high="red4")
      
    } 
    
    else {
      # if categorical 
      p = p + geom_nodelabel(data = g,
                             aes(fill = as.factor(get(color_nodes_by)),
                                 label = gene_short_name),
                             size = node_size) + 
        labs(fill = color_nodes_by)
      
    }
    
  } else {
    p = p + geom_nodelabel(data = g,
                   aes(label = gene_short_name),
                   size = node_size)
  }

  p = p + scale_size_identity() +
      monocle3:::monocle_theme_opts() +
      theme_blank() 
  return(p)

}


all_dir_edges = all_dir_edges %>% 
  filter(edge_type != "undirected")

plot_grn(all_dir_edges, 
         path = fm_path,
         color_nodes_by = "regulator_score",
         num_levels = 3)


search_edge = function(all_dir_edges, name) {
  all_dir_edges %>% 
    left_join(ensembl_to_shortname, by = c("to"="id")) %>%
    left_join(ensembl_to_shortname, by = c("from"="id"), suffix = c(".from", ".to")) %>% 
    select(from, to, gene_short_name.to, gene_short_name.from, edge_type) %>% 
    filter(from == name | to == name |
             gene_short_name.from == name | gene_short_name.to == name) 
}

search_edge(all_dir_edges, "hey1") %>% 
  filter(edge_type != "undirected") %>% View()

# check for genes
genes_to_check = c("hey1", "myorg", "myod1", "mef2ca", "mef2cb", "smad1", "hif1aa", "hif1an", "hif1ab")

rowData(meso_mt_cds) %>% as.data.frame() %>% 
  filter(grepl("hif", gene_short_name))

# MYORG not in regulators
shortname_to_id(meso_mt_cds, genes_to_check) %in% fm_regulators$gene_id
shortname_to_id(meso_mt_cds, "hey1")


  

