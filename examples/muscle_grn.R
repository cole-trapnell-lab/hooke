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

# load files ------------------------------------------------------------------

meso_cds <- readRDS("examples/R_objects/meso_18_24_cds.rds")


meso_ccs = new_cell_count_set(meso_cds,
                              sample_group = "embryo",
                              cell_group = "cluster")

paga_graph = hooke:::get_paga_graph(meso_ccs@cds)
edge_whitelist = igraph::as_data_frame(paga_graph)

meso_ccm = new_cell_count_model(meso_ccs,
                                model_formula_str = "~ as.factor(timepoint)", 
                                whitelist=edge_whitelist,
                                base_penalty=30)


time_18 = estimate_abundances(meso_ccm, data.frame(timepoint="18"))
time_20 = estimate_abundances(meso_ccm, data.frame(timepoint="20"))
time_22 = estimate_abundances(meso_ccm, data.frame(timepoint="22"))
time_24 = estimate_abundances(meso_ccm, data.frame(timepoint="24"))

# perform a contrast between 18 and 24 hours ----------------------------------

cond_18_vs_24_tbl = compare_abundances(meso_ccm, time_18, time_24)


neg_rec_edges = hooke:::collect_pln_graph_edges(meso_ccm, cond_18_vs_24_tbl) %>% 
  as_tibble %>%
  filter(edge_type != "undirected" &
           to_delta_p_value < 1 &
           from_delta_p_value < 1) 

get_positive_edges = function(ccm, cond_b_vs_a_tbl, p_value_threshold = 1.0) {
  
  positive_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
    as_tibble %>%
    filter(pcor > 0 &
             to_delta_p_value < p_value_threshold &
             from_delta_p_value < p_value_threshold)
  
}

pos_edges = get_positive_edges(meso_ccm, cond_18_vs_24_tbl)
weighted_edges = get_weighted_edges(meso_ccm, pos_edges)

saveRDS(weighted_edges, "examples/muscle_grn/weighted_edges.rds")

green_edges = neg_rec_edges %>%
  dplyr::mutate(shortest_path = purrr::map2(.f = 
    get_shortest_path, .x = from, .y = to, weighted_edges)) %>% 
  select(shortest_path) %>%
  tidyr::unnest(shortest_path) %>% 
  select(-weight) %>%
  distinct()


green_edges = green_edges %>% 
    dplyr::mutate(states = purrr::map2(.f = purrr::possibly(
    collect_transition_states, NA_real_), .x = from,
    .y = to ))

meso_pb_cds = pseudobulk(meso_ccs)

# calculate degs between WT ---------------------------------------------------
green_edges = green_edges %>% 
    dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_), 
                                    .x = states, 
                                    model_formula_str = "~ as.factor(cell_group)", 
                                    meso_pb_cds))


saveRDS(green_edges, "examples/muscle_grn/wt_pos_path.rds")
green_edges <- readRDS("examples/muscle_grn/wt_pos_path.rds")


# do a single path first ------------------------------------------------------

# pos_path = get_weighted_edges(meso_ccm, pos_edges) %>% 
#   get_shortest_path(from="4", to="9") %>% as_tibble()
# 
# pos_path = pos_path %>% 
#   dplyr::mutate(states = purrr::map2(.f = purrr::possibly(
#   collect_transition_states, NA_real_), .x = from,
#   .y = to ))
# 
# # do degs along this path
# pos_path = pos_path %>% 
#   dplyr::mutate(degs = purrr::map(.f = purrr::possibly(
#   find_degs_between_states, NA_real_), .x = states, meso_pseudobulk_cds))

# saveRDS(pos_path, "examples/muscle_grn/wt_pos_path_4_to_9.rds")
# pos_path <- readRDS("examples/muscle_grn/wt_pos_path_4_to_9.rds")


# get mutant data also --------------------------------------------------------

meso_mt_cds <- readRDS("examples/R_objects/meso_tbx16-like_18_24_proj_cds.rds")

colData(meso_mt_cds)$cluster = gsub("cluster_", "", colData(meso_mt_cds)$new_cluster) %>% as.numeric()
names(colData(meso_mt_cds)$cluster) = rownames(colData(meso_mt_cds))
meso_mt_ccs = new_cell_count_set(meso_mt_cds,
                                 sample_group = "embryo",
                                 cell_group = "cluster")

# meso_mt_ccm = new_cell_count_model(meso_mt_ccs,
#                                    model_formula_str = "~ as.factor(timepoint) + gene_target",
#                                    # whitelist = whitelist_edges,
#                                    # blacklist = green_edge_blacklist,
#                                    # base_penalty = 30
#                                    )


meso_mt_pb_cds = pseudobulk(meso_mt_ccs, gene_ids = c("tbx16", "msgn1"))

# currently use the single pos path
# mt_pos_path = pos_path %>% select(from, to, states)
mt_pos_path = green_edges %>% select(from, to, states)

# calculate mutant degs using adjacent pairs ----------------------------------

mt_pos_path = mt_pos_path %>% 
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_), 
                                  .x = states, 
                                  model_formula_str = "~ as.factor(cell_group)", 
                                  meso_mt_pb_cds))

# saveRDS(mt_pos_path, "examples/muscle_grn/mt_pos_path_4_to_9.rds")
saveRDS(mt_pos_path, "examples/muscle_grn/mt_pos_path.rds")
mt_pos_path <- readRDS("examples/muscle_grn/mt_pos_path.rds")


# to do : dig into why did the first one fail?  from = 2, to = 8

# debug(build_pln_model_on_genes)
# pln_model_2_to_8 = build_pln_model_on_genes(genes = mt_pos_path$degs[[1]],
#                                             knockout_genes = knockouts, 
#                                             states = mt_pos_path$states[[1]],
#                                             cds = meso_mt_pb_cds)
# 
# Error in lm.wfit(covariates, log(1 + responses[, j]), weights, offset = offsets[, : 
#                                                                                   NA/NaN/Inf in 'y'
# this is potentially not a super interesting contrast anyways? 


# Now let's make a whitelist of regulators -------------------------------------

# tbx16 and msgn1
knockouts = c("ENSDARG00000007329", "ENSDARG00000070546")

allowed_reg_ids = get_allowed_regulators(meso_pb_cds, knockouts) 
mt_pos_path_gene_mod = get_gene_modules(meso_pb_cds, unlist(mt_pos_path$degs), allowed_reg_ids)

# calculate mutant DEGs along adjacent pairs 

mt_pos_path = mt_pos_path %>% 
  dplyr::mutate(gene_count_model = purrr::map2(.f = purrr::possibly(
    build_pln_model_on_genes, NA_real_),
    .x = degs,
    .y = states,
    cds = meso_mt_pb_cds,
    knockout_genes = knockouts, 
    model_formula_str = "~cell_state",
    regulatory_genes = allowed_reg_ids,
    gene_module_df = mt_pos_path_gene_mod)) 

# score these regulators 
mt_pos_path = mt_pos_path %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(rank_regulators, NA_real_), 
                                         .x = states,
                                         .y = gene_count_model, 
                                         meso_mt_pb_cds))


saveRDS(mt_pos_path, "examples/muscle_grn/transcription_kinases/mt_pos_path_regulators.rds")
mt_pos_path <- readRDS("examples/muscle_grn/transcription_kinases/mt_pos_path_regulators.rds")

# instead of adjacent pairs do DEG along the two paths ------------------------
# now instead of adjacent pairs, try to do DEG on the entire two branches -----


# library(tidygraph)
# g = green_edges %>% 
#   select(from, to) %>% 
#   igraph::graph_from_data_frame()
# 
# node_info = g %>% 
#   as_tbl_graph() %>% 
#   tidygraph::activate(nodes) %>%
#   mutate(is_leaf = node_is_leaf(), 
#          is_root = node_is_root()) %>% 
#   as.data.frame()
# 
# sp_4_to_8 = igraph::shortest_paths(g, 4, to = 8)
# sp_4_to_9 = igraph::shortest_paths(g, 4, to = 9)


two_branch_df = data.frame("from" = c(4, 4), 
                           "to" = c(8, 9)) %>% 
  mutate(path = row_number()) %>% 
  dplyr::mutate(shortest_path = purrr::map2(.f = get_shortest_path, 
                                            .x = from, .y = to, 
                                            weighted_edges)) %>% 
  dplyr::rename("root" = from, "leaf" = to) %>%
  select(root, leaf, shortest_path) %>% 
  dplyr::mutate(states = purrr::map(.f = collect_between_transition_states, 
                                    .x = shortest_path)) %>% 
  mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_),
                           .x = states,
                           model_formula_str = "~ as.factor(cell_group)",
                           meso_mt_pb_cds))

# saveRDS(two_branch_df, "examples/muscle_grn/two_branch_df.rds")
two_branch_df <- readRDS("examples/muscle_grn/two_branch_df.rds")


two_branch_all_gene_mods = get_gene_modules(meso_pb_cds, unlist(two_branch_df$degs), allowed_reg_ids)


# two_branch_df = two_branch_df %>% 
#   dplyr::mutate(gene_modules = purrr::map(.f = get_gene_modules, 
#                                          .x = degs, 
#                                          pb_cds = meso_mt_pb_cds, 
#                                          allowed_regulator_ids = allowed_reg_ids)) 
two_branch_df = two_branch_df %>% 
  dplyr::mutate(gene_count_model = purrr::map2(.f = 
    build_pln_model_on_genes,
    .x = degs,
    .y = states,
    cds = meso_mt_pb_cds,
    model_formula_str = "~cell_state",
    knockout_genes = knockouts, 
    regulatory_genes = allowed_reg_ids,
    gene_module_df = gene_modules,
    sparsity_factor = 1,
    pln_min_ratio = 1e-3)) 



# idk why this failed in the purrr
# gene_count_model = build_pln_model_on_genes(two_branch_df %>% slice(2) %>% pull(degs) %>% unlist(), 
#                          two_branch_df %>% slice(2) %>% pull(states) %>% unlist(), 
#                          cds = meso_mt_pb_cds,
#                          model_formula_str = "~cell_state",
#                          knockout_genes = knockouts, 
#                          regulatory_genes = allowed_reg_ids,
#                          gene_module_df = two_branch_all_gene_mods)
# 
# two_branch_df$gene_count_model[[2]] = gene_count_model
# saveRDS(gene_count_model, "examples/muscle_grn/fm_all_pln_model.rds")


two_branch_df = two_branch_df %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(rank_regulators, NA_real_), 
                                         .x = states,
                                         .y = gene_count_model, 
                                         meso_mt_pb_cds))




saveRDS(two_branch_df, "examples/muscle_grn/two_branch_df_regulators.rds")
two_branch_df <- readRDS("examples/muscle_grn/two_branch_df_regulators.rds")


# plot top regulators for these

for (i in 1:nrow(two_branch_df)) {
  from = two_branch_df[i,]$root
  to = two_branch_df[i,]$leaf
  file_name =  paste("meso_tbx16msgn1",from,to, sep="_")
  print(file_name)
  plot_top_regulators(two_branch_df$regulators[[i]], meso_mt_pb_cds) %>% 
    ggsave(filename = paste0("examples/muscle_grn/plots/",file_name,"two_branch.png"))
}
two_branch_df %>% colnames()


# helper functions ------------------------------------------------------------


comp_abund <- function(to, from, ccm) {
  to_abund = estimate_abundances(ccm, tibble("cell_state" = to))
  from_abund = estimate_abundances(ccm, tibble("cell_state" = from))
  return(compare_abundances(ccm, to_abund, from_abund))
}


id_to_shortname = function(gene_ids) {
  rowData(meso_mt_cds) %>% 
    as.data.frame() %>% 
    select(gene_short_name, id) %>% 
    filter(id %in% gene_ids) %>% 
    pull(gene_short_name)
}

shortname_to_id = function(short_names) {
  rowData(meso_mt_cds) %>% 
    as.data.frame() %>% 
    select(gene_short_name, id) %>% 
    filter(gene_short_name %in% short_names) %>% 
    pull(id)
}


add_plot <- function(edges, full_regulators) {
  
  edges = edges %>% dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))
  n = network(edges, directed = F, loops = F)
  n = add_node_attribute(n, full_regulators, "regulator_score")
  n %v% "gene_short_name" = id_to_shortname(network.vertex.names(n))
  
  g = ggnetwork(n, arrow.gap=0.02) 
  
  p = ggplot(mapping=aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
    geom_nodes(data = g, aes(fill = regulator_score),
               size = 3, colour = "black", shape = 21) + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high="red4") + 
    geom_nodelabel_repel(data = g, aes(label = gene_short_name)) +
    geom_edges(data = g %>% filter(edge_type != "undirected"), color="gray") +
    # draw activator edges
    geom_edges(data = g %>% filter(edge_type == "activator"),
               arrow = arrow(length = unit(2, "pt"), type="closed")) +
    # draw repressor edges 
    geom_edges(data = g %>% filter(edge_type == "repressor"),
               arrow = arrow(angle = 90, length = unit(.075, "cm"))) + 
    scale_size_identity() + 
    monocle3:::monocle_theme_opts() + theme_blank()
  
  return(p)
}

# get branch and distance info in there
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
  # find the regulators 
  
  
# inner_join(edges, adj_path, by = c("to", "from"))
  

edges %>% slice(1) %>% mutate(purrr::map(.f = find_all_edges,
                                         .x = sub_states, 
                                         path = adj_path))


# FC with respect to start


# this gives us for each ending cell state
# what is the delta_log_abundance change

delta_change_from_root = edges %>%  
  dplyr::mutate(cond_b_vs_a_tbl = 
                  purrr::map2(.f = comp_abund, 
                              .x = as.character(root),
                              .y = to,
                              ccm = adj_reg_pln_model)) %>% 
  tidyr::unnest(cond_b_vs_a_tbl) %>% 
  select(cell_group, "cell_state" = cell_state_y, 
         log_abund_y, log_abund_sd_y, 
         delta_log_abund, delta_p_value)

# get state order


  
# color nodes by this 
add_delta_change_from_root = function(n, 
                                      state, 
                                      delta_change_from_root, 
                                      p_value_threshold= 1.0) {
  
  delta_change_from_root = delta_change_from_root %>% 
    filter(cell_stat == state) %>% 
    column_to_rownames("cell_group")
  
  # if below pvalue threshold, no change
  delta_change_from_root = delta_change_from_root %>%
    mutate(delta_log_abund = ifelse(delta_p_value < p_value_threshold, delta_log_abund, 0))
  
  n %v% "delta_log_abund" = delta_change_from_root[network.vertex.names(n),"delta_log_abund"]
  
  # delta_change_from_root[network.vertex.names(n),"log_abund_y"]
  
  return(n)
  
}







# what happens if we pool the adj clusters along branches ---------------------
# let's just pool the fast muscle branch for now 

fm_branch = two_branch_df %>% 
  filter(leaf == "9") %>% 
  select(root, leaf, states)

fm_states = unlist(fm_branch$states)

fm_degs = mt_pos_path %>% 
  select(to, from, states, degs) %>% 
  filter(to %in% fm_states, from %in% fm_states) %>% 
  pull(degs) %>% 
  unlist()

mt_pos_path %>% 
  filter(to %in% fm_states, from %in% fm_states) %>% 
  select(regulators) %>% 
  tidyr::unnest(regulators) %>% 
  ggplot(aes(abs(regulator_score))) + geom_density() 


fm_regulators = mt_pos_path %>% 
  filter(to %in% fm_states, from %in% fm_states) %>% 
  select(regulators) %>% 
  tidyr::unnest(regulators) %>%
  filter(abs(regulator_score) > 0) %>% 
  pull(gene_id) %>% 
  unique() # 615 regulators , 182


# fm_pln_model = build_pln_model_on_genes(genes = fm_regulators, 
#                                         states = all_states,
#                                         cds = meso_mt_pb_cds,
#                                         knockout_genes = knockouts, 
#                                         regulatory_genes = allowed_regulator_ids,
#                                         model_formula_str = "~cell_state",
#                                         gene_module_df = mt_pos_path_gene_mod)




# saveRDS(fm_pln_model, "examples/muscle_grn/fm_pln_model_0.1.rds")
fm_pln_model <- readRDS("examples/muscle_grn/fm_pln_model.rds")


fm_path = mt_pos_path %>% 
  filter(to %in% fm_states, from %in% fm_states)

fm_path = fm_path %>% 
  dplyr::mutate(full_regulators = purrr::map(.f = purrr::possibly(rank_regulators, NA_real_), 
                                           .x = states,
                                           gene_model_ccm = fm_pln_model, 
                                           cds = meso_mt_pb_cds)) 

fm_path = fm_path %>% 
  dplyr::mutate(filter_regulators = purrr::map(full_regulators, ~filter(., regulator_score!=0)))

fm_path = fm_path %>% 
  dplyr::mutate(cond_b_vs_a_tbl = 
         purrr::map2(.f = comp_abund, 
                     .x = from,
                     .y = to,
                     ccm = fm_pln_model)) %>% 
  dplyr::mutate(filter_cond_b_vs_a_tbl = purrr::map(cond_b_vs_a_tbl, ~filter(., delta_p_value < 0.05))) %>% 
  dplyr::mutate(edges = 
         purrr::map(.f = collect_pln_gene_edges, 
                      .x = filter_cond_b_vs_a_tbl, 
                      ccm = fm_pln_model)) %>% 
  dplyr::mutate(dir_edges = purrr::map(edges, ~filter(., edge_type != "undirected")))


fm_path = fm_path %>% 
  dplyr::mutate(plot = purrr::map2(.f = add_plot, 
                                   .x = edges, 
                                   .y = full_regulators))


saveRDS(fm_path, "examples/muscle_grn/fm_path.rds")
fm_path <- readRDS("examples/muscle_grn/fm_path.rds")



  

# plot these regulators -------------------------------------------------------



# we can now pool the regulators and build a single pln model on that ---------
# currently we are just including regulators that are non-zero

reg_threshold = 0

# this should be the same for both options 
all_states = mt_pos_path %>% pull(states) %>% unlist() %>% unique()


# option 1: use adjacent DEGs

adj_regulators = mt_pos_path %>% 
                  filter(!is.na(regulators)) %>% 
                  select(regulators) %>% 
                  tidyr::unnest(regulators) %>%
                  filter(abs(regulator_score) > 0.1) %>% 
                  pull(gene_id) %>% 
                  unique()

# 215 with 0.01

adj_reg_pln_model = build_pln_model_on_genes(genes = adj_regulators, 
                                             states = all_states,
                                             cds = meso_mt_pb_cds,
                                             knockout_genes = knockouts, 
                                             regulatory_genes = allowed_regulator_ids,
                                             model_formula_str = "~cell_state",
                                             gene_module_df = mt_pos_path_gene_mod)


saveRDS(adj_reg_pln_model, "examples/muscle_grn/adj_reg_pln_model_0.1.rds")
adj_reg_pln_model <- readRDS("examples/muscle_grn/adj_reg_pln_model_0.1.rds")


adj_reg_pln_model = select_model(adj_reg_pln_model, sparsity_factor = 10)

cov_graph = adj_reg_pln_model@best_model %>% return_igraph()
cov_edges <- igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(abs(weight) > 0.01)



# go along the adjacent paths and score these + direct edges

adj_path = mt_pos_path %>% select(to, from, states) %>% 
  dplyr::mutate(full_regulators = purrr::map(.f = purrr::possibly(rank_regulators, NA_real_), 
                                             .x = states,
                                             gene_model_ccm = adj_reg_pln_model, 
                                             cds = meso_mt_pb_cds)) 

adj_path = adj_path %>% 
  dplyr::mutate(cond_b_vs_a_tbl = 
                purrr::map2(.f = comp_abund, 
                            .x = from,
                            .y = to,
                            ccm = adj_reg_pln_model)) %>% 
  dplyr::mutate(edges = 
                  purrr::map(.f = collect_pln_gene_edges, 
                             .x = cond_b_vs_a_tbl, 
                             ccm = adj_reg_pln_model)) %>%
  dplyr::mutate(filter_edges = purrr::map(edges, ~filter(., abs(pcor) > 0.01)))



rowData(meso_mt_cds) %>% 
  as.data.frame %>% filter(grepl("myod", gene_short_name))

# "ENSDARG00000010192" pax3a
# "ENSDARG00000028348" pax3b
# "ENSDARG00000030110" myod1

check_ids = c("ENSDARG00000010192", "ENSDARG00000028348", "ENSDARG00000030110")

adj_regulators %>% 
  filter(to %in% check_ids)


log_abund_threshold = 0
delta_p_value = 1.0


# gene A (from) goes up from 4 to 5

upreg_4_to_5 = adj_path %>% 
  slice(6) %>% 
  select(filter_edges) %>% 
  tidyr::unnest(filter_edges) %>% 
  filter(from_delta_log_abund > log_abund_threshold) %>% 
  filter(from_delta_p_value < delta_p_value) %>%
  select(from, to, pcor, 
         from_cell_state_x, from_cell_state_y, 
         from_log_abund_x, from_log_abund_sd_x, 
         from_log_abund_y, from_log_abund_x, 
         from_delta_log_abund, from_delta_p_value)

downreg_4_to_5 = adj_path %>% 
  slice(6) %>% 
  select(filter_edges) %>% 
  tidyr::unnest(filter_edges) %>% 
  filter(from_delta_log_abund < log_abund_threshold) %>% 
  filter(from_delta_p_value < delta_p_value) %>%
  select(from, to, pcor, 
         from_cell_state_x, from_cell_state_y, 
         from_log_abund_x, from_log_abund_sd_x, 
         from_log_abund_y, from_log_abund_x, 
         from_delta_log_abund, from_delta_p_value)

# gene B (to)
upreg_5_to_3 = adj_path %>% 
  slice(8) %>% 
  select(filter_edges) %>% 
  tidyr::unnest(filter_edges) %>% 
  filter(to_delta_log_abund > log_abund_threshold) %>% 
  filter(to_delta_p_value < delta_p_value) %>%
  select(from, to, pcor, 
         to_cell_state_x, to_cell_state_y, 
         to_log_abund_x, to_log_abund_sd_x, 
         to_log_abund_y, to_log_abund_x, 
         to_delta_log_abund, to_delta_p_value)

downreg_5_to_3 = adj_path %>% 
  slice(8) %>% 
  select(filter_edges) %>% 
  tidyr::unnest(filter_edges) %>% 
  filter(to_delta_log_abund < log_abund_threshold) %>% 
  filter(to_delta_p_value < delta_p_value) %>%
  select(from, to, pcor, 
         to_cell_state_x, to_cell_state_y, 
         to_log_abund_x, to_log_abund_sd_x, 
         to_log_abund_y, to_log_abund_x, 
         to_delta_log_abund, to_delta_p_value)


find_reg_states <- function(path, 
                            from_state, 
                            to_state, 
                            upreg = T,
                            state_1 = T, 
                            log_abund_threshold = 0,
                            delta_p_value_threshold = 1.0) {
  

  path_edges = path %>% 
    filter(from == from_state, to == to_state) %>% 
    select(filter_edges) %>% 
    tidyr::unnest(filter_edges) %>% 
    filter(to_delta_p_value < delta_p_value_threshold | 
           from_delta_p_value < delta_p_value_threshold)
  
  # either get upreg or downreg states
  if (upreg) {
    path_edges = path_edges %>% 
      filter(to_delta_log_abund > log_abund_threshold | 
             from_delta_log_abund > log_abund_threshold)
  } else {
    path_edges = path_edges %>% 
      filter(to_delta_log_abund < log_abund_threshold | 
             from_delta_log_abund < log_abund_threshold)
  }
  
  # is this state the first state (from) or the second state (to)
  if (state_1) {
    path_edges = path_edges %>%
      select_if(grepl("from", colnames(path_edges)) | colnames(path_edges) %in% c("from", "to", "pcor"))
  } else {
    path_edges = path_edges %>%
      select_if(grepl("to", colnames(path_edges)) | colnames(path_edges) %in% c("from", "to", "pcor"))
  }

  return(path_edges)
  
}

find_all_edges = function(path, sub_states, ...) {
  
  # direct all upreg states
  upreg_1_to_2 = find_reg_states(path, from = sub_states[1], to = sub_states[2], upreg = T, state_1 = T, ...)
  upreg_2_to_3 = find_reg_states(path, from = sub_states[2], to = sub_states[3], upreg = T, state_1 = F, ...)
  
  upreg_edges = inner_join(upreg_1_to_2, upreg_2_to_3, 
                           by = c("from", "to"), 
                           suffix = c(".from", ".to")) %>%
                rename("pcor" = "pcor.from") %>% 
                select(-c(pcor.to)) %>% 
                mutate(edge_type = "activator")
  
  # direct all downreg states
  downreg_1_to_2 = find_reg_states(path, from = sub_states[1], to = sub_states[2], upreg = F, state_1 = T, ...)
  downreg_2_to_3 = find_reg_states(path, from = sub_states[2], to = sub_states[3], upreg = F, state_1 = F, ...)

  downreg_edges = inner_join(downreg_1_to_2, downreg_2_to_3,
                             by = c("from", "to"),
                             suffix = c(".from", ".to")) %>%
                  rename("pcor" = "pcor.from") %>%
                  select(-c(pcor.to)) %>%
                  mutate(edge_type = "repressor")

  # stick them together
  all_edges = rbind(upreg_edges,
                    downreg_edges) %>%
    dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))

  return(all_edges)
}


# pcors should be the same -- so can drop one

upreg_edges = inner_join(upreg_4_to_5, upreg_5_to_3, 
           by = c("from", "to"), 
           suffix = c(".from", ".to")) %>%
  rename("pcor" = "pcor.from") %>% select(-c(pcor.to)) %>% 
  mutate(edge_type = "activator")

downreg_edges = inner_join(downreg_4_to_5, downreg_5_to_3, 
                         by = c("from", "to"), 
                         suffix = c(".from", ".to")) %>%
  rename("pcor" = "pcor.from") %>% select(-c(pcor.to)) %>% 
  mutate(edge_type = "repressor")

# score regulators? 
all_edges = rbind(upreg_edges, 
              downreg_edges) %>% 
        dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))

# some sort of timing of nodes
# all_edges = find_all_edges()





# n = network(all_edges, directed = F, loops = F)
#n = add_node_attribute(n, full_regulators, "regulator_score")
n %v% "gene_short_name" = id_to_shortname(network.vertex.names(n))


# some of the edges are 4 to 5
# from_cell_state_x from_cell_state_y to_cell_state_x to_cell_state_y
# 4 5 5 3
#  A   B
# color gene A node by 5? 
# color gene B node by 3? 

# or for each gene order by which state it is most highly expressed? 

# time or stage of trajectory

# which state does each gene belong to
gene_state = all_edges %>% select("gene" = from, "state" = from_cell_state_y)
      all_edges %>% select("gene" = to, "state" = to_cell_state_y)) 

gene_state = inner_join(gene_state, delta_change_from_root, by = c("gene" = "cell_group","state" = "cell_state_y"))

gene_state  = gene_state %>% group_by(gene) %>% 
  summarize(state = mean(as.numeric(state)),
            log_abund_y = mean(log_abund_y), 
            delta_log_abund = mean(delta_log_abund))

gene_state = gene_state %>% column_to_rownames("gene")

n %v% "delta_log_abund" = gene_state[network.vertex.names(n),"delta_log_abund"]




# g = fortify(n)
# fruchtermanreingold
# kamadakawai
# spring
# different layouts https://rdrr.io/cran/sna/man/gplot.layout.html
# maybe use mds with some distance ? that is distance from root? 

# practice w a mini network! 


set.seed(2022)
g = ggnetwork(n,
              layout = "mds",
              arrow.gap=0.05)
ggplot(mapping = aes(x, y, xend = xend, yend = yend, size = scaled_weight)) +
  # draw activator edges
  geom_edges(data = g %>% filter(edge_type == "activator"),
             arrow = arrow(length = unit(2, "pt"), type="closed")) +
  # draw repressor edges 
  geom_edges(data = g %>% filter(edge_type == "repressor"),
             arrow = arrow(angle = 90, length = unit(.075, "cm"))) +
  # geom_nodes(data = g,
  #           aes(color = delta_log_abund),
  #            size=3) +
  # geom_nodetext(data = g,
  #                aes(label = gene_short_name),
  #                size=3) +
  geom_nodelabel(data = g,
                 aes(
                     label = gene_short_name),
                 size=3) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high="red4") + 
  scale_size_identity() + 
  monocle3:::monocle_theme_opts() + 
  theme_blank() #+ 
  # facet_wrap(~ Frequency) +
  # theme_facet()




# for (i in nrow(adj_path)) {
#   to = adj_path %>% slice(i) %>% pull(to)
#   from = adj_path %>% slice(i) %>% pull(from)
#   file_name =  paste("meso_tbx16msgn1",from,to, sep="_")
#   p = adj_path %>% slice(i) %>% pull(plot) %>% unlist()
#   ggsave(p, filename = paste0("examples/muscle_grn/plots/",file_name,"adj_path.png"))
# }

# option 2: use 2 branch DEGs

two_branch_regulators = two_branch_df %>% 
                          filter(!is.na(regulators)) %>% 
                          select(regulators) %>% 
                          tidyr::unnest(regulators) %>%
                          filter(abs(regulator_score) > reg_threshold) %>% 
                          pull(gene_id) %>% 
                          unique()

two_branch_reg_pln_model = build_pln_model_on_genes(genes = two_branch_regulators, 
                                             states = all_states,
                                             cds = meso_mt_pb_cds,
                                             knockout_genes = knockouts, 
                                             regulatory_genes = allowed_regulator_ids,
                                             model_formula_str = "~cell_state",
                                             gene_module_df = two_branch_all_gene_mods)
saveRDS(two_branch_reg_pln_model, "examples/muscle_grn/two_branch_reg_pln_model.rds")
two_branch_reg_pln_model <- readRDS("examples/muscle_grn/two_branch_reg_pln_model.rds")


# score and direct edges between adjacent

two_branch_path =  mt_pos_path %>% select(to, from, states) %>% 
  dplyr::mutate(full_regulators = purrr::map(.f = purrr::possibly(rank_regulators, NA_real_), 
                                           .x = states,
                                           gene_model_ccm = two_branch_reg_pln_model, 
                                           cds = meso_mt_pb_cds))  
two_branch_path = two_branch_path %>% 
  dplyr::mutate(cond_b_vs_a_tbl = 
                  purrr::map2(.f = comp_abund, 
                              .x = from,
                              .y = to,
                              ccm = two_branch_reg_pln_model))

two_branch_path = two_branch_path %>%  
  dplyr::mutate(filter_cond_b_vs_a_tbl = purrr::map(cond_b_vs_a_tbl, 
                                                    ~filter(., delta_p_value < 0.05))) %>% 
  dplyr::mutate(edges = 
                  purrr::map(.f = collect_pln_gene_edges, 
                             .x = filter_cond_b_vs_a_tbl, 
                             ccm = two_branch_reg_pln_model)) %>% 
  dplyr::mutate(plot = purrr::map2(.f = add_plot, 
                                   .x = edges, 
                                   .y = full_regulators))


saveRDS(two_branch_path, "examples/muscle_grn/two_branch_path.rds")
two_branch_path <- readRDS("examples/muscle_grn/two_branch_path.rds")


for (i in nrow(two_branch_path)) {
  # p = two_branch_path %>% slice(i) %>% pull(plot) %>% unlist()
  p = two_branch_path$plot[[i]]
  to = two_branch_path %>% slice(i) %>% pull(to)
  from = two_branch_path %>% slice(i) %>% pull(from)
  file_name =  paste("meso_tbx16msgn1",from,to, sep="_")
  ggsave(p, filename = paste0("examples/muscle_grn/plots/",file_name,"two_branch_path.png"))
}



# test orienting edges --------------------------------------------------------

# 4 --> 5 -- > 3

abund_4 = estimate_abundances(fm_pln_model, tibble("cell_state" = "4"))
abund_5 = estimate_abundances(fm_pln_model, tibble("cell_state" = "5"))
abund_3 = estimate_abundances(fm_pln_model, tibble("cell_state" = "3"))

cond_4_v_5_tbl = compare_abundances(fm_pln_model, abund_4, abund_5)
cond_5_v_3_tbl = compare_abundances(fm_pln_model, abund_3, abund_5)

upreg_4_v_5 = cond_4_v_5_tbl %>% filter(delta_log_abund > 0) #, delta_p_value < 0.05)
upreg_5_v_3 = cond_5_v_3_tbl %>% filter(delta_log_abund > 0) #, delta_p_value < 0.05)

upreg_edges_4_to_5 = collect_pln_gene_edges(fm_pln_model, upreg_4_v_5) %>% filter(pcor > 0)
upreg_edges_5_to_3 = collect_pln_gene_edges(fm_pln_model, upreg_5_v_3) %>% filter(pcor > 0)

merge(upreg_edges_4_to_5, 
      upreg_edges_5_to_3,
      by = c("from", "to")) %>%
  select(from, to) %>% mutate(edge_type = "activator")



downreg_4_v_5 = cond_4_v_5_tbl %>% filter(delta_log_abund < 0) #, delta_p_value < 0.05)
downreg_5_v_3 = cond_5_v_3_tbl %>% filter(delta_log_abund < 0) #, delta_p_value < 0.05)

downreg_edges_4_to_5 = collect_pln_gene_edges(fm_pln_model, downreg_4_v_5) %>% filter(pcor > 0)
downreg_edges_5_to_3 = collect_pln_gene_edges(fm_pln_model, downreg_5_v_3) %>% filter(pcor > 0)

merge(downreg_edges_4_to_5, 
      downreg_edges_5_to_3,
      by = c("from", "to")) %>%
  select(from, to) %>% mutate(edge_type = "repressor")





# genes that go up in 4
state_1_up = edges_4_to_5 %>% filter(from_delta_log_abund > 0, from_delta_p_value < 0.05, pcor > 0)
# genes that go up in 5
state_2_up = edges_5_to_3 %>% filter(from_delta_log_abund > 0, from_delta_p_value < 0.05, pcor > 0)


left_join(state_1_up, state_2_up, by = c("from", "to")) %>% head()


# where does it go up? 
# from 4 to 5
state_1_2 = edges_4_to_5 %>% 
  filter(from_delta_log_abund > 0 & pcor > 0)


state_1_2 %>% select(from, to, pcor)

state_2_3 = edges_5_to_3 %>% 
  filter(to_delta_log_abund > 0 & pcor > 0)

# gene A goes up in 4 to 5
# gene B goes up in 5 to 3



# to_delta_log_abund goes up 


  
edges_5_to_3 %>% head(1) 


# from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor < 0 ~ "directed_to_from",
# from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor < 0 ~ "directed_from_to"

# where does it go down 
edges_4_to_5 %>% filter()


# -----------------------------------------------------------------------------







# given a node and a graph, travel to root



reg_threshold = 0 

two_branch_df %>% 
  filter(!is.na(regulators)) %>% 
  select(regulators) %>% 
  tidyr::unnest(regulators) %>% 
  ggplot(aes(abs(regulator_score))) + geom_density() + geom_vline(xintercept = 0.25)

all_regulators = two_branch_df %>% 
  filter(!is.na(regulators)) %>% 
  select(regulators) %>% 
  tidyr::unnest(regulators) %>%
  filter(abs(regulator_score) > 0.25) %>% 
  pull(gene_id) %>% 
  unique()

# 1348, 280 filtered at 0.25

all_states = two_branch_df %>% pull(states) %>% unlist() %>% unique()
all_gene_modules = get_gene_modules(meso_pb_cds, 
                                     unlist(two_branch_df$degs), allowed_regulator_ids)


all_reg_pln_model = build_pln_model_on_genes(genes = all_regulators, 
                                             states = all_states,
                                             cds = meso_mt_pb_cds,
                                             model_formula_str = "~cell_state",
                                             knockout_genes = c("ENSDARG00000007329", 
                                                                "ENSDARG00000070546"), 
                                             regulatory_genes = allowed_regulator_ids,
                                             gene_module_df = all_gene_modules,
                                             sparsity_factor = 1,
                                             pln_min_ratio = 1e-3)

saveRDS(all_reg_pln_model, "examples/muscle_grn/all_reg_pln_model_2_branch.rds")

# now rank the regulators in this model 
# direct edges 


all_reg_pln_model


two_branch_reg = two_branch_df %>% 
  mutate(path = row_number()) %>%
  select(path, "root" = neg_from, "leaf" = neg_to, 
        degs, regulators, shortest_path) %>% 
  tidyr::unnest(shortest_path) %>% 
  dplyr::mutate(adj_states = purrr::map2(.f = purrr::possibly(
                collect_transition_states, NA_real_), 
                .x = from,
                .y = to )) %>% 
  dplyr::mutate(full_regulators = purrr::map(.f = purrr::possibly(rank_regulators, NA_real_), 
                                             .x = adj_states,
                                             gene_model_ccm = all_reg_pln_model, 
                                             cds = meso_mt_pb_cds))


two_branch_reg = two_branch_reg %>% 
  dplyr::mutate(filter_regulators = purrr::map(full_regulators, ~filter(., regulator_score!=0)))




two_branch_reg = two_branch_reg %>% 
  mutate(cond_b_vs_a_tbl = 
      purrr::map2(.f =  comp_abund, 
            .x = from,
            .y = to,
            ccm = all_reg_pln_model)) 

two_branch_reg = two_branch_reg %>% 
  mutate(dir_edges = 
      purrr::map(.f = collect_pln_gene_edges, 
             .x = cond_b_vs_a_tbl, 
             ccm = all_reg_pln_model))

gene_edges = two_branch_reg %>% 
  slice(i) %>% 
  select(dir_edges) %>% 
  tidyr::unnest(dir_edges) 

reg_score_df = two_branch_reg %>% 
  slice(i) %>% 
  tidyr::unnest(full_regulators) %>% 
  select(gene_id, regulator_score)

n = network(gene_edges, directed = F, loops = F)
n = add_node_attribute(n, reg_score_df, "regulator_score")
n %v% "gene_short_name" = id_to_shortname(network.vertex.names(n))

g = ggnetwork(n, arrow.gap=0.02) 



ggplot(mapping=aes(x, y, xend = xend, yend = yend)) +
    geom_nodes(data = g, aes(fill = regulator_score),
               size = 3, colour = "black", shape = 21) + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high="red4") + 
    geom_nodelabel_repel(data = g, aes(label = gene_short_name)) + 
    geom_edges(data = g %>% filter(edge_type != "undirected"), color="gray") +
    # draw activator edges
    geom_edges(data = g %>% filter(edge_type == "activator"),
               arrow = arrow(length = unit(2, "pt"), type="closed")) +
    # draw repressor edges 
    geom_edges(data = g %>% filter(edge_type == "repressor"),
               arrow = arrow(angle = 90, length = unit(.075, "cm"))) + 
    scale_size_identity() + 
    monocle3:::monocle_theme_opts() + theme_blank()


# rank regulators along the path 

# saveRDS(two_branch_reg, "examples/muscle_grn/two_branch_reg.rds")
# two_branch_reg <- readRDS("examples/muscle_grn/two_branch_reg.rds")



























# include all states along the path instead ------------------------------------

edges_along_path = neg_rec_edges %>% 
  dplyr::mutate(shortest_path = purrr::map2(.f = get_shortest_path, 
                                            .x = from, .y = to, 
                                            weighted_edges)) %>%
  dplyr::rename("neg_from" = from, "neg_to" = to) %>%
  select(neg_from, neg_to, shortest_path) %>% 
  dplyr::mutate(states = purrr::map(.f = collect_between_transition_states, 
                                    .x = shortest_path))

edge_mt_degs =  edges_along_path %>%
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_),
                                 .x = states,
                                 model_formula_str = "~ as.factor(cell_group)",
                                 meso_mt_pb_cds))


edge_mt_degs <- readRDS("examples/muscle_grn/mt_edges_along_path.rds")
# saveRDS(edge_mt_degs, "examples/muscle_grn/mt_edges_along_path.rds")
# saveRDS(edge_mt_degs, "examples/muscle_grn/mt_edges_along_path_gene_target.rds")

edge_mt_degs = edge_mt_degs %>%
  dplyr::mutate(gene_count_model = purrr::map2(.f = purrr::possibly(build_pln_model_on_genes, NA_real_),
                .x = degs,
                .y = states,
                meso_mt_pb_cds,
                "~cell_state",
                regulatory_genes = allowed_regulator_ids,
                gene_module_df = all_gene_modules,
                sparsity_factor=1,
                pln_min_ratio=1e-3)) 

edge_mt_degs = edge_mt_degs %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(rank_regulators, NA_real_), 
                                         .x = states,
                                         .y = gene_count_model, 
                                         meso_mt_pb_cds))

saveRDS(edge_mt_degs, "examples/muscle_grn/mt_edges_along_path_regulators.rds")

for (i in 1:nrow(edge_mt_degs)) {
  from = edge_mt_degs[i,]$neg_from
  to = edge_mt_degs[i,]$neg_to
  file_name =  paste("meso_tbx16msgn1_all_",from,to, sep="_")
  plot_top_regulators(edge_mt_degs$regulators[[i]], meso_mt_pb_cds) %>% 
    ggsave(filename = paste0("examples/muscle_grn/plots/",file_name,".png"))
}


# build model on all 
# regulators_2 = edge_mt_degs %>% pull(degs) %>% unlist() %>% unique()
# states_2 = edge_mt_degs %>% pull(states) %>% unlist() %>% unique()
# 
# all_reg_pln_model_2 = build_pln_model_on_genes(genes = regulators_2,
#                                              states = states_2,
#                                              cds = meso_mt_pb_cds)
#                                           
# 
# saveRDS(all_reg_pln_model_2, "all_reg_pln_model_2.rds")




# # Debugging results -----------------------------------------------------------
# 
# # 1. do the FC estimates agree
# 
# find_coefs_between_states = function(states,
#                                     cds,
#                                     model_formula_str = "~cell_group",
#                                     gene_whitelist = NULL,
#                                     q_value_thresh = 0.05,
#                                     effect_thresh = 2,
#                                     cores = 1) {
#   cds_states_tested = cds[,as.character(colData(cds)$cell_group) %in% states]
#   models = monocle3::fit_models(cds_states_tested,
#                                 model_formula_str=model_formula_str,
#                                 weights=colData(cds_states_tested)$num_cells_in_group,
#                                 cores=cores)
#   coefs = monocle3::coefficient_table(models) %>%
#     filter(
#       !grepl("Intercept", term) &
#              q_value < q_value_thresh &
#              abs(normalized_effect) > effect_thresh) %>%
#     select(gene_id, gene_short_name, term, q_value, normalized_effect, estimate)
#   
#   return(coefs)
# }
# 
# # let's just to from 3 to 1 for speed
# 
# green_edge_debug = green_edges %>% 
#   filter(from == 3, to == 1 ) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_pb_cds, 
#                                    q_value_thresh = 1.0)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# # saveRDS(green_edge_debug, "examples/muscle_grn/green_edge_debug.rds")
# # green_edge_debug <- readRDS("examples/muscle_grn/green_edge_debug.rds")
# 
# mt_pos_path_debug = mt_pos_path %>% 
#   filter(from == 3, to == 1 ) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_mt_pb_cds, 
#                                    q_value_thresh = 0.05, 
#                                    effect_thresh = 0)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# # saveRDS(mt_pos_path_debug, "examples/muscle_grn/mt_pos_path_debug.rds")
# # mt_pos_path_debug <- readRDS("examples/muscle_grn/mt_pos_path_debug.rds")
# 
# 
# wt_df = green_edge_debug %>% 
#   select(from, to, gene_short_name, estimate, q_value) %>% 
#   tidyr::unnest(c(gene_short_name, estimate, q_value))
# 
# mt_df = mt_pos_path_debug %>% 
#   select(from, to, gene_short_name, estimate, q_value) %>% 
#   tidyr::unnest(c(gene_short_name, estimate, q_value))
# 
# wt_mt_df = full_join(wt_df, mt_df, 
#                      by = c("from", "to", "gene_short_name"), 
#                      suffix = c(".wt", ".mt"))
# 
# wt_mt_df %>%
#   filter(!is.na(estimate.mt) , !is.na(estimate.wt)) %>% 
#   filter(q_value.wt < 0.05) %>%
#   ggplot(aes(estimate.wt, estimate.mt, color=q_value.mt)) + 
#     geom_point() + monocle3:::monocle_theme_opts() 
# 
# 
# # 2. how many pseudobulks are there -------------------------------------------
# 
# g1 = colData(meso_mt_pb_cds) %>% 
#   as.data.frame() %>% 
#   group_by(group_id) %>% 
#   filter(sum(num_cells_in_group) >= 25) %>% 
#   ungroup() %>% 
#   pull(group_id)
# 
# 
# colData(meso_pb_cds) %>% 
#   as.data.frame() %>% 
#   group_by(group_id) %>% 
#   summarise(total_cells = sum(num_cells_in_group))
# 
# # these are filtered now 
# meso_pb_cds = pseudobulk(meso_ccs)
# meso_mt_pb_cds = pseudobulk(meso_mt_ccs)
# 
# green_edges_debug_2 = green_edges %>% 
#   filter(from == 3, to == 1 ) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_pb_cds, 
#                                    model_formula_str = "~ as.factor(cell_group)", 
#                                    q_value_thresh = 0.05)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# 
# colData(meso_mt_cds) %>% 
#   as.data.frame() %>% 
#   group_by(embryo,cluster, gene_target)  
#   
# 
# # add gene_target to meso_mt_pb_cds
# meso_mt_pb_cds <- add_covariate()
# 
# mt_pos_path_debug_1 = mt_pos_path %>% 
#   filter(from == 3, to == 1) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_mt_pb_cds, 
#                                    model_formula_str = "~ as.factor(cell_group)", 
#                                    q_value_thresh = 0.05, 
#                                    effect_thresh = 2)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# mt_pos_path_debug_2 = mt_pos_path %>% 
#   filter(from == 3, to == 1) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_mt_pb_cds, 
#                                    model_formula_str = "~ as.factor(cell_group) + gene_target", 
#                                    q_value_thresh = 0.05, 
#                                    effect_thresh = 2)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# 
# a = green_edges_debug_2 %>% select(coefs) %>% 
#   tidyr::unnest(coefs) %>% 
#   select(term, gene_short_name, estimate, normalized_effect, q_value)
# 
# b = mt_pos_path_debug_2 %>% select(coefs) %>% 
#   tidyr::unnest(coefs) %>% 
#   select(term, gene_short_name, estimate, normalized_effect, q_value)
# 
# b %>% 
#   filter(grepl("cell_group", term) & 
#            q_value < 0.05 &
#            abs(normalized_effect) > 2) %>% nrow()
# 
# b %>%
#   filter(grepl("gene_target", term) &
#            q_value < 0.05 &
#   abs(normalized_effect) > 2) %>% nrow()
# 
# 
# full_join(a, b, by = c("gene_short_name", "term"), suffix = c(".wt", ".mt")) %>% 
#   filter(grepl("Intercept", term)) %>% 
#   filter(q_value.wt < 0.05, q_value.mt < 0.05) %>%  
#   ggplot() + 
#   geom_point(aes(normalized_effect.wt, normalized_effect.mt)) 
# 
# 
# full_join(a, b, by = "gene_short_name", suffix = c(".wt", ".mt")) %>% 
#   filter(!is.na(estimate.wt)) %>% 
#   ggplot() + 
#   geom_point(aes(normalized_effect.wt, normalized_effect.mt)) + 
#   geom_hline(yintercept = 2) + 
#   geom_hline(yintercept = -2) + 
#   geom_vline(xintercept = 2) + 
#   geom_vline(xintercept = -2) + 
#   ylim(-3,3) + xlim(-3,3)
# 
# 
# # 3. what happens when you run WT + ctrl-inj/MT together
# 
# colnames(meso_mt_cds) = gsub("_1|_2", "", colnames(meso_mt_cds))
# meso_comb_cds <- combine_cds(list(meso_cds[,colData(meso_cds)$gene_target == "ctrl-uninj"], 
#                                   meso_mt_cds), cell_names_unique = T, keep_reduced_dims = T)
# colData(meso_comb_cds)$sample = NULL
# 
# meso_comb_ccs = new_cell_count_set(meso_comb_cds,
#                                    sample_group = "embryo",
#                                    cell_group = "cluster")
# 
# meso_comb_pb_cds <- pseudobulk(meso_comb_ccs, gene_ids = c("tbx16", "msgn1"))
# 
# 
# meso_comb_pos_path = green_edges %>% 
#   select(from, to, states) %>% 
#   filter(from == 3, to == 1 ) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    cds = meso_comb_pb_cds, 
#                                    q_value_thresh = 0.05)) %>% 
#   dplyr::mutate(gene_short_name = purrr::map(coefs, ~pull(., gene_short_name))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate))) %>% 
#   dplyr::mutate(q_value = purrr::map(coefs, ~pull(., q_value)))
# 
# comb_df = meso_comb_pos_path %>% 
#   select(from, to, gene_short_name, estimate, q_value) %>% 
#   tidyr::unnest(c(gene_short_name, estimate, q_value))
# 
# wt_mt_comb_df = full_join(wt_mt_df, comb_df, 
#                      by = c("from", "to", "gene_short_name"))
# 
# wt_mt_comb_df %>% 
#   filter(q_value.wt < 0.05) %>% 
#   ggplot(aes(estimate.wt, estimate, color=q_value.mt)) + 
#   geom_point() + monocle3:::monocle_theme_opts() 
# 
# 
# wt_mt_comb_df %>% filter(q_value.wt < 0.05) %>% pull(gene_short_name) %>% length()
# comb_genes = wt_mt_comb_df %>% filter(q_value < 0.05) %>% pull(gene_short_name) %>% length()
# 
# c = meso_comb_pos_path %>% select(coefs) %>% 
#   tidyr::unnest(coefs) %>% 
#   select(gene_short_name, estimate, normalized_effect, q_value)
# 
#   
# # 4. for each neg pair, test over complete pos path --------------------------
# 
# meso_mt_pb_cds = add_covariate(meso_mt_ccs, meso_mt_pb_cds, "gene_target")
# edges_along_path = neg_rec_edges %>% 
#   filter(from == 4, to == 9) %>% 
#   dplyr::mutate(shortest_path = purrr::map2(.f = get_shortest_path, 
#                                             .x = from, .y = to, 
#                                             weighted_edges)) %>%
#   dplyr::rename("neg_from" = from, "neg_to" = to) %>%
#   select(neg_from, neg_to, shortest_path) %>% 
#   dplyr::mutate(states = purrr::map(.f = collect_between_transition_states, 
#                                     .x = shortest_path))
# 
# edge_mt_degs =  edges_along_path %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    q_value_thresh = 0.05, 
#                                    model_formula_str = "~ as.factor(cell_group)", 
#                                    effect_thresh = 2, 
#                                    meso_mt_pb_cds)) %>% 
#   dplyr::mutate(gene_id = purrr::map(coefs, ~pull(., gene_id))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate)))
# 
# edge_mt_degs_2 =  edges_along_path %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    q_value_thresh = 0.05, 
#                                    model_formula_str = "~ as.factor(cell_group) + gene_target", 
#                                    effect_thresh = 2, 
#                                    meso_mt_pb_cds)) %>% 
#   dplyr::mutate(gene_id = purrr::map(coefs, ~pull(., gene_id))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate)))
# 
# edge_mt_degs_2 %>% 
#   select(coefs) %>% 
#   tidyr::unnest(coefs) %>% filter(grepl("cell_group", term)) %>% 
#   filter(abs(normalized_effect) > 2)
# 
#   
# # saveRDS(edge_mt_degs, "examples/muscle_grn/mt_pos_path_between.rds")
# # edge_mt_degs <- readRDS("examples/muscle_grn/mt_pos_path_between.rds")
# 
# edge_wt_degs = edges_along_path %>%
#   dplyr::rename("neg_from" = from, "neg_to" = to) %>%
#   select(neg_from, neg_to, shortest_path) %>% 
#   dplyr::mutate(states = purrr::map(.f = collect_between_transition_states, 
#                                     .x = shortest_path)) %>% 
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    meso_pb_cds)) %>% 
#   dplyr::mutate(gene_id = purrr::map(coefs, ~pull(., gene_id))) %>% 
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate)))
# 
# saveRDS(edge_wt_degs, "examples/muscle_grn/wt_pos_path_between.rds")
# # edge_wt_degs <- readRDS("examples/muscle_grn/wt_pos_path_between.rds")
# 
# edge_combo_degs = neg_rec_edges %>%
#   filter(from == 4, to == 9) %>%
#   dplyr::mutate(shortest_path = purrr::map2(.f = get_shortest_path,
#                                             .x = from, .y = to,
#                                             weighted_edges)) %>%
#   dplyr::rename("neg_from" = from, "neg_to" = to) %>%
#   select(neg_from, neg_to, shortest_path) %>%
#   dplyr::mutate(states = purrr::map(.f = collect_between_transition_states,
#                                     .x = shortest_path)) %>%
#   dplyr::mutate(coefs = purrr::map(.f = purrr::possibly(find_coefs_between_states, NA_real_), 
#                                    .x = states, 
#                                    meso_comb_pb_cds)) %>%
#   dplyr::mutate(gene_id = purrr::map(coefs, ~pull(., gene_id))) %>%
#   dplyr::mutate(estimate = purrr::map(coefs, ~pull(., estimate)))
# 
# # saveRDS(edge_combo_degs, "examples/muscle_grn/comb_pos_path_between.rds")
# 
# 
# plot_cells(meso_comb_cds[,colData(meso_comb_cds)$gene_target == c("ctrl-uninj", "ctrl-inj")],
#            color_cells_by = "cluster") + facet_wrap(~gene_target)
# 
# 
# # ------------------------------------------------------------------------------
# 
# neg_rec_edges = hooke:::collect_pln_graph_edges(meccm, cond_b_vs_a_tbl) %>% 
#   as_tibble %>%
#   filter(edge_type != "undirected" &
#            to_delta_p_value < p_value_threshold &
#            from_delta_p_value < p_value_threshold) 



