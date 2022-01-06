library(monocle3)
library(msigdbr)
# library(hooke)

devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")

setwd("~/OneDrive/UW/Trapnell/hooke/")

if (!file.exists("examples/R_objects/meso_18_24_cds.rd")) {
  drive_download("https://drive.google.com/file/d/1QMeWa2Pzg9coQNaybTMELwFH6MzR58uo/")
}

if (!file.exists("examples/R_objects/meso_tbx16-like_18_24_proj_cds.rds")) {
  drive_download("https://drive.google.com/file/d/1P-W-HXD6f7BhAqD5W-53DcHYbh4pLeuS/")
}

# ------------------------------------------------------------------------------

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

# perform a contrast between 18 and 24 hours

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

mt_pos_path = mt_pos_path %>% 
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(find_degs_between_states, NA_real_), 
                                  .x = states, 
                                  model_formula_str = "~ as.factor(cell_group)", 
                                  meso_mt_pb_cds))

# saveRDS(mt_pos_path, "examples/muscle_grn/mt_pos_path_4_to_9.rds")
saveRDS(mt_pos_path, "examples/muscle_grn/mt_pos_path.rds")
mt_pos_path <- readRDS("examples/muscle_grn/mt_pos_path.rds")


plot_cells(meso_mt_cds, color_cells_by = "cluster_transfer") + theme(legend.position = "bottom")

# Now let's make a whitelist of regulators -------------------------------------

gene_set = msigdbr(species = "Danio rerio", subcategory = "GO:MF")

transcription_regulators = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("_TRANSCRIPTION", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

kinases = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("_KINASE", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

allowed_regulator_symbols = c(transcription_regulators,
                              kinases)

allowed_regulator_ids = rowData(meso_mt_pb_cds) %>% as.data.frame %>% filter(gene_short_name %in% allowed_regulator_symbols) %>% pull(id)
allowed_regulator_ids = c(allowed_regulator_ids, c("ENSDARG00000007329","ENSDARG00000070546"))


# Project genes into 2D UMAP space. We'll use this to set up penalties between genes
all_degs = c(unlist(mt_pos_path$degs), allowed_regulator_ids) %>% unique
all_gene_modules = find_gene_modules(meso_mt_pb_cds[all_degs,], resolution=1e-2)
all_gene_modules = left_join(all_gene_modules, rowData(meso_mt_pb_cds) %>% as.data.frame, by="id")
qplot(dim_1, dim_2, color=module, data=all_gene_modules) + monocle3:::monocle_theme_opts()


mt_pos_path = mt_pos_path %>%
  dplyr::mutate(gene_count_model = purrr::map2(.f = purrr::possibly(
    build_pln_model_on_genes, NA_real_),
    .x = degs,
    .y = states,
    meso_mt_pb_cds,
    "~cell_state",
    regulatory_genes = allowed_regulator_ids,
    gene_module_df = all_gene_modules,
    sparsity_factor=1,
    pln_min_ratio=1e-3)) 

mt_pos_path = mt_pos_path %>%
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(rank_regulators, NA_real_), 
                                         .x = states,
                                         .y = gene_count_model, 
                                         meso_mt_pb_cds))


for (i in 1:nrow(mt_pos_path)) {
  from = mt_pos_path[i,]$from
  to = mt_pos_path[i,]$to
  file_name =  paste("meso_tbx16msgn1",from,to, sep="_")
  plot_top_regulators(mt_pos_path$regulators[[i]], meso_mt_pb_cds) %>% 
    ggsave(filename = paste0("examples/muscle_grn/plots/",file_name,".png"))
}

# saveRDS(edge_mt_degs, "examples/muscle_grn/mt_edges_along_path_gene_target.rds")

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



