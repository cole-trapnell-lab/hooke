library(monocle3)
library(msigdbr)
# library(hooke)

devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")

setwd("~/OneDrive/UW/Trapnell/hooke/")

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
  

# do a single path first --------

pos_path = get_weighted_edges(meso_ccm, pos_edges) %>% 
  get_shortest_path(from="4", to="9") %>% as_tibble()

pos_path = pos_path %>% 
  dplyr::mutate(states = purrr::map2(.f = purrr::possibly(
  collect_transition_states, NA_real_), .x = from,
  .y = to ))

# do degs along this path

pos_path = pos_path %>% 
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(
  find_degs_between_states, NA_real_), .x = states, meso_pseudobulk_cds))

# saveRDS(pos_path, "examples/muscle_grn/wt_pos_path.rds")
pos_path <- readRDS("examples/muscle_grn/wt_pos_path.rds")

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


genes = c("tbx16", "msgn1")


meso_mt_pb_cds = pseudobulk(meso_mt_ccs, gene_ids = c("tbx16", "msgn1"))

# currently use the single pos path
mt_pos_path = pos_path %>% select(from, to, states)

mt_pos_path = mt_pos_path %>% 
  dplyr::mutate(degs = purrr::map(.f = purrr::possibly(
    find_degs_between_states, NA_real_), .x = states, meso_mt_pb_cds))

# saveRDS(mt_pos_path, "examples/muscle_grn/mt_pos_path.rds")
mt_pos_path <- readRDS("examples/muscle_grn/mt_pos_path.rds")

# ----------

# Now let's make a whitelist of regulators:

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
qplot(dim_1, dim_2, color=module, data=all_gene_modules)


states = mt_pos_path %>% filter(from == 3) %>% tidyr::unnest(states) %>% pull(states)
genes = mt_pos_path %>% filter(from == 3) %>% tidyr::unnest(degs) %>% pull(degs)


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
  dplyr::mutate(regulators = purrr::map2(.f = purrr::possibly(
    rank_regulators, NA_real_), .x = states,
    .y = gene_count_model, pseudobulk_cds))

plot_top_regulators(mt_pos_path$regulators[[2]], meso_mt_pb_cds)

# not enough for model? 
# plot_top_regulators(mt_pos_path$regulators[[1]], meso_mt_pb_cds)
# plot_top_regulators(mt_pos_path$regulators[[3]], meso_mt_pb_cds)
# plot_top_regulators(mt_pos_path$regulators[[4]], meso_mt_pb_cds)

