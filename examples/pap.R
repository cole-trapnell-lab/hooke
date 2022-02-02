library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

#pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/pap.al_2021-09-15.cds.RDS")


pap_cds = detect_genes(pap_cds)

# assign best celltype column and reduce dims
colData(pap_cds)$cell_type = colData(pap_cds)$CW_assignedCellType
colData(pap_cds)$cluster = monocle3::clusters(pap_cds)

colData(pap_cds)$Size_Factor = size_factors(pap_cds)


colData(pap_cds)$experiment = colData(pap_cds)$sample
colData(pap_cds)$sample = NULL

plot_cells(pap_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("pap_cell_types.png", width=3, height=3)

plot_cells(pap_cds, color_cells_by="Genotype", show_trajectory_graph=FALSE) + facet_wrap(~Genotype)


plot_cells(pap_cds, genes=c("Csf2",
                            "Csf2ra",
                            "Csf2rb",
                            "Chil3",
                            "Lpl",
                            "Car4",
                            "Apoe",
                            "Fabp4",
                            "Fabp5",
                            "F13a1"),
           show_trajectory_graph=FALSE)

#plot_cells(cds)

ccs = new_cell_count_set(pap_cds,
                         sample_group = "sampleName",
                         cell_group = "cell_type")


ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "Genotype",
                            nuisance_model_formula_str = "batch")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=10)
plot(model(ccm, "reduced"), output="corrplot")

cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_csf2rb = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))


cond_ra_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)

plot_contrast(ccm, cond_ra_vs_wt_tbl, scale_shifts_by="sender", q_value_thresh=0.01) +
  theme_minimal() + monocle3:::monocle_theme_opts() +
  ggsave("pap_mac_shift.png", width=7, height=6)


full_edges = hooke:::collect_pln_graph_edges(ccm, cond_ra_vs_wt_tbl, model_for_pcors = "full") %>% dplyr::select(from, to, pcor)
reduced_edges = hooke:::collect_pln_graph_edges(ccm, cond_ra_vs_wt_tbl, model_for_pcors = "reduced") %>% dplyr::select(from, to, pcor)

full_vs_reduced = full_join(full_edges, reduced_edges, by=c("from", "to")) %>% tidyr::replace_na(list(pcor.x = 0, pcor.y = 0))
qplot(pcor.x, pcor.y, data=full_vs_reduced) + geom_abline() + geom_smooth(method="lm") + xlab("Partial correlation (full model)") + ylab("Partial correlation (reduced model)")



cond_rb_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2rb)

plot_contrast(ccm, cond_rb_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.05)

pred_abund_mat = cbind(cond_wt$log_abund, cond_csf2ra$log_abund, cond_csf2rb$log_abund)
colnames(pred_abund_mat) = c("WT", "Csfr2ra-/-", "Csfr2rb-/-")

pheatmap::pheatmap(pred_abund_mat, scale="row")

######

# Test DE analysis:

agg_expr_mat = monocle3::aggregate_gene_expression(ccs@cds,
                                                   cell_group_df = tibble::rownames_to_column(ccs@metadata[["cell_group_assignments"]]),
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0,
                                                   cell_agg_fun="mean")
agg_coldata = ccs@metadata[["cell_group_assignments"]] %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
row.names(agg_coldata) = colnames(agg_expr_mat)

gene_set = msigdbr(species = "Mus musculus", subcategory = "GO:MF")
transcription_regulators = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("Transcription", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique %>% sort

pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
pseudobulk_cds = preprocess_cds(pseudobulk_cds)
pseudobulk_cds = reduce_dimension(pseudobulk_cds)
rowData(pseudobulk_cds)$TF = rowData(pseudobulk_cds)$gene_short_name %in% transcription_regulators

plot_cells(pseudobulk_cds, color_cells_by="cell_group")

pseudo_fit = fit_models(pseudobulk_cds,
                            "~cell_group",
                            weights=colData(pseudobulk_cds)$num_cells_in_group,
                            cores=4
)
pseudo_coefs = coefficient_table(pseudo_fit)

expressed_genes = rowData(pap_cds)[rowData(pap_cds)$num_cells_expressed > 100,] %>% as.data.frame %>% pull(id)

all_gene_modules = find_gene_modules(pseudobulk_cds[expressed_genes,], resolution=1e-2)
all_gene_modules = left_join(all_gene_modules, rowData(pseudobulk_cds) %>% as.data.frame, by="id")
qplot(dim_1, dim_2, color=module, data=all_gene_modules)


heathy_vs_pap_macs = pseudobulk_cds[,colData(pseudobulk_cds)$cell_group %in% c("Healthy alveolar macrophages","hPAP alveolar macrophages")]

mac_pseudo_fit = fit_models(heathy_vs_pap_macs,
                            "~cell_group",
                            weights=colData(heathy_vs_pap_macs)$num_cells_in_group,
                            cores=4
)

mac_pseudo_coefs = coefficient_table(mac_pseudo_fit)
#mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(grepl("cell_group", term) & q_value < 0.05)
#mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(gene_short_name %in% c("Pparg", "Abca1", "Abcg1", "Csf2ra"))
#mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(term != "(Intercept)" & gene_short_name %in% c("Lpl", "Siglecf", "Apoe", "Car4", "F13a1", "Fabp4", "Fabp5"))

mac_gene_modules = monocle3::find_gene_modules(heathy_vs_pap_macs)
qplot(dim_1, dim_2, color=module, data=mac_gene_modules)


