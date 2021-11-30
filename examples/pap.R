library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/pap.al_2021-09-15.cds.RDS")


# assign best celltype column and reduce dims
colData(pap_cds)$cell_type = colData(pap_cds)$CW_assignedCellType
colData(pap_cds)$cluster = monocle3::clusters(pap_cds)

colData(pap_cds)$Size_Factor = size_factors(pap_cds)


colData(pap_cds)$experiment = colData(pap_cds)$sample
colData(pap_cds)$sample = NULL

plot_cells(pap_cds, color_cells_by="cell_type")
#plot_cells(cds)

ccs = new_cell_count_set(pap_cds,
                         sample_group = "sampleName",
                         cell_group = "cell_type")


ccm  = new_cell_count_model(ccs,
                            model_formula_str = "~Genotype + Age + batch + experiment")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=1)


cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_csf2rb = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))


cond_ra_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)

plot_contrast(ccm, cond_ra_vs_wt_tbl, scale_shifts_by="none", p_value_thresh=0.1)


cond_rb_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)

plot_contrast(ccm, cond_rb_vs_wt_tbl, scale_shifts_by="none")

pred_abund_mat = cbind(cond_wt$log_abund, cond_csf2ra$log_abund, cond_csf2rb$log_abund)
colnames(pred_abund_mat) = c("WT", "Csfr2ra-/-", "Csfr2rb-/-")

pheatmap::pheatmap(pred_abund_mat, scale="row")


# Test simulations:

test_cond_a = tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1, Offset=1)
test_cond_b = tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1, Offset=1)

pred_mu_cond_a = as.vector(predict(model(ccm), newdata=test_cond_a))
pred_mu_cond_b = as.vector(predict(model(ccm), newdata=test_cond_b))

pred_sigma  = coef(model(ccm), type="covariance")


sim_reps_a = rPLN(n=100, mu=pred_mu_cond_a, Sigma=pred_sigma, depths=mean(colSums(counts(ccs))))
colnames(sim_reps_a) = colnames(pred_sigma)

sim_reps_b = rPLN(n=100, mu=pred_mu_cond_b, Sigma=pred_sigma, depths=mean(colSums(counts(ccs))))
colnames(sim_reps_b) = colnames(pred_sigma)


ggplot2::qplot(colMeans(sim_reps_a), colMeans(sim_reps_b), label=colnames(sim_reps_b), log="xy") +
  ggplot2::geom_text() +
  ggplot2::geom_abline() +
  ggplot2::xlab("Csf2ra-/-") +
  ggplot2::ylab("WT")


ggplot2::qplot(exp(cond_csf2ra$log_abund), colMeans(sim_reps_b), label=colnames(sim_reps_b), log="xy") +
  ggplot2::geom_text() +
  ggplot2::geom_abline() +
  ggplot2::xlab("Csf2ra-/- empirical") +
  ggplot2::ylab("Csf2ra-/- simulated")

qplot(condition, cells, geom="boxplot", data=rbind(data.frame(condition="WT", cells=sim_reps_a[,"T cells"]),
           data.frame(condition="Csf2ra-/-", cells=sim_reps_b[,"T cells"])))


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

pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccs@cds) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
pseudobulk_cds = preprocess_cds(pseudobulk_cds)
pseudobulk_cds = reduce_dimension(pseudobulk_cds)

plot_cells(pseudobulk_cds, color_cells_by="cell_group")
plot_genes_violin(pseudobulk_cds[rowData(pseudobulk_cds)$gene_short_name %in% c("Cd68", "Pparg", "Sftpc"),], group_cells_by="cell_group") + coord_flip()

all_gene_modules = find_gene_modules(pseudobulk_cds)
all_gene_modules = left_join(all_gene_modules, rowData(pseudobulk_cds) %>% as.data.frame, by="id")
qplot(dim_1, dim_2, color=module, data=all_gene_modules)


heathy_vs_pap_macs = pseudobulk_cds[,colData(pseudobulk_cds)$cell_group %in% c("Healthy alveolar macrophages","hPAP alveolar macrophages")]
heathy_vs_pap_macs = detect_genes(heathy_vs_pap_macs)
heathy_vs_pap_macs = heathy_vs_pap_macs[rowData(heathy_vs_pap_macs)$num_cells_expressed > 1,]
#monocle3:::size_factors(heathy_vs_pap_macs) = colData(heathy_vs_pap_macs)$num_cells_in_group

# Look at mean-variance properties in the pseudobulks:

norm_agg_mat = normalized_counts(pseudobulk_cds, norm_method="size_only", pseudocount = 0)
pb_means = rowMeans(norm_agg_mat)
pb_vars = rowVars(norm_agg_mat)
pb_sds = rowSds(norm_agg_mat)
pb_cvs = pb_sds / pb_means
qplot(pb_means, pb_cvs, log="xy")
qplot(pb_means, pb_vars, log="xy") +geom_abline(color="blue")

mac_pseudo_fit = fit_models(heathy_vs_pap_macs,
                            "~cell_group",
                            weights=colData(heathy_vs_pap_macs)$num_cells_in_group,
                            cores=4
                            )

mac_pseudo_coefs = coefficient_table(mac_pseudo_fit)
mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(grepl("cell_group", term) & q_value < 0.05)
mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(gene_short_name %in% c("Pparg", "Abca1", "Abcg1", "Csf2ra"))
mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value) %>% dplyr::filter(term != "(Intercept)" & gene_short_name %in% c("Lpl", "Siglecf", "Apoe", "Car4", "F13a1", "Fabp4", "Fabp5"))

mac_gene_modules = monocle3::find_gene_modules(heathy_vs_pap_macs)
qplot(dim_1, dim_2, color=module, data=mac_gene_modules)


# Some tests for genes we know should be DE:
plot_genes_violin(heathy_vs_pap_macs[rowData(heathy_vs_pap_macs)$gene_short_name %in% c("Pparg", "Abca1", "Abcg1")], group_cells_by = "cell_group")

xxx_df = data.frame(expression = counts(heathy_vs_pap_macs)[rowData(heathy_vs_pap_macs)$gene_short_name %in% c("Csf2rb"),],
           weight = colData(heathy_vs_pap_macs)$num_cells_in_group,
           cell_group = colData(heathy_vs_pap_macs)$cell_group,
           size_factors = size_factors(heathy_vs_pap_macs))
qplot(cell_group, expression, size=weight, data=xxx_df, log="y")
qplot(cell_group, expression/size_factors, size=weight, data=xxx_df, log="y")
qplot(size_factors, weight, data=xxx_df, log="xy")


test_gene_cds = heathy_vs_pap_macs[rowData(heathy_vs_pap_macs)$gene_short_name %in% c("Pparg", "Abca1", "Abcg1")]
mac_pseudo_fit = fit_models(test_gene_cds,
                            "~cell_group")
mac_pseudo_coefs = coefficient_table(mac_pseudo_fit)
without_weights = mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect, q_value)

mac_pseudo_fit = fit_models(test_gene_cds,
                            "~cell_group",
                            weights=colData(heathy_vs_pap_macs)$num_cells_in_group)
mac_pseudo_coefs = coefficient_table(mac_pseudo_fit)
with_weights = mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect,q_value)


mac_pseudo_fit = fit_models(test_gene_cds,
                            "~cell_group",
                            weights=size_factors(heathy_vs_pap_macs))
mac_pseudo_coefs = coefficient_table(mac_pseudo_fit)
with_sfs = mac_pseudo_coefs %>% dplyr::select(gene_short_name, term, estimate, std_err, normalized_effect,q_value)

##########

gene_set = msigdbr(species = "Mus musculus", subcategory = "GO:MF")
transcription_regulators = gene_set %>%
  dplyr::select(gs_id, gene_symbol, gs_name) %>%
  dplyr::filter(grepl("Transcription", gs_name, ignore.case=TRUE)) %>%
  pull(gene_symbol) %>% unique
all_gene_modules$is_TF = all_gene_modules$gene_short_name %in% transcription_regulators

qplot(dim_1, dim_2, color=module, data=all_gene_modules) + facet_wrap(~is_TF) + theme(legend.position = "none")

mac_fit = fit_models(heathy_vs_pap_macs[rowData(heathy_vs_pap_macs)$gene_short_name %in% transcription_regulators], "~cell_group + num_cells_in_group")
mac_coefs = coefficient_table(mac_fit)
mac_coefs %>% dplyr::select(gene_short_name, term, estimate, p_value) %>% dplyr::filter(grepl("cell_group", term) & p_value < 0.05)


##### Cell-level DEG comparison between healthy and PAP macs:

mac_cds = ccs@cds[,colData(ccs@cds)$CW_assignedCellType %in% c("Healthy alveolar macrophages","hPAP alveolar macrophages")]

mac_cds_fit = fit_models(mac_cds, "~CW_assignedCellType + sampleName + Age + Genotype", cores=8)
mac_cds_coefs = coefficient_table(mac_cds_fit)
mac_cds_coefs %>% dplyr::select(gene_short_name, term, estimate, q_value) %>% dplyr::filter(grepl("CW_assignedCellType", term) & q_value < 0.05)




