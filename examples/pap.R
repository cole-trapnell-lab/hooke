library(ggplot2)
library(monocle3)
library(hooke)

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




