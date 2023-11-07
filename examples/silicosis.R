cds = readRDS("silicosis_cds.cds")

# for simplicity, we are lumping together pre and post i.t. silica
colData(cds)$exposed = ifelse(colData(cds)$Timepoint == 0, "not exposed", "exposed")
colData(cds)$Rep = as.factor(colData(cds)$Rep)

plot_cells(cds, color_cells_by = "fine_annotation", label_groups_by_cluster = F)

ccs = new_cell_count_set(cds,
                         sample_group = "ID",
                         cell_group = "fine_annotation")

ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~ exposed")


cond_exp = estimate_abundances(ccm, tibble::tibble(exposed = "exposed"))
cond_not_exp = estimate_abundances(ccm, tibble::tibble(exposed = "not exposed"))

cond_ne_v_e_tbl = compare_abundances(ccm, cond_not_exp, cond_exp)

cond_ne_v_e_tbl %>% select(cell_group, perturbation_x, perturbation_y,
                           delta_log_abund, delta_log_abund_se, delta_q_value)

plot_contrast(ccm, cond_ne_v_e_tbl, q_value_threshold = 0.05)


# controlling for batch

ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~ exposed",
                            nuisance_model_formula_str = "~ Rep")

cond_exp = estimate_abundances(ccm, tibble::tibble(exposed = "exposed", Rep = "3"))
cond_not_exp = estimate_abundances(ccm, tibble::tibble(exposed = "not exposed", Rep = "3"))
cond_ne_v_e_tbl = compare_abundances(ccm, cond_not_exp, cond_exp)

plot_contrast(ccm, cond_ne_v_e_tbl, q_value_thresh = 0.05)
