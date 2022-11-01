

pap_cds = readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pap.al_2021-09-15.cds.RDS")
pap_cds = detect_genes(pap_cds)
# assign best celltype column and reduce dims
colData(pap_cds)$cell_type = colData(pap_cds)$CW_assignedCellType
colData(pap_cds)$cluster = monocle3::clusters(pap_cds)
colData(pap_cds)$Size_Factor = size_factors(pap_cds)
colData(pap_cds)$experiment = colData(pap_cds)$sample
colData(pap_cds)$sample = NULL

ccs = new_cell_count_set(pap_cds,
                         sample_group = "sampleName",
                         cell_group = "cell_type")

ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "Genotype",
                            nuisance_model_formula_str = "batch", 
                            vhat_method = "wald")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=5)

cond_csf2ra = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_wt = estimate_abundances(ccm, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_csf2rb = estimate_abundances(ccm, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_ra_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2ra)
cond_rb_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_csf2rb)

plot_contrast(ccm, cond_ra_vs_wt_tbl, 0.01)
plot_contrast(ccm, cond_rb_vs_wt_tbl, 0.01)

cond_ra_vs_wt_tbl %>% filter(delta_q_value < 0.05) %>% nrow()
cond_rb_vs_wt_tbl %>% filter(delta_q_value < 0.05) %>% nrow()

# boot ------------------------------------------------------------------------

# debug(new_cell_count_model)
ccm_boot  = new_cell_count_model(ccs,
                            main_model_formula_str = "Genotype",
                            vhat_method = "bootstrap",
                            num_bootstraps = 100,
                            nuisance_model_formula_str = "batch")

ccm_boot = select_model(ccm_boot, criterion="StARS", sparsity_factor=5)

cond_boot_csf2ra = estimate_abundances(ccm_boot, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_boot_wt = estimate_abundances(ccm_boot, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_boot_csf2rb = estimate_abundances(ccm_boot, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_boot_ra_vs_wt_tbl = compare_abundances(ccm_boot, cond_boot_wt, cond_boot_csf2ra)
cond_boot_rb_vs_wt_tbl = compare_abundances(ccm_boot, cond_boot_wt, cond_boot_csf2rb)

cond_boot_ra_vs_wt_tbl %>% filter(delta_q_value < 0.05) #%>% nrow()
cond_boot_rb_vs_wt_tbl %>% filter(delta_q_value < 0.05) %>% nrow()

cond_boot_ra_vs_wt_tbl %>% filter(grepl("macrophages", cell_group))

plot_contrast(ccm_boot, cond_boot_ra_vs_wt_tbl, 0.01)
plot_contrast(ccm_boot, cond_boot_rb_vs_wt_tbl, 0.01)


# sandwich ---------------------------------------------------------------------

ccm_sand  = new_cell_count_model(ccs,
                                 main_model_formula_str = "Genotype",
                                 vhat_method = "sandwich",
                                 nuisance_model_formula_str = "batch")
ccm_sand = select_model(ccm_sand, criterion="StARS", sparsity_factor=5)

cond_sand_csf2ra = estimate_abundances(ccm_sand, tibble::tibble(Genotype="Csf2ra-/-", Age="12-13 weeks", batch="A", experiment=1))
cond_sand_wt = estimate_abundances(ccm_sand, tibble::tibble(Genotype="WT", Age="12-13 weeks", batch="A", experiment=1))
cond_sand_csf2rb = estimate_abundances(ccm_sand, tibble::tibble(Genotype="Csf2rb-/-", Age="12-13 weeks", batch="A", experiment=1))

cond_sand_ra_vs_wt_tbl = compare_abundances(ccm_sand, cond_sand_wt, cond_sand_csf2ra)
cond_sand_rb_vs_wt_tbl = compare_abundances(ccm_sand, cond_sand_wt, cond_sand_csf2rb)

cond_sand_ra_vs_wt_tbl %>% filter(delta_q_value < 0.05) %>% nrow()
cond_sand_rb_vs_wt_tbl %>% filter(delta_q_value < 0.05) %>% nrow()

plot_contrast(ccm_sand, cond_sand_ra_vs_wt_tbl)
plot_contrast(ccm, cond_sand_rb_vs_wt_tbl)


# compare values ---------------------------------------------------------------

rbind(cond_ra_vs_wt_tbl %>% mutate(method = "wald"), 
      cond_boot_ra_vs_wt_tbl %>% mutate(method = "bootstrap"), 
      cond_sand_ra_vs_wt_tbl %>% mutate(method = "sandwich")) %>% 
  select(cell_group, delta_log_abund, method) %>%
  pivot_wider(names_from = method, values_from = c(delta_log_abund, delta_q_value))

rbind(cond_rb_vs_wt_tbl %>% mutate(method = "wald"), 
      cond_boot_rb_vs_wt_tbl %>% mutate(method = "bootstrap"), 
      cond_sand_rb_vs_wt_tbl %>% mutate(method = "sandwich"))


# what is predict cond mean 


data(trichoptera)
trichoptera_prep <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLN(Abundance ~ Temperature + Wind, trichoptera_prep)


#Condition on the set of the first two species in the dataset (Hym, Hys) at the ten first sites
Yc <- trichoptera$Abundance[1:10, c(1, 2,3), drop=FALSE]
newX <- cbind(1, trichoptera$Covariate[1:10, c("Temperature", "Wind")])
pred <- predict_cond(myPLN, newX, Yc, type = "response")

pred[1:5,]
