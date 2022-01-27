library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(monocle3)
library(googledrive)
library(splines)


devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")
setwd("~/OneDrive/UW/Trapnell/hooke/examples/")

# devtools::load_all("~insert_path_here/hooke-master/")


if (!file.exists("R_objects/cds_metadata_proj_eliza.RDS")) {
  drive_download("https://drive.google.com/file/d/1lnCtTZig1inhu3KCNOjjqIUp7aZnOJNx/", 
                 path = "R_objects/cds_metadata_proj_eliza.RDS")
}


chemfish_cds = readRDS("R_objects/cds_metadata_proj_eliza.RDS")

# colData(chemfish_cds)$drug %>% unique()
# [1] "DMSO" "CP"   "EtOH" "SB" 

# together --------------------------------------------------------------------
chemfish_cds = chemfish_cds[,!is.na(colData(chemfish_cds)$major_group)]
chemfish_cds = chemfish_cds[,!is.na(colData(chemfish_cds)$cell_type_broad)]
cf_ccs = new_cell_count_set(chemfish_cds,
                            sample_group = "embryo",
                            cell_group = "cell_type_broad")

cf_ccm = new_cell_count_model(cf_ccs, 
                              model_formula_str = "~ drug + splines::ns(timepoint, df=3) ")


cond_sb = estimate_abundances(cf_ccm, tibble::tibble(drug="SB", timepoint = "24"))
cond_dmso = estimate_abundances(cf_ccm, tibble::tibble(drug="DMSO", timepoint = "24"))
cond_sb_vs_dmso_tbl = compare_abundances(cf_ccm, cond_dmso, cond_sb)

plot_contrast(cf_ccm, cond_sb_vs_dmso_tbl, scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

cond_cp = estimate_abundances(cf_ccm, tibble::tibble(drug="CP", timepoint = "24"))
cond_etoh = estimate_abundances(cf_ccm, tibble::tibble(drug="EtOH",timepoint = "24"))
cond_cp_vs_etoh_tbl = compare_abundances(cf_ccm, cond_etoh, cond_cp)


plot_contrast(cf_ccm, 
              cond_cp_vs_etoh_tbl, 
              scale_shifts_by="none", 
              plot_labels = "significant", 
              plot_edges = T, 
              q_value_thresh = 0.05) 





# figure out how to do the normalization? 
# normalize by timepoint + drug ? 

# SB --------------------------------------------------------------------------
sb_cds = chemfish_cds[,colData(chemfish_cds)$drug %in% c("SB", "DMSO")]
sb_cds = sb_cds[, !is.na(colData(sb_cds)$cell_type_broad)]
sb_ccs = new_cell_count_set(sb_cds,
                            sample_group = "embryo",
                            cell_group = "cell_type_broad")


sb_ccm = new_cell_count_model(sb_ccs, 
                              model_formula_str = "~ drug + as.factor(timepoint)")

# 24, 36, 60
cond_sb = estimate_abundances(sb_ccm, tibble::tibble(drug="SB", timepoint ="60"))
cond_dmso = estimate_abundances(sb_ccm, tibble::tibble(drug="DMSO", timepoint ="60"))

cond_sb_vs_dmso_tbl = compare_abundances(sb_ccm, cond_dmso, cond_sb)

cond_sb_vs_dmso_tbl %>% 
  filter(delta_q_value < 0.05) %>% 
  pull(cell_group)


plot_contrast(sb_ccm, cond_sb_vs_dmso_tbl, 
              scale_shifts_by="none", 
              plot_labels = "significant", 
              plot_edges = F, 
              q_value_thresh = 0.05)



# CP --------------------------------------------------------------------------
cp_cds = chemfish_cds[,colData(chemfish_cds)$drug %in%  c("CP", "EtOH")]
cp_cds = cp_cds[, !is.na(colData(cp_cds)$cell_type_broad)]
cp_ccs = new_cell_count_set(cp_cds,
                            sample_group = "embryo",
                            cell_group = "cell_type_broad")

cp_ccm = new_cell_count_model(cp_ccs, 
                              model_formula_str = "~ drug + as.factor(timepoint)")

# 36, 48, 72
cond_cp = estimate_abundances(cp_ccm, tibble::tibble(drug="CP", timepoint = "36"))
cond_etoh = estimate_abundances(cp_ccm, tibble::tibble(drug="EtOH", timepoint = "36"))

cond_cp_vs_etoh_tbl = compare_abundances(cp_ccm, cond_etoh, cond_cp)

cond_cp_vs_etoh_tbl %>% 
  filter(delta_q_value < 0.05) %>% 
  pull(cell_group)

plot_contrast(cp_ccm, 
              cond_cp_vs_etoh_tbl, 
              scale_shifts_by="none", 
              plot_labels = "significant",
              plot_edges = T, 
              q_value_thresh = 0.05)







