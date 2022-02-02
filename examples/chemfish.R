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

# together --------------------------------------------------------------------
chemfish_cds = chemfish_cds[,!is.na(colData(chemfish_cds)$major_group)]
chemfish_cds = chemfish_cds[,!is.na(colData(chemfish_cds)$cell_type_broad)]
cf_ccs = new_cell_count_set(chemfish_cds,
                            sample_group = "embryo",
                            cell_group = "cell_type_sub")

cf_ccm = new_cell_count_model(cf_ccs, 
                              model_formula_str = "~ splines::ns(timepoint, df=3) + drug")


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


cond_cp_36 = estimate_abundances(cf_ccm, tibble::tibble(drug="CP", timepoint = "36"))
cond_etoh_36 = estimate_abundances(cf_ccm, tibble::tibble(drug="EtOH", timepoint = "36"))
cond_cp_vs_etoh_tbl_36 = compare_abundances(cf_ccm, cond_etoh_36, cond_cp_36)

cond_cp_48 = estimate_abundances(cf_ccm, tibble::tibble(drug="CP", timepoint = "48"))
cond_etoh_48 = estimate_abundances(cf_ccm, tibble::tibble(drug="EtOH", timepoint = "48"))
cond_cp_vs_etoh_tbl_48 = compare_abundances(cf_ccm, cond_etoh_48, cond_cp_48)


cond_cp_72 = estimate_abundances(cf_ccm, tibble::tibble(drug="CP", timepoint = "72"))
cond_etoh_72 = estimate_abundances(cf_ccm, tibble::tibble(drug="EtOH", timepoint = "72"))
cond_cp_vs_etoh_tbl_72 = compare_abundances(cf_ccm, cond_etoh_72, cond_cp_72)

rbind(cond_cp_vs_etoh_tbl_36,
      cond_cp_vs_etoh_tbl_48,
      cond_cp_vs_etoh_tbl_72) %>% 
  filter(cell_group == "pronephric podocyte")


plot_contrast(cf_ccm,
              cond_cp_vs_etoh_tbl_36,
              scale_shifts_by="none",
              plot_labels = "significant",
              plot_edges = T,
              q_value_thresh = 0.05)


cond_cp_48 = estimate_abundances(cp_ccm, tibble::tibble(drug="CP", timepoint = "48"))
cond_etoh_48 = estimate_abundances(cp_ccm, tibble::tibble(drug="EtOH", timepoint = "48"))

cond_cp_vs_etoh_tbl_48 = compare_abundances(cp_ccm, cond_etoh_48, cond_cp_48)

plot_contrast(cf_ccm,
              cond_cp_vs_etoh_tbl_48,
              scale_shifts_by="none",
              plot_labels = "significant",
              plot_edges = T,
              q_value_thresh = 0.05)

cond_cp_72 = estimate_abundances(cp_ccm, tibble::tibble(drug="CP", timepoint = "72"))
cond_etoh_72 = estimate_abundances(cp_ccm, tibble::tibble(drug="EtOH", timepoint = "72"))

cond_cp_vs_etoh_tbl_72 = compare_abundances(cp_ccm, cond_etoh_72, cond_cp_72)

plot_contrast(cf_ccm,
              cond_cp_vs_etoh_tbl_72,
              scale_shifts_by="none",
              plot_labels = "significant",
              plot_edges = T,
              q_value_thresh = 0.05)


rbind(cond_cp_vs_etoh_tbl_36,
      cond_cp_vs_etoh_tbl_48,
      cond_cp_vs_etoh_tbl_72) %>% 
  filter(cell_group == "pronephric podocyte")


plot_sub_contrast(cf_ccm, 
                  cond_b_vs_a_tbl = cond_cp_vs_etoh_tbl_36, 
                  cell_group = "cell_type_sub",
                  scale_shifts_by="none", 
                  plot_labels = "significant", 
                  plot_edges = F, 
                  q_value_thresh = 0.05)

plot_sub_contrast(cf_ccm, 
                  cond_b_vs_a_tbl = cond_cp_vs_etoh_tbl_36, 
                  cell_group = "cell_type_sub",
                  select_group = "CNS",
                  scale_shifts_by="none", 
                  plot_labels = "significant", 
                  plot_edges = F, 
                  q_value_thresh = 0.05)


# filter(cell_group == "xanthophore") 