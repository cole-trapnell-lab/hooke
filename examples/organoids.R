library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(monocle3)
library(googledrive)
library(splines)


devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")
setwd("~/OneDrive/UW/Trapnell/hooke/examples/")


if (!file.exists("R_objects/processed_retina_brain.RDS")) {
  drive_download("https://drive.google.com/file/d/1MMcdH77Xepz1s4_gyNlaH1Ax_SbJokGs/", 
                 path = "R_objects/processed_retina_brain.RDS")
}


if (!file.exists("R_objects/211221_pulse-screen_celltype_annotated_from_combined.RDS")) {
  drive_download("https://drive.google.com/file/d/1ZLRrOXFEI_3yC_Tmzn4ZDspTUPrNuGGV/", 
                 path = "R_objects/211221_pulse-screen_celltype_annotated_from_combined.RDS")
}


# trying out hooke for amy's data 

# ------------------------------------------------------------------------------

pulse_cds <- readRDS("R_objects/211221_pulse-screen_celltype_annotated_from_combined.RDS")
colData(pulse_cds)$sample = NULL
colData(pulse_cds)$treatment = paste0(pulse_cds@colData$molecule, " ", pulse_cds@colData$dose, "uM")

colData(pulse_cds)$treatment = pulse_cds@colData %>% as.data.frame() %>%
  mutate(treatment = case_when(grepl("5", days_on) ~ gsub("Ascl1", "Ascl1 OE", treatment),
                               grepl("0", days_on) ~ gsub("Ascl1", "No Ascl1", treatment),
                               TRUE ~ treatment)) %>%
  pull(treatment)


pulse_ccs = new_cell_count_set(pulse_cds,
                               sample_group = "top_oligo",
                               cell_group = "cell_type")

# colData(pulse_cds)$molecule %>% unique()
# "SAHA"                      "CPI 203"                   "Ascl1"                     "Purmorphamine"            
# [5] "RepSox"                    "I-BET 151 dihydrochloride" "LY 294002 hydrochloride"   "O4I2"  


pulse_ccm = new_cell_count_model(pulse_ccs, 
                                 model_formula_str = "~ molecule")


cond_repsox = estimate_abundances(pulse_ccm, tibble::tibble(molecule="RepSox"))
cond_cpi = estimate_abundances(pulse_ccm, tibble::tibble(molecule="CPI 203"))
cond_ibet = estimate_abundances(pulse_ccm, tibble::tibble(molecule="I-BET 151 dihydrochloride"))
cond_ly = estimate_abundances(pulse_ccm, tibble::tibble(molecule="LY 294002 hydrochloride"))
cond_purm = estimate_abundances(pulse_ccm, tibble::tibble(molecule="Purmorphamine"))
cond_O4I2 = estimate_abundances(pulse_ccm, tibble::tibble(molecule="O4I2"))
cond_ascl1 = estimate_abundances(pulse_ccm, tibble::tibble(molecule="Ascl1"))


# plot_contrast(pulse_ccm, 
#               cond_sb_vs_dmso_tbl, 
#               scale_shifts_by="none", 
#               plot_edges = F, 
#               plot_labels = "significant", 
#               q_value_thresh = 0.05)

# ------------------------------------------------------------------------------

retina_cds <- readRDS("R_objects/processed_retina_brain.RDS")
colData(retina_cds)$sample = NULL
retina_ccs = new_cell_count_set(retina_cds,
                   sample_group = "top_oligo",
                   cell_group = "cell_type")


colData(retina_ccs)$BMP = ifelse(colData(retina_ccs)$treatment == "BMP", T, F)
my_plot_cells(retina_ccs, color_cells_by="BMP")

colData(retina_cds)$treatment %>% unique()

retina_ccm = new_cell_count_model(retina_ccs, 
                                  model_formula_str = "~ treatment ")


cond_ct = estimate_abundances(retina_ccm, tibble::tibble(treatment="CONTROL"))
cond_bs = estimate_abundances(retina_ccm, tibble::tibble(treatment="BMP:SU5402"))
cond_bsc = estimate_abundances(retina_ccm, tibble::tibble(treatment="BMP:BMP:SU5402:CHIR"))
cond_bmp = estimate_abundances(retina_ccm, tibble::tibble(treatment="BMP"))
cond_bc = estimate_abundances(retina_ccm, tibble::tibble(treatment="BMP:CHIR"))
cond_scbbs = estimate_abundances(retina_ccm, tibble::tibble(treatment="SAG:CHIR:BMP:BMP:SU5402"))
cond_sag = estimate_abundances(retina_ccm, tibble::tibble(treatment="SAG"))


# "BMP:SU5402"
cond_bs_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_bs)

plot_contrast(retina_ccm, 
              cond_bs_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

# "BMP:BMP:SU5402:CHIR"
cond_bsc_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_bsc)

plot_contrast(retina_ccm, 
              cond_bsc_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

# "BMP"
cond_bmp_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_bmp)

plot_contrast(retina_ccm, 
              cond_bmp_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)
# "BMP:CHIR"

cond_bc_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_bc)

plot_contrast(retina_ccm, 
              cond_bc_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

# "SAG:CHIR:BMP:BMP:SU5402"

cond_scbbs_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_scbbs)

plot_contrast(retina_ccm, 
              cond_scbbs_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

# "SAG"

cond_sag_vs_ct_tbl = compare_abundances(retina_ccm, cond_ct, cond_sag)

plot_contrast(retina_ccm, 
              cond_sag_vs_ct_tbl, 
              scale_shifts_by="none", 
              plot_edges = F, 
              plot_labels = "significant", 
              q_value_thresh = 0.05)

# rbind(cond_bs_vs_ct_tbl$log_abund_y,
#       cond_bsc_vs_ct_tbl$log_abund_y,
#       cond_bmp_vs_ct_tbl$log_abund_y,
#       cond_bc_vs_ct_tbl$log_abund_y,
#       cond_scbbs_vs_ct_tbl$log_abund_y,
#       cond_sag_vs_ct_tbl$log_abund_y) %>% pheatmap::pheatmap()
