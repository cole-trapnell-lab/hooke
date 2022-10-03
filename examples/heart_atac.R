library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

#pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

heart_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/heart_atac/PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ArchRActivharmonyAligned_cds.rds")

sample_metadata = colData(heart_cds)$sampleName
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.heart", "")
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.Left.Vent", "\\.LV")
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.Right.Vent", "\\.RV")
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.s1", "")
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.apex", "\\.Apex")
sample_metadata = stringr::str_replace_all(sample_metadata, "\\.septum", "\\.Septum")
colData(heart_cds)$sampleName = sample_metadata

sample_metadata_cols = stringr::str_split_fixed(sample_metadata, "\\.", 2)
colData(heart_cds)$Anatomical_Site = sample_metadata_cols[,2]
colData(heart_cds)$Donor = sample_metadata_cols[,1]


hardAssignDonorAges <- function(inputCDS){
  colData(inputCDS)$Age = 0
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W134", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W135", 60, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W136", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W137", 49, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W139", 45, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W142", 55, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W144", 53, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W145", 51, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W146", 25, colData(inputCDS)$Age)

  colData(inputCDS)$Log10Age = log10(colData(inputCDS)$Age)

  return(inputCDS)
}

hardAssignDonorSexes <- function(inputCDS){
  colData(inputCDS)$Sex = "Not_Set"
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W134", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W135", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W136", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W137", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W139", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W142", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W144", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W145", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W146", "F", colData(inputCDS)$Sex)

  return(inputCDS)
}

heart_cds = hardAssignDonorSexes(heart_cds)
heart_cds = hardAssignDonorAges(heart_cds)

#heart_cds = cluster_cells(heart_cds)

plot_cells(heart_cds, color_cells_by="harmonyKNN_type")


# assign best celltype column and reduce dims
colData(heart_cds)$cell_type = colData(heart_cds)$harmonyKNN_type
#colData(heart_cds)$cluster = monocle3::clusters(heart_cds)

colData(heart_cds)$Size_Factor = size_factors(heart_cds)


colData(heart_cds)$tissue_sample = colData(heart_cds)$sampleName
colData(heart_cds)$sample = NULL

plot_cells(heart_cds, color_cells_by="cell_type", show_trajectory_graph=FALSE)

plot_cells(heart_cds, color_cells_by="tech", show_trajectory_graph=FALSE)


plot_cells(heart_cds, color_cells_by="Anatomical_Site", show_trajectory_graph=FALSE) + facet_wrap(~Anatomical_Site)

ccs = new_cell_count_set(heart_cds,
                         sample_group = "sampleName",
                         cell_group = "cell_type")


ccm  = new_cell_count_model(ccs,
                            model_formula_str = "~tech * (Age + Sex)")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=1)


cond_atac_lv = estimate_abundances(ccm, tibble::tibble(tech="ATAC", Anatomical_Site="LV", Age=55, Sex="M"))
cond_rna_lv = estimate_abundances(ccm, tibble::tibble(tech="RNA", Anatomical_Site="LV", Age=55, Sex="M"))

cond_atac_vs_rna_lv = compare_abundances(ccm, cond_atac_lv, cond_rna_lv)

plot_contrast(ccm, cond_atac_vs_rna_lv, scale_shifts_by="none", p_value_thresh=1)


# cond_atac_lv = estimate_abundances(ccm, tibble::tibble(tech="ATAC", Anatomical_Site="LV"))
# cond_atac_septum = estimate_abundances(ccm, tibble::tibble(tech="ATAC", Anatomical_Site="Septum"))
#
# cond_atac_lv_vs_septum = compare_abundances(ccm, cond_atac_lv, cond_atac_septum)
#
# plot_contrast(ccm, cond_atac_lv_vs_septum, scale_shifts_by="none", p_value_thresh=1)

cond_rna_lv_m_young = estimate_abundances(ccm, tibble::tibble(tech="RNA", Anatomical_Site="LV", Age=25, Sex="M"))
cond_rna_lv_m_old = estimate_abundances(ccm, tibble::tibble(tech="RNA", Anatomical_Site="LV", Age=60, Sex="M"))

cond_rna_lv_old_vs_young = compare_abundances(ccm, cond_rna_lv_m_young, cond_rna_lv_m_old)

plot_contrast(ccm, cond_rna_lv_old_vs_young, scale_shifts_by="none", p_value_thresh=1)

cond_rna_lv_m_old = estimate_abundances(ccm, tibble::tibble(tech="RNA", Anatomical_Site="LV", Age=60, Sex="M"))
cond_rna_lv_f_old = estimate_abundances(ccm, tibble::tibble(tech="RNA", Anatomical_Site="LV", Age=60, Sex="F"))

cond_rna_lv_m_vs_f = compare_abundances(ccm, cond_rna_lv_m_old, cond_rna_lv_f_old)

plot_contrast(ccm, cond_rna_lv_m_vs_f, scale_shifts_by="none", p_value_thresh=1)


cond_atac_lv_m_old = estimate_abundances(ccm, tibble::tibble(tech="ATAC", Anatomical_Site="LV", Age=60, Sex="M"))
cond_atac_lv_f_old = estimate_abundances(ccm, tibble::tibble(tech="ATAC", Anatomical_Site="LV", Age=60, Sex="F"))

cond_atac_lv_m_vs_f = compare_abundances(ccm, cond_atac_lv_m_old, cond_atac_lv_f_old)

plot_contrast(ccm, cond_atac_lv_m_vs_f, scale_shifts_by="none", p_value_thresh=0.05)



