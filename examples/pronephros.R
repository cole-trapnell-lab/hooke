library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

#kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

kidney_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/kidney/kidney.cds.cole.RDS")

kidney_cds = detect_genes(kidney_cds)

# assign best celltype column and reduce dims
colData(kidney_cds)$cell_type = colData(kidney_cds)$kidney.celltype
colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
colData(kidney_cds)$genotype = colData(kidney_cds)$gene_target1
colData(kidney_cds)$genotype[colData(kidney_cds)$genotype == "ctrl"] = "wt"
colData(kidney_cds)$Size_Factor = size_factors(kidney_cds)
colData(kidney_cds)$cell_type_dk = case_when(
  colData(kidney_cds)$cluster == 1 ~ "Proximal Convoluted Tubule (1)",
  colData(kidney_cds)$cluster == 2 ~ "Late Distal (late)",
  colData(kidney_cds)$cluster == 3 ~ "Distal Early",
  colData(kidney_cds)$cluster == 4 ~ "Early Distal (late)",
  colData(kidney_cds)$cluster == 5 ~ "Proximal Straight Tubule",
  colData(kidney_cds)$cluster == 6 ~ "Early podocyte",
  colData(kidney_cds)$cluster == 7 ~ "Early duct",
  colData(kidney_cds)$cluster == 8 ~ "Cloaca",
  colData(kidney_cds)$cluster == 9 ~ "Late neck",
  colData(kidney_cds)$cluster == 10 ~ "Proximal Convoluted Tubule (10)",
  colData(kidney_cds)$cluster == 11 ~ "Podocyte",
  colData(kidney_cds)$cluster == 12 ~ "Corpuscles of Stannius",
  colData(kidney_cds)$cluster == 13 ~ "Early Neck",
  colData(kidney_cds)$cluster == 14 ~ "Multiciliated cells",
  colData(kidney_cds)$cluster == 15 ~ "Proximal Convoluted Tubule (15)",
  TRUE ~ "Unknown"
)

colData(kidney_cds)$experiment = colData(kidney_cds)$expt
colData(kidney_cds)$sample = NULL

plot_cells(kidney_cds, color_cells_by="cluster", show_trajectory_graph=FALSE)


plot_cells(kidney_cds, color_cells_by="cell_type_dk", show_trajectory_graph=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_cell_types.png", width=3, height=3)


plot_cells(kidney_cds, color_cells_by="gene_target1", show_trajectory_graph=FALSE) + facet_wrap(~gene_target1)

#plot_cells(cds)

kidney_cds = kidney_cds[,is.na(colData(kidney_cds)$Oligo) == FALSE & is.na(colData(kidney_cds)$timepoint.1) == FALSE & colData(kidney_cds)$timepoint <= 48]


plot_cells(kidney_cds, color_cells_by="timepoint.1", show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #theme(legend.position="none") +
  ggsave("kidney_time.png", width=3, height=3)


noto_smo_cds = kidney_cds[,colData(kidney_cds)$genotype %in% c("wt", "noto", "noto-mut", "smo")]


noto_smo_gap16_cds = kidney_cds[,colData(kidney_cds)$genotype %in% c("wt", "noto", "smo", "noto-mut") &
                                  colData(kidney_cds)$experiment == "GAP16" &
                                  colData(kidney_cds)$timepoint %in% c(18,24,36)]

gap16_cds = kidney_cds[,colData(kidney_cds)$experiment == "GAP16"]



plot_cells(noto_smo_gap16_cds, color_cells_by="timepoint.1", show_trajectory_graph=FALSE, label_cell_groups=FALSE) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  ggsave("kidney_time.png", width=3, height=3)



#subset_kidney_cds = kidney_cds[,colData(kidney_cds)$genotype]

ccs = new_cell_count_set(kidney_cds,
                         sample_group = "Oligo",
                         cell_group = "cell_type_dk")


ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~splines::ns(timepoint,df=3) + genotype",
                            nuisance_model_formula_str = "experiment")

ccm = select_model(ccm, criterion="StARS", sparsity_factor=0.1)


cond_wt_18 = estimate_abundances(ccm, tibble::tibble(timepoint=18, genotype="wt", experiment="GAP16"))
cond_wt_24 = estimate_abundances(ccm, tibble::tibble(timepoint=24, genotype="wt", experiment="GAP16"))
cond_wt_36 = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="wt", experiment="GAP16"))
cond_18_vs_24_wt = compare_abundances(ccm, cond_wt_18, cond_wt_24)
plot_contrast(ccm, cond_18_vs_24_wt, scale_shifts_by="none", q_value_thresh=0.01)

cond_24_vs_36_wt = compare_abundances(ccm, cond_wt_24, cond_wt_36)
plot_contrast(ccm, cond_24_vs_36_wt, scale_shifts_by="none", q_value_thresh=0.01)


cond_smo_18 = estimate_abundances(ccm, tibble::tibble(timepoint=18, genotype="smo", experiment="GAP16"))

cond_smo = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="smo", experiment="GAP16"))
cond_wt = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="wt", experiment="GAP16"))
cond_noto = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="noto", experiment="GAP16"))
cond_noto_mut = estimate_abundances(ccm, tibble::tibble(timepoint=36,genotype="noto-mut", experiment="GAP16"))

#cond_tbxta = estimate_abundances(ccm, tibble::tibble(timepoint.1="36", genotype="tbxta", experiment="GAP16"))

cond_met = estimate_abundances(ccm, tibble::tibble(timepoint=36, genotype="met", experiment="GAP16"))


# cond_ctrl_vs_wt = compare_abundances(ccm, cond_wt, cond_ctrl_inj)
#
# plot_contrast(ccm, cond_ctrl_vs_wt, scale_shifts_by="none", p_value_thresh=0.05, plot_labels="none") +
#   theme_minimal() + monocle3:::monocle_theme_opts() +
#   ggsave("kidney_wt_vs_ctrl_inj.png", width=7, height=6)

cond_smo_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_smo)

plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01)


plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01, plot_labels="none", receiver_cell_groups=c("none")) +
  ggsave("smo_vs_wt_36h_fc.png", width=4.5, height=3)


#plot_contrast(ccm, cond_smo_vs_wt_tbl, scale_shifts_by="none", p_value_thresh=1, receiver_cell_groups = c("Neck Segment"))

cond_noto_vs_wt_tbl = compare_abundances(ccm, cond_wt, cond_noto)

plot_contrast(ccm, cond_noto_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01)


plot_contrast(ccm, cond_noto_vs_wt_tbl, scale_shifts_by="none", q_value_thresh=0.01, plot_labels="none", receiver_cell_groups=c("none"))+
  ggsave("noto_vs_wt_36h_fc.png", width=4.5, height=3)

#########
# Plot changes over time in wild type

cond_wt_timeseries = estimate_abundances(ccm, tibble::tibble(timepoint=seq(18,48),genotype="wt", experiment="GAP16"))



