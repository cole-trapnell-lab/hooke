library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(garnett)
library(msigdbr)
library(tidyr)

# kidney_cds = readRDS("/Users/maddyduran/OneDrive/UW/Trapnell/gap-notebook-ct-1/R_objects/kidney.cds.cole.RDS")
# kidney_cds = detect_genes(kidney_cds)
#
# # assign best celltype column and reduce dims
# colData(kidney_cds)$cell_type = colData(kidney_cds)$kidney.celltype
# colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
# colData(kidney_cds)$genotype = colData(kidney_cds)$gene_target1
# colData(kidney_cds)$genotype[colData(kidney_cds)$genotype == "ctrl"] = "wt"
# colData(kidney_cds)$Size_Factor = size_factors(kidney_cds)
# colData(kidney_cds)$cell_type_dk = case_when(
#   colData(kidney_cds)$cluster == 1 ~ "Proximal Convoluted Tubule (1)",
#   colData(kidney_cds)$cluster == 2 ~ "Mature Distal late",
#   colData(kidney_cds)$cluster == 3 ~ "Immature Proximal Tubule" , #Proximal Convoluted Tubule (1)",
#   colData(kidney_cds)$cluster == 4 ~ "Mature Distal Early",
#   colData(kidney_cds)$cluster == 5 ~ "Immature Distal late", #
#   colData(kidney_cds)$cluster == 6 ~ "Proximal Straight Tubule", #"Early duct",
#   colData(kidney_cds)$cluster == 7 ~ "Immature podocyte",
#   colData(kidney_cds)$cluster == 8 ~ "Cloaca",
#   colData(kidney_cds)$cluster == 9 ~ "Mature neck",
#   colData(kidney_cds)$cluster == 10 ~ "Proximal Convoluted Tubule (10)",
#   colData(kidney_cds)$cluster == 11 ~ "Podocyte",
#   colData(kidney_cds)$cluster == 12 ~ "Immature Neck", #
#   colData(kidney_cds)$cluster == 13 ~ "Corpuscles of Stannius",
#   colData(kidney_cds)$cluster == 14 ~ "Multiciliated cells",
#   colData(kidney_cds)$cluster == 15 ~ "Proximal Convoluted Tubule (15)",
#   colData(kidney_cds)$cluster == 16 ~ "Unknown",
#   TRUE ~ "Unknown"
# )
# kidney_cds = cluster_cells(kidney_cds, resolution=1e-4)
#
# colData(kidney_cds)$cluster = monocle3::clusters(kidney_cds)
# colData(kidney_cds)$cell_type_ct = case_when(
#   colData(kidney_cds)$cluster %in% c(12,1,3,5,26,21) ~ "Proximal Convoluted Tubule",
#   colData(kidney_cds)$cluster %in% c(7,16,13) ~ "Distal Early",
#   colData(kidney_cds)$cluster %in% c(6,9,4) ~ "Distal Late", #
#   colData(kidney_cds)$cluster %in% c(2,25) ~ "Proximal Straight Tubule", #"Early duct",
#   colData(kidney_cds)$cluster %in% c(29,14) ~ "Cloaca",
#   colData(kidney_cds)$cluster %in% c(19,31,23,11) ~ "Podocyte",
#   colData(kidney_cds)$cluster %in% c(36,22,39,8) ~ "Neck", #
#   colData(kidney_cds)$cluster %in% c(20) ~ "Corpuscles of Stannius",
#   colData(kidney_cds)$cluster %in% c(28) ~ "Multiciliated cells",
#   colData(kidney_cds)$cluster %in% c(10,15,27) ~ "Renal progenitors",
#   #colData(kidney_cds)$cluster == 16 ~ "Unknown",
#   TRUE ~ "Unknown"
# )
#
# pronephros_classifier <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pronephros_classifier.RDS")
#
# kidney_cds = classify_cells(kidney_cds, pronephros_classifier, db="none")
# #
# colData(kidney_cds)$segment = case_when(
#   colData(kidney_cds)$cell_type %in% c("Renal progenitors") ~ -1,
#   colData(kidney_cds)$cell_type %in% c("Podocyte") ~ 0,
#   colData(kidney_cds)$cell_type %in% c("Neck") ~ 1,
#   colData(kidney_cds)$cell_type %in% c("Proximal Convoluted Tubule") ~ 2,
#   colData(kidney_cds)$cell_type %in% c("Proximal Straight Tubule") ~ 3,
#   colData(kidney_cds)$cell_type %in% c("Distal Early") ~ 4,
#   colData(kidney_cds)$cell_type %in% c("Corpuscles of Stannius") ~ 5,
#   colData(kidney_cds)$cell_type %in% c("Distal Late") ~ 6,
#   colData(kidney_cds)$cell_type %in% c("Multiciliated cells")  ~ 7,
#   colData(kidney_cds)$cell_type %in% c("Cloaca")  ~ 8
# )
# colData(kidney_cds)$experiment = colData(kidney_cds)$expt
# colData(kidney_cds)$sample = NULL
# kidney_cds = kidney_cds[,is.na(colData(kidney_cds)$Oligo) == FALSE & is.na(colData(kidney_cds)$timepoint.1) == FALSE & colData(kidney_cds)$timepoint <= 48]

# saveRDS(kidney_cds, "~/OneDrive/UW/Trapnell/hooke/examples/R_objects/maddy_kidney_cds.rds")

# this has the similar clustering resolution as cole's i think
kidney_cds <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/maddy_kidney_cds.rds")


plot_cells(kidney_cds, color_cells_by = "cluster")

# -------

wt_cds = kidney_cds[,colData(kidney_cds)$gene_target %in% c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met") &
                      colData(kidney_cds)$experiment %in% c("GAP13", "GAP14", "GAP18", "HF4")  ]


colData(wt_cds)$cluster = clusters(wt_cds)


# cds from cole
# wt_cds <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/wt_cds.RDS")


wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "Oligo",
                            cell_group = "cluster")


wl = get_paga_graph(wt_ccs@cds) %>% igraph::as_data_frame()

wt_ccm  = new_cell_count_model(wt_ccs,
                               main_model_formula_str = "~splines::ns(timepoint,df=4)",
                               nuisance_model_formula_str = "~experiment",
                               whitelist = NULL )

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 main_model_formula_str = "~splines::ns(timepoint,df=4)",
                                 nuisance_model_formula_str = "~experiment",
                                 whitelist = wl )

wt_ccm_01 = select_model(wt_ccm, criterion = "EBIC", sparsity_factor=0.1)
wt_ccm_02 = select_model(wt_ccm, criterion = "EBIC", sparsity_factor=0.2)

wt_ccm_wl_01 = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.1)
wt_ccm_wl_02 = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=0.2)


get_timepoint_pred <- function(wt_ccm) {

  timepoint_pred_df = tibble(timepoint=seq(18,48))

  timepoint_pred_df = timepoint_pred_df %>%
    dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
      function(tp){ estimate_abundances(wt_ccm, tibble(timepoint=tp, experiment="GAP14"))}, NA_real_),
      .x = timepoint)) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  cell_type_assignments = colData(wt_ccm@ccs@cds) %>%
    as.data.frame %>%
    dplyr::count(cluster, cell_type) %>%
    group_by(cluster) %>% slice_max(n) %>%
    dplyr::select(cell_group=cluster, cell_type)

  timepoint_pred_df = left_join(timepoint_pred_df, cell_type_assignments)
  timepoint_pred_df = timepoint_pred_df %>% mutate(cell_group_label = paste(cell_type, " (", cell_group, ")", sep=""))

  return(timepoint_pred_df)

}


timepoint_pred = get_timepoint_pred(wt_ccm)
timepoint_pred_01 = get_timepoint_pred(wt_ccm_01)
timepoint_pred_02 = get_timepoint_pred(wt_ccm_02)

timepoint_pred_wl = get_timepoint_pred(wt_ccm_wl)
timepoint_pred_wl_01 = get_timepoint_pred(wt_ccm_wl_01)
timepoint_pred_wl_02 = get_timepoint_pred(wt_ccm_wl_02)


# -----------------------------------------------------------------------------


plot_contrast_wrapper <- function(wt_ccm, t1, t2) {

  timepoint_pred_df = get_timepoint_pred(wt_ccm)

  plot_contrast(wt_ccm, compare_abundances(wt_ccm,
                                           timepoint_pred_df %>% filter(timepoint == t1),
                                           timepoint_pred_df %>% filter(timepoint == t2)),
                q_value_thresh = 0.01)

}

t1 = 18
t2 = 30
plot_contrast_wrapper(wt_ccm, t1, t2)
plot_contrast_wrapper(wt_ccm_01, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_02, t1, t2)

plot_contrast_wrapper(wt_ccm_wl, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_01, t1, t2)
plot_contrast_wrapper(wt_ccm_wl_02, t1, t2)

# -----------------------------------------------------------------------------


find_union_path <- function(wt_ccm, q_val=0.01, time_diff = 4) {

  timepoints = seq(18, 48, 2)

  timepoint_pred_df = get_timepoint_pred(wt_ccm)

  select_timepoints <- function(timepoint_pred_df, t1, t2)  {
    cond_x = timepoint_pred_df %>% filter(timepoint == t1)
    cond_y = timepoint_pred_df %>% filter(timepoint == t2)
    return(compare_abundances(wt_ccm, cond_x, cond_y))
  }

  times = expand.grid("t1" = timepoints, "t2" = timepoints) %>%
    filter(t1 < t2, (t2-t1) >= time_diff) %>%
    mutate(comp_abund = purrr::map2(.f = select_timepoints,
                                    .x = t1,
                                    .y = t2,
                                    timepoint_pred_df = timepoint_pred_df)) %>%
    mutate(neg_rec_edges = purrr::map(.f = get_neg_dir_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(pos_rec_edges = purrr::map(.f = get_positive_edges,
                                      .x = comp_abund,
                                      ccm = wt_ccm,
                                      q_value_threshold = q_val)) %>%
    mutate(path = purrr::map(.f = purrr::possibly(get_path, NA_real_),
                             .x = comp_abund,
                             ccm = wt_ccm,
                             q_value_threshold = q_val))

  distinct_edges = times %>%
    select(neg_rec_edges) %>%
    tidyr::unnest(neg_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  pos_edges = times %>%
    select(pos_rec_edges) %>%
    tidyr::unnest(pos_rec_edges) %>%
    group_by(from, to) %>%
    summarise(pcor = abs(sum(pcor))) %>%
    mutate(scaled_weight = abs(pcor) / max(abs(pcor)))

  paths = times %>% select(path) %>%
    filter(!is.na(path)) %>%
    tidyr::unnest(path) %>%
    group_by(from, to) %>%
    tally() %>%
    mutate(scaled_weight = abs(n) / max(abs(n)))

  return(list(neg_rec_edges = distinct_edges,
              pos_edges = pos_edges,
              paths = paths))

}

results_00 = find_union_path(wt_ccm)
results_01 = find_union_path(wt_ccm_01)
results_02 = find_union_path(wt_ccm_02)

results_wl = find_union_path(wt_ccm_wl)
results_wl_01 = find_union_path(wt_ccm_wl_01)
results_wl_02 = find_union_path(wt_ccm_wl_02)


plot_path(wt_ccm, path_df = results_00$paths)
plot_path(wt_ccm, path_df = results_01$paths)
plot_path(wt_ccm, path_df = results_02$paths)

plot_path(wt_ccm_01, path_df = results_wl$paths)
plot_path(wt_ccm_01, path_df = results_wl_01$paths)
plot_path(wt_ccm_01, path_df = results_wl_02$paths)




# how many clusters are connected
union(results_wl$paths$to, results_wl$paths$from) %>% unique %>% length()

# ----------------------------------------------------------------------------



#gp = my_plot_cells(kidney_ccs, color_cells_by = "timepoint")

# edges = distinct_edges %>%
#   add_umap_coords(centroids(wt_ccs))
#
# edges = pos_edges %>%
#   add_umap_coords(centroids(wt_ccs))

my_plot_path <- function(wt_ccs, edges, x = 1, y = 2) {

  gp = my_plot_cells(wt_ccs, color_cells_by = "timepoint", x=x, y=y)

  edges = edges %>%
      add_umap_coords(centroids(wt_ccs))

  gp +
    geom_segment(data = edges,
                 aes(x = get(paste0("umap_to_", x)),
                     y = get(paste0("umap_to_", y)),
                     xend =  get(paste0("umap_from_", x)),
                     yend =  get(paste0("umap_from_", y)),
                     size = scaled_from_weight,
                 ),
                 color="black") +
    geom_segment(data = edges,
                 aes(x = get(paste0("umap_from_", x)),
                     y = get(paste0("umap_from_", y)),
                     xend = ( get(paste0("umap_to_", x))+ get(paste0("umap_from_", x)) )/2,
                     yend = ( get(paste0("umap_to_", y))+ get(paste0("umap_from_", y)) )/2,
                     size = scaled_from_weight,
                 ),
                 color="black",
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
    scale_size_identity()

  # return(p)

}



# gp = my_plot_cells(wt_ccs, color_cells_by = "timepoint")
my_plot_path(wt_ccs, results_wl$paths) %>%
  ggsave(filename = "kidney_wt_path.png")

my_plot_path(wt_ccs, results_wl$paths, x = 2, y = 3) %>%
  ggsave(filename = "kidney_wt_path_x2y3.png")

# my_plot_path(gp, results_wl$pos_edges)
# my_plot_path(gp, results_wl$neg_rec_edges)


# plot by segment -------------------------------------------------------------



space_time_summary = colData(wt_ccm_01@ccs@cds) %>% as.data.frame %>% group_by(cluster) %>% summarise(mean_hpf = mean(timepoint),
                                                                                                mean_segment = mean(segment, na.rm=TRUE))
cell_type_assignments = colData(wt_ccm_01@ccs@cds) %>% as.data.frame %>% dplyr::count(cluster, cell_type) %>% group_by(cluster) %>% slice_max(n)
space_time_summary = left_join(space_time_summary, cell_type_assignments)
space_time_summary$cell_group = space_time_summary$cluster
space_time_summary$cluster = NULL

gp = ggplot(aes(x=mean_hpf, y=mean_segment), data=space_time_summary)


directed_edge_df = results_wl$paths %>%
  left_join(space_time_summary, by = c("from" = "cell_group")) %>%
  dplyr::rename("from_mean_hpf" = "mean_hpf", "from_mean_segment" = "mean_segment", "from_cell_type" = cell_type) %>%
  left_join(space_time_summary, by = c("to" = "cell_group")) %>%
  dplyr::rename("to_mean_hpf" = "mean_hpf", "to_mean_segment" = "mean_segment", "to_cell_type" = cell_type)


label_df = space_time_summary #%>% filter(cell_group %in% union(directed_edge_df$from, directed_edge_df$to))

segment_plot_all_labels = ggplot() +
  geom_segment(data = directed_edge_df ,
               aes(x = from_mean_hpf,
                   y = from_mean_segment,
                   xend=to_mean_hpf,
                   yend = to_mean_segment),
               color="black") +
  geom_segment(data = directed_edge_df ,
               aes(x = from_mean_hpf,
                   y = from_mean_segment,
                   xend=(from_mean_hpf+to_mean_hpf)/2,
                   yend = (from_mean_segment+to_mean_segment)/2),
               color="black",
               linejoin='mitre',
               arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
  scale_size_identity() +
  ggrepel::geom_label_repel(data = label_df,
                                       mapping = aes(mean_hpf, mean_segment, label=paste(cell_group," - ", cell_type)),
                                       # size=I(group_label_size),
                                       fill = "white") +
  monocle3:::monocle_theme_opts()

segment_plot = ggplot() +
                geom_segment(data = directed_edge_df ,
                             aes(x = from_mean_hpf,
                                 y = from_mean_segment,
                                 xend=to_mean_hpf,
                                 yend = to_mean_segment),
                             color="black") +
                geom_segment(data = directed_edge_df ,
                             aes(x = from_mean_hpf,
                                 y = from_mean_segment,
                                 xend=(from_mean_hpf+to_mean_hpf)/2,
                                 yend = (from_mean_segment+to_mean_segment)/2),
                             color="black",
                             linejoin='mitre',
                             arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
                scale_size_identity() +
                ggrepel::geom_label_repel(data = space_time_summary %>%
                                            filter(cell_group %in% union(directed_edge_df$from, directed_edge_df$to)),
                                          mapping = aes(mean_hpf, mean_segment, label=paste(cell_group," - ", cell_type)),
                                          # size=I(group_label_size),
                                          fill = "white") +
                monocle3:::monocle_theme_opts()
segment_plot


# -----------------------------------------------------------------------------

kidney_ccs = new_cell_count_set(kidney_cds,
                                sample_group = "Oligo",
                                cell_group = "cluster")

my_plot_cells(kidney_ccs, color_cells_by = "cluster") %>%
  ggsave(filename = "kidney_cluster.png")

my_plot_cells(kidney_ccs, color_cells_by = "timepoint")# %>%
  ggsave(filename = "kidney_timepoint.png")


colData(kidney_ccs@cds)$segment = as.numeric(colData(kidney_ccs@cds)$segment)

cell_size = 1

colData(kidney_ccs@cds)$umap2D_1 = reducedDim(kidney_ccs@cds, type="UMAP")[,1]
colData(kidney_ccs@cds)$umap2D_2 = reducedDim(kidney_ccs@cds, type="UMAP")[,2]
plot_df = colData(kidney_ccs@cds) %>%
  as.data.frame()

gp = ggplot() +
  geom_point(data = plot_df,
    aes(umap2D_1, umap2D_2),
    color = "black",
    size = 1.5 * cell_size,
    stroke = 0
  ) +
  geom_point(data = plot_df,
    aes(umap2D_1, umap2D_2, color = as.numeric(segment)),
    size = cell_size,
    stroke = 0
  ) +
  theme_void() +
  viridis::scale_color_viridis(option = "C")+
  theme(legend.position = "right") +
  labs(color="segment") +
  monocle3:::monocle_theme_opts()

ggsave(gp, filename = "kidney_segment.png")

my_plot_cells(kidney_ccs, color_cells_by = "cell_type", cell_size = 4, legend_position = "right") #%>%
  ggsave(filename = "kidney_cell_type.png")

my_plot_cells(kidney_ccs, color_cells_by = "cell_type", x = 2, y = 3) %>%
  ggsave(filename = "kidney_cell_type_x2_y3.png")

my_plot_cells(kidney_ccs, color_cells_by = "cluster", x = 2, y = 3)

gp = my_plot_cells(kidney_ccs, color_cells_by = "aberant.noto.smo.cells",
              x = 2, y = 3) +
  scale_color_manual(values = c("royalblue3", "orangered3"))
ggsave(gp, filename = "kidney_aberant.png")

space_time_summary %>%
  filter(cell_group %in% c(1, 33, 5, 13, 4))


# mutant ----------------------------------------------------------------------

WT_whitelist = results_wl$paths %>% select(from, to)



ccs = new_cell_count_set(kidney_cds[,colData(kidney_cds)$experiment %in% c("GAP14", "GAP18", "GAP16")],
                         sample_group = "Oligo",
                         cell_group = "cluster")

fit_genotype_ccm = function(genotype,
                            ccs,
                            ctrl_ids=c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                            num_time_breaks=3){
  subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]

  colData(subset_ccs)$knockout = colData(subset_ccs)$gene_target == genotype
  knockout_time_start = min(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  knockout_time_stop = max(colData(subset_ccs)$timepoint[colData(subset_ccs)$knockout])
  subset_ccs = subset_ccs[,colData(subset_ccs)$timepoint >= knockout_time_start & colData(subset_ccs)$timepoint <= knockout_time_stop]
  time_breakpoints = c()
  if (num_time_breaks > 2 & knockout_time_stop > knockout_time_start){
    time_breakpoints = seq(knockout_time_start, knockout_time_stop, length.out=num_time_breaks)
    time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
    main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ") + knockout")
  }else{
    main_model_formula_str = "~ knockout"
  }

  if (length(unique(colData(subset_ccs)$experiment)) > 1)
    nuisance_model_formula_str = "~experiment"
  else
    nuisance_model_formula_str = "~1"

  # FIXME: Should prob use splines
  #  compound_ccm = new_cell_count_model(subset_ccs,
  #                                      main_model_formula_str = "~ splines::ns(compound_dose, knots=c(1,3))")
  genotype_ccm = suppressWarnings(new_cell_count_model(subset_ccs,
                                                       main_model_formula_str = main_model_formula_str,
                                                       #main_model_formula_str = "~ splines::ns(timepoint, knots=c(24, 30, 36)) + knockout",
                                                       #main_model_formula_str = "~ as.factor(timepoint) + knockout",
                                                       nuisance_model_formula_str = nuisance_model_formula_str,
                                                       whitelist = WT_whitelist
  ))
  genotype_ccm = select_model(genotype_ccm, sparsity_factor = 0.1)
  return(genotype_ccm)
}

genotype_df = colData(ccs@cds) %>% as_tibble() %>% dplyr::select(gene_target) %>% distinct()



crispant_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("smo", "noto", "tbxta", "cdx4", "egr2b", "mafba", "epha4a"))

crispant_models_tbl = crispant_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, ctrl_ids=c("wt", "ctrl-inj")))

noto_mut_genotype_df = dplyr::filter(genotype_df, gene_target %in% c("noto-mut"))
noto_models_tbl = noto_mut_genotype_df %>%
  dplyr::mutate(genotype_ccm = purrr::map(.f = purrr::possibly(
    fit_genotype_ccm, NA_real_), .x = gene_target, ccs, ctrl_ids=c("wt", "ctrl-noto")))

genotype_models_tbl = rbind(noto_models_tbl, crispant_models_tbl)

collect_genotype_effects = function(ccm, timepoint=24, experiment="GAP16"){
  control_abund = estimate_abundances(ccm, tibble(knockout=FALSE, timepoint=timepoint, experiment=experiment))
  knockout_abund = estimate_abundances(ccm, tibble(knockout=TRUE, timepoint=timepoint, experiment=experiment))
  genotype_comparison_tbl = compare_abundances(ccm, control_abund, knockout_abund)
}

genotype_models_tbl_18 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=18))

genotype_models_tbl_24 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=24))

genotype_models_tbl_30 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=30))

genotype_models_tbl_36 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=36))

genotype_models_tbl_42 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=42))

genotype_models_tbl_48 = genotype_models_tbl %>%
  dplyr::mutate(genotype_eff = purrr::map(.f = purrr::possibly(
    collect_genotype_effects, NA_real_), .x = genotype_ccm, timepoint=48))



collect_state_transitions = function(ccm, dose_contrast){
  transitions = hooke:::collect_pln_graph_edges(ccm, dose_contrast)
  transitions = transitions %>% dplyr::filter(edge_type != "undirected")
  return(transitions)
}



plot_contrast(genotype_models_tbl_24$genotype_ccm[[4]],
              genotype_models_tbl_24$genotype_eff[[4]], q_value_thresh = 0.01, x=2, y=3)

plot_contrast(genotype_models_tbl_24$genotype_ccm[[4]],
              genotype_models_tbl_24$genotype_eff[[4]],
              x = 2, y=3, q_value_thresh = 0.01)
ggsave(p, filename = "smo24_x2y3.png", width = 10, height=7)

# collect edges


get_neg_dir_edges(genotype_models_tbl_24$genotype_ccm[[4]],
                  genotype_models_tbl_24$genotype_eff[[4]],
                  0.01)

# cluster from 1 to 33
# cluster from 13 to 4


p = plot_contrast(genotype_models_tbl_24$genotype_ccm[[8]],
              genotype_models_tbl_24$genotype_eff[[8]],
              x = 2, y=3, q_value_thresh = 0.01)
ggsave(p, filename = "noto24_x2y3.png", width = 10, height=7)



get_neg_dir_edges(genotype_models_tbl_24$genotype_ccm[[8]],
                  genotype_models_tbl_24$genotype_eff[[8]],
                  0.01)


p = plot_contrast(genotype_models_tbl_24$genotype_ccm[[1]],
                  genotype_models_tbl_24$genotype_eff[[1]],
                  x = 2, y=3, q_value_thresh = 0.01)
ggsave(p, filename = "noto-mut_x2y3.png", width = 10, height=7)

p = plot_contrast(genotype_models_tbl_24$genotype_ccm[[7]],
                  genotype_models_tbl_24$genotype_eff[[7]],
                  x = 2, y=3, q_value_thresh = 0.01)
ggsave(p, filename = "cdx4_x2y3.png", width = 10, height=7)

plot_contrast(genotype_models_tbl_30$genotype_ccm[[4]],
              genotype_models_tbl_30$genotype_eff[[4]],
              x = 2, y=3, q_value_thresh = 0.01)

plot_contrast(genotype_models_tbl_30$genotype_ccm[[8]],
              genotype_models_tbl_30$genotype_eff[[8]],
              x = 2, y=3, q_value_thresh = 0.01)

plot_contrast(genotype_models_tbl_18$genotype_ccm[[4]],
              genotype_models_tbl_18$genotype_eff[[4]],
              x = 2, y=3, q_value_thresh = 0.01)

plot_contrast(genotype_models_tbl_18$genotype_ccm[[8]],
              genotype_models_tbl_18$genotype_eff[[8]],
              x = 2, y=3, q_value_thresh = 0.01)



# ------------------------------------------------------------------------------

gen_mat = genotype_models_tbl_24 %>%
  tidyr::unnest(genotype_eff) %>%
  select(gene_target, cell_group, delta_log_abund) %>%
  pivot_wider(names_from = "cell_group", values_from = delta_log_abund) %>%
  column_to_rownames("gene_target")

gen_mat %>% pheatmap::pheatmap()


genotype_paths_18 = genotype_models_tbl_18 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                            .x = genotype_ccm,
                            .y = genotype_eff,
                            q_value_threshold = 0.01))

genotype_paths_24 = genotype_models_tbl_24 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                           .x = genotype_ccm,
                           .y = genotype_eff,
                            q_value_threshold = 0.01))

genotype_paths_30 = genotype_models_tbl_30 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                            .x = genotype_ccm,
                            .y = genotype_eff,
                            q_value_threshold = 0.01))


genotype_paths_36 = genotype_models_tbl_36 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                            .x = genotype_ccm,
                            .y = genotype_eff,
                            q_value_threshold = 0.01))

genotype_paths_42 = genotype_models_tbl_42 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                            .x = genotype_ccm,
                            .y = genotype_eff,
                            q_value_threshold = 0.01))

genotype_paths_48 = genotype_models_tbl_48 %>%
  mutate(path = purrr::map2(.f = purrr::possibly(get_path, NA_real_),
                            .x = genotype_ccm,
                            .y = genotype_eff,
                            q_value_threshold = 0.01))



colData(kidney_cds) %>%
  as.data.frame() %>%
  filter(gene_target == "egfrb") %>% pull(timepoint) %>% unique()


all_path_df = rbind(genotype_paths_18,
                genotype_paths_24,
                genotype_paths_30,
                genotype_paths_36,
                genotype_paths_48) %>%
            filter(!is.na(path)) %>%
            tidyr::unnest(genotype_eff) %>%
            filter(delta_q_value < 0.01)


smo18 = genotype_paths_18 %>% slice(4)
noto18 = genotype_paths_18 %>% slice(8)

smo24 = genotype_paths_24 %>% slice(4)
noto24 = genotype_paths_24 %>% slice(8)


smo30 = genotype_paths_30 %>% slice(4)
noto30 = genotype_paths_30 %>% slice(8)


smo36 = genotype_paths_36 %>% slice(4)
noto36 = genotype_paths_36 %>% slice(8)

all_paths = rbind(genotype_paths_18 %>% mutate(timepoint = 18),
                  genotype_paths_24 %>% mutate(timepoint = 24),
                  genotype_paths_30 %>% mutate(timepoint = 30),
                  genotype_paths_36 %>% mutate(timepoint = 36),
                  genotype_paths_42 %>% mutate(timepoint = 42),
                  genotype_paths_48 %>% mutate(timepoint = 48)) %>%
              filter(!is.na(path)) %>%
              tidyr::unnest(path)

smo_path = all_paths %>%
  filter(gene_target == "smo") #%>%
  filter(timepoint == 24)


noto_path = all_paths %>%
  filter(gene_target == "noto") %>%
  filter(timepoint == 24)


plot_mt_edges = function(segment_plot,
                         path,
                         space_time_summary,
                         color="red") {

  mt_path_df = path %>%
    left_join(space_time_summary, by = c("from" = "cell_group")) %>%
    dplyr::rename("from_mean_hpf" = "mean_hpf", "from_mean_segment" = "mean_segment", "from_cell_type" = cell_type) %>%
    left_join(space_time_summary, by = c("to" = "cell_group")) %>%
    dplyr::rename("to_mean_hpf" = "mean_hpf", "to_mean_segment" = "mean_segment", "to_cell_type" = cell_type)


  label_df = space_time_summary %>%
              filter(cell_group %in% union(mt_path_df$from, mt_path_df$to))

  hp = segment_plot +
    geom_segment(data = mt_path_df ,
                 aes(x = from_mean_hpf,
                     y = from_mean_segment,
                     xend=to_mean_hpf,
                     yend = to_mean_segment),
                 color=color) +
    geom_segment(data = mt_path_df ,
                 aes(x = from_mean_hpf,
                     y = from_mean_segment,
                     xend=(from_mean_hpf+to_mean_hpf)/2,
                     yend = (from_mean_segment+to_mean_segment)/2),
                 color=color,
                 linejoin='mitre',
                 arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) #+
    ggrepel::geom_label_repel(data = label_df,
                              mapping = aes(mean_hpf, mean_segment, label=paste(cell_group," - ", cell_type)),
                              fill = "red")
  return(hp)
}

plot_mt_edges(segment_plot_all_labels, smo_path, space_time_summary)

plot_mt_edges(segment_plot_all_labels, noto_path, space_time_summary)



# seg_levels = space_time_summary %>% arrange(cell_group) %>%
#   ungroup() %>% mutate(new_cell_group = row_number()) %>%
#   tibble::column_to_rownames("new_cell_group")
#
# seg_levels = seg_levels %>%
#   arrange(mean_segment) %>%
#   mutate(rn = row_number()) %>%
#   mutate(group = cut(rn, 5, labels=F))
#
# G = igraph::graph_from_data_frame(data.frame(results_wl$paths))
#
# lay1 = get_layout(G, seg_levels)
#
# g = ggnetwork::ggnetwork(G,
#                          layout = lay1$layout,
#                          arrow.gap = 0.02, multiple = T, loops=T)
#
# # use the x coords from this
# x_coords = g %>% select(x, name) %>% distinct()
#
# space_time_summary_x = space_time_summary %>% left_join(x_coords, by = c("cell_group" = "name"))
#
#
# direct_edge_x = results_wl$paths %>%
#   left_join(space_time_summary_x, by = c("from" = "cell_group")) %>%
#   dplyr::rename("from_mean_hpf" = "mean_hpf", "from_mean_segment" = "mean_segment",
#                 "from_cell_type" = cell_type, "from_x" = x) %>%
#   left_join(space_time_summary_x, by = c("to" = "cell_group")) %>%
#   dplyr::rename("to_mean_hpf" = "mean_hpf", "to_mean_segment" = "mean_segment",
#                 "to_cell_type" = cell_type, "to_x" = x)
#
#
#
# ggplot() +
#   geom_segment(data = direct_edge_x ,
#                aes(x = from_x,
#                    y = from_mean_segment,
#                    xend=to_mean_hpf,
#                    yend = to_mean_segment),
#                color="black") +
#   geom_segment(data = direct_edge_x ,
#                aes(x = from_x,
#                    y = from_mean_segment,
#                    xend=(from_x+to_x)/2,
#                    yend = (from_mean_segment+to_mean_segment)/2),
#                color="black",
#                linejoin='mitre',
#                arrow = arrow(type="closed", angle=30, length=unit(1, "mm"))) +
#   scale_size_identity() #+
#   ggrepel::geom_label_repel(data = label_df,
#                             mapping = aes(mean_hpf, mean_segment, label=paste(cell_group," - ", cell_type)),
#                             # size=I(group_label_size),
#                             fill = "white")
#
# # ggnetwork::ggnetwork(G) %>%
#   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#   ggnetwork::geom_edges()+
#   ggnetwork::geom_nodes( size = 7,colour="black",shape=21)


hp18 = plot_mt_edges(segment_plot, smo18$path[[1]], space_time_summary)
hp18_2 = plot_mt_edges(segment_plot, noto18$path[[1]], space_time_summary, color = "blue")
hp18_2

hp24 = plot_mt_edges(segment_plot, smo24$path[[1]], space_time_summary)
hp24_2 = plot_mt_edges(hp24, noto24$path[[1]], space_time_summary, color = "blue")
hp24_2

hp30 = plot_mt_edges(segment_plot, smo30$path[[1]], space_time_summary)
hp30_2 = plot_mt_edges(hp30, noto30$path[[1]], space_time_summary, color = "blue")
hp30_2

hp36 = plot_mt_edges(segment_plot, smo36$path[[1]], space_time_summary)
hp36_2 = plot_mt_edges(hp36, noto36$path[[1]], space_time_summary, color = "blue")
hp36_2




# where is the mainly smo/noto cluster

plot_cells_3d(kidney_cds, color_cells_by = "cluster")
plot_cells_3d(kidney_cds, color_cells_by = "aberant.noto.smo.cells")


colData(kidney_cds) %>%
  as.data.frame() %>%
  filter(aberant.noto.smo.cells) %>%
  group_by(cluster) %>% tally() %>% arrange(-n)

##33





# # what are the negative relationships -----------------------------------------
#
#   # computed segment id for each clsuter, na.rm =
#
#   wt_ccm@best_full_model %>% return_igraph()
#
#
# neg_pcor = wt_ccm@best_full_model %>%
#   return_igraph() %>%
#   igraph::as_data_frame()
# filter(weight < 0)
#
# # where are these contrasts changing
# # find max abundance at each timepoint
# # for a given cell group, where are they changing
# max_abund = timepoint_pred_df %>%
#   group_by(cell_group) %>%
#   top_n(n = 1, wt = max(log_abund))
#
# min_abund = timepoint_pred_df %>%
#   group_by(cell_group) %>%
#   top_n(n = 1, wt = max(log_abund)) %>% ungroup
#
# min_abund %>% as.data.frame() %>%
#   group_by(cell_group) %>%
#   top_n(n=1, wt=min(timepoint)) %>%
#   filter(cell_group == 1)
#
# pos_pcor = wt_ccm@best_full_model %>%
#   return_igraph() %>%
#   igraph::as_data_frame() %>%
#   filter(weight > 0)
#
#
# # plot contrasts ---------------------------------------------------------------
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 18),
#                                          timepoint_pred_df %>% filter(timepoint == 24)),
#               q_value_thresh = 0.05)
#
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 20),
#                                          timepoint_pred_df %>% filter(timepoint == 26)),
#               q_value_thresh = 0.01)
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 24),
#                                          timepoint_pred_df %>% filter(timepoint == 30)),
#               q_value_thresh = 0.01)
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 30),
#                                          timepoint_pred_df %>% filter(timepoint == 36)),
#               q_value_thresh = 0.01)
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 36),
#                                          timepoint_pred_df %>% filter(timepoint == 42)),
#               q_value_thresh = 0.01)
#
# plot_contrast(wt_ccm, compare_abundances(wt_ccm,
#                                          timepoint_pred_df %>% filter(timepoint == 42),
#                                          timepoint_pred_df %>% filter(timepoint == 48)),
#               q_value_thresh = 0.01)
#
#
# comp_24_30 = compare_abundances(wt_ccm,
#                                 timepoint_pred_df %>% filter(timepoint == 24),
#                                 timepoint_pred_df %>% filter(timepoint == 30))
#
#
# # let's start with a mature cell type at 48 hrs and keep predicting back ------
#
# comp_30_36 = compare_abundances(wt_ccm,
#                                 timepoint_pred_df %>% filter(timepoint == 30),
#                                 timepoint_pred_df %>% filter(timepoint == 36))
#
# q_value_threshold = 0.05
# ccm = wt_ccm
# cond_b_vs_a_tbl = comp_30_36
#
# pos_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
#   as_tibble %>%
#   filter(pcor > 0 &
#            to_delta_q_value < q_value_threshold &
#            from_delta_q_value < q_value_threshold)
#
# weighted_edges = get_weighted_edges(ccm, pos_edges)
#
# neg_rec_edges = hooke:::collect_pln_graph_edges(ccm, cond_b_vs_a_tbl) %>%
#   as_tibble %>%
#   filter(edge_type != "undirected" &
#            to_delta_q_value < q_value_threshold &
#            from_delta_q_value < q_value_threshold)
#
# edge_path = neg_rec_edges %>%
#   dplyr::mutate(shortest_path = purrr::map2(.f =
#                                               purrr::possibly(get_shortest_path, NA_real_),
#                                             .x = from, .y = to,
#                                             weighted_edges)) %>%
#   select(shortest_path) %>%
#   tidyr::unnest(shortest_path) #%>%
# # select(-weight) %>%
# # distinct()
#
#
# # if there are no shortest paths, just use the neg reciprocal path
# neg_rec_edges %>% select(from, to, pcor, to_timepoint_x, to_timepoint_y)
#
# get_neg_dir_edges(wt_ccm, comp_30_36, q_value_threshold=0.00) %>% select(from, to, pcor, to_timepoint_x, to_timepoint_y)
#
#
# ### aggregate across timepoint windows ----------------------------------------
#
# # currently let's just get the neg reciprocal edges across all timepoint contrasts
#
#
