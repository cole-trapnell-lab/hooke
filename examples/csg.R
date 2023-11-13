library(monocle3)
library(hooke)

cds = readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/all-geno_sensory-cranial-ganglion_neuron_29k_cds.RDS")

ganglia_colors =
  c("cranial ganglion progenitor" = "#A29ADE",
    "trigeminal ganglion" = "#B4286C",
    "epibranchial ganglion" = "#2E74BC",
    "lateral line ganglion" = "#5ABC6B",
    "statoacoustic ganglion" = "#5ED3F3",
    "unknown sensory ganglion" = "#ADADAD",
    "rohon-beard neuron" = "#C9AE56")

plot_cells(cds, x = 1, y = 3, color_cells_by = "cell_type_sub",
           label_groups_by_cluster = F,  show_trajectory_graph = F) + ggtitle("test")+
  scale_color_manual(values = ganglia_colors) + theme(plot.title  = element_text(size=20))


control_ids = c("ctrl-uninj", "ctrl-inj", "ctrl-hgfa", "ctrl-met", "ctrl-mafba", "ctrl-noto", "ctrl-tbx16")
cds = cds[, colData(cds)$gene_target %in% c("phox2a", "foxi1", control_ids)]

colData(cds)$perturbation = ifelse(colData(cds)$gene_target %in% control_ids,
                                   "control",
                                   colData(cds)$gene_target)

ccs = new_cell_count_set(cds,
                         sample_group = "embryo",
                         cell_group = "cell_type_sub")

start_time = 18
stop_time = 72
time_formula = build_interval_formula(ccs, num_breaks = 3, interval_start = 18, interval_stop = 72)

ccm = new_cell_count_model(ccs,
                           main_model_formula_str = paste0("perturbation +",  time_formula),
                           nuissance_model_formula_str = "~ expt")

# predict for 48 hpf
cond_wt = estimate_abundances(ccm, tibble(timepoint = 48, perturbation = "control"))
cond_phox2a = estimate_abundances(ccm, tibble(timepoint = 48, perturbation = "phox2a"))
cond_foxi1 = estimate_abundances(ccm, tibble(timepoint = 48, perturbation = "foxi1"))

wt_v_phox2a_tbl = compare_abundances(ccm, cond_wt, cond_phox2a)
wt_v_foxi1_tbl = compare_abundances(ccm, cond_wt, cond_foxi1)


plot_contrast(ccm, wt_v_phox2a_tbl, x=1, y=3, q_value_threshold = 0.05)
plot_contrast(ccm, wt_v_foxi1_tbl, x=1, y=3, q_value_threshold = 0.05)



# wt kinetics

wt_cds = cds[, colData(cds)$gene_target %in% control_ids]
wt_ccs = new_cell_count_set(wt_cds,
                            sample_group = "embryo",
                            cell_group = "cell_type_sub")

colData(wt_ccs)$timepoint = as.numeric(colData(wt_ccs)$timepoint)
wt_ccm = new_cell_count_model(wt_ccs,
                              main_model_formula_str = "ns(timepoint, df=3)")

wt_timepoint_pred_df = estimate_abundances_over_interval(wt_ccm,
                                                         interval_start=18,
                                                         interval_stop=72,
                                                         interval_col="timepoint",
                                                         interval_step=2)

log_abund_detection_thresh=-3

ggplot(wt_timepoint_pred_df, aes(x = timepoint)) +
  geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh))) +
  facet_wrap(~cell_group, scales="free_y") + monocle3:::monocle_theme_opts()

# wt kinetics w batch

wt_expt_ccm = new_cell_count_model(wt_ccs,
                                   main_model_formula_str = "ns(timepoint, df=3)",
                                   nuisance_model_formula_str = "~ expt")

batches = batches %>% mutate(tp_preds = purrr::map(.f = function(batch) {
  estimate_abundances_over_interval(wt_expt_ccm,
                                    start_time,
                                    stop_time,
                                    knockout=FALSE,
                                    interval_col="timepoint",
                                    interval_step=2,
                                    expt = batch)
}, .x=batch))

wt_timepoint_pred_df = batches %>% select(tp_preds) %>% tidyr::unnest(tp_preds)

ggplot(wt_timepoint_pred_df, aes(x = timepoint)) +
  geom_line(aes(y = exp(log_abund) + exp(log_abund_detection_thresh), color=expt)) +
  facet_wrap(~cell_group, scales="free_y", nrow = 2) + monocle3:::monocle_theme_opts() +
  ggtitle("wild-type kinetics by expt")


# perturbation kinetics

foxi1_cds = cds[, colData(cds)$gene_target %in% c("foxi1", control_ids)]
foxi1_ccs = new_cell_count_set(foxi1_cds,
                               sample_group = "embryo",
                               cell_group = "cell_type_sub")

colData(foxi1_ccs)$timepoint = as.numeric(colData(foxi1_ccs)$timepoint)


start_time = 18
stop_time = 72

time_formula = build_interval_formula(foxi1_ccs, num_breaks = 3, interval_start = 18, interval_stop = 72)

foxi1_ccm = new_cell_count_model(foxi1_ccs,
                                 main_model_formula_str = past0("perturbation" + time_formula))

wt_timepoint_pred_df = estimate_abundances_over_interval(foxi1_ccm,
                                                         interval_start=start_time,
                                                         interval_stop=stop_time,
                                                         interval_col="timepoint",
                                                         interval_step=2,
                                                         perturbation = "control")

ko_timepoint_pred_df = estimate_abundances_over_interval(foxi1_ccm,
                                                         interval_start=start_time,
                                                         interval_stop=stop_time,
                                                         interval_col="timepoint",
                                                         interval_step=2,
                                                         perturbation = "foxi1")
timepoints = seq(start_time, stop_time, 2)
perturb_vs_wt_nodes = tibble(t1=timepoints) %>%
  mutate(comp_abund = purrr::map(.f = compare_ko_to_wt_at_timepoint,
                                 .x = t1,
                                 perturbation_ccm=foxi1_ccm,
                                 interval_col="timepoint",
                                 wt_pred_df = wt_timepoint_pred_df,
                                 ko_pred_df = ko_timepoint_pred_df)) %>%
  tidyr::unnest(comp_abund)


ggplot(perturb_vs_wt_nodes, aes(x = t1)) +
  geom_line(aes(y = exp(log_abund_x) + exp(log_abund_detection_thresh), linetype = "Wild-type")) +
  geom_line(aes(y = exp(log_abund_y) + exp(log_abund_detection_thresh), linetype = "Knockout")) +
  ggh4x::stat_difference(aes(ymin = exp(log_abund_x)+exp(log_abund_detection_thresh), ymax = exp(log_abund_y) +exp(log_abund_detection_thresh)), alpha=0.3) +
  facet_wrap(~cell_group, scales="free_y", nrow=2) + monocle3:::monocle_theme_opts() + ggtitle("foxi1 kinetics")


