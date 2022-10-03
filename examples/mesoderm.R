library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
#library(hooke)
library(garnett)
library(msigdbr)


setwd("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/mesoderm")

meso_cds = readRDS("final-ref_mesoderm-PA_350k_saved-model-algn_anno_cds.RDS")

#meso_cds = cluster_cells(meso_cds, resolution=1e-4, random_seed=42)

colData(meso_cds)$cluster = as.character(clusters(meso_cds))
meso_cds = meso_cds[, is.na(colData(meso_cds)$embryo) == FALSE]
colData(meso_cds)$cell_type = colData(meso_cds)$cell_type_sub


ggplot(aes(x=expt, y=mean_nn_time),
       data=colData(meso_cds) %>% as.data.frame) +
  geom_boxplot() + facet_wrap(~timepoint, scale="free") +
  monocle3:::monocle_theme_opts()

ggplot(aes(x=as.numeric(timepoint), y=mean_nn_time, color=temp),
       data=colData(meso_cds) %>% as.data.frame) + geom_jitter() +
  facet_wrap(~expt) + geom_abline(color="grey") + geom_smooth(method="lm", color="black") +
  ylim(c(16,100)) + monocle3:::monocle_theme_opts()

staging_df = colData(meso_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))

staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
staging_df$predicted_timepoint = predict(staging_model)

colData(meso_cds)$adjusted_timepoint = staging_df$predicted_timepoint

# JUST FOR DEBUGGING:
colData(meso_cds)$timepoint = NULL

wt_ccs = new_cell_count_set(meso_cds,
                            sample_group = "embryo",
                            cell_group = "cluster")

#
# wt_start = 18
# wt_stop = 96
# num_time_breaks = 5
# time_breakpoints = seq(wt_start, wt_stop, length.out=num_time_breaks)
# time_breakpoints = time_breakpoints[2:(length(time_breakpoints) - 1)] #exclude the first and last entry as these will become boundary knots
# wt_main_model_formula_str = paste("~ splines::ns(timepoint, knots=", paste("c(",paste(time_breakpoints, collapse=","), ")", sep=""), ")")


# Use custom breaks because the sampling over the time range is so uneven
#wt_main_model_formula_str = build_interval_formula(wt_ccs, interval_var="adjusted_timepoint", interval_start=18, interval_stop=96, num_breaks=6)
wt_main_model_formula_str = "~ splines::ns( adjusted_timepoint , knots= c(24,48,60), Boundary.knots=c(20,92) )"

wt_ccm_wl = new_cell_count_model(wt_ccs,
                                 #main_model_formula_str = "~expt",
                                 main_model_formula_str = wt_main_model_formula_str,
                                 #nuisance_model_formula_str = "~1",
                                 #nuisance_model_formula_str = "~expt",
                                 whitelist = initial_pcor_graph(wt_ccs))

wt_ccm_wl = select_model(wt_ccm_wl, criterion = "EBIC", sparsity_factor=1)

# xxx_interval_abundances = estimate_abundances_over_interval(wt_ccm_wl, 18, 96, interval_col="adjusted_timepoint")
# qplot(adjusted_timepoint, log_abund, color = log_abund - 3*log_abund_se > 0, geom="point", data=xxx_interval_abundances) +
#   #geom_ribbon(aes(ymin=log_abund - 2*log_abund_se, ymax=log_abund + 2*log_abund_se), alpha=0.25) +
#   facet_wrap(~cell_group)

plot_contrast_wrapper <- function(ccm, t1, t2, q_val=0.01) {

  timepoint_pred_df = estimate_abundances_over_interval(ccm, t1, t2, interval_col="adjusted_timepoint", expt="GAP16")

  plot_contrast(ccm, compare_abundances(ccm,
                                           timepoint_pred_df %>% filter(adjusted_timepoint == t1),
                                           timepoint_pred_df %>% filter(adjusted_timepoint == t2)),
                scale_shifts_by = "none",
                q_value_thresh = q_val)

}

#t1 = 18
#t2 = 22
#debug(plot_contrast_wrapper)
plot_contrast_wrapper(wt_ccm_wl, 18, 22)


# Reboot
# -----------------------------------------------------------------------------


state_transition_graph = assemble_timeseries_transitions(wt_ccm_wl,
                                                            start_time=18, stop_time=96,
                                                            interval_col="adjusted_timepoint",
                                                            min_interval = 2,
                                                            log_abund_detection_thresh=-2,
                                                            experiment="GAP14")


plot_state_transition_graph(wt_ccm_wl, state_transition_graph %>% igraph::as_data_frame(),
                            color_nodes_by = "cell_type_broad", group_nodes_by="cell_type_broad")


