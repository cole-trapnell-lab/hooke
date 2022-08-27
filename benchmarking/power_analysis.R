library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(monocle3)
library(googledrive)
library(splines)
library(VGAM)
library(data.table)
library(hooke)
library(purrr)

setwd("~/OneDrive/UW/Trapnell/hooke_split_model/benchmarking/")
# setwd("/net/trapnell/vol1/home/duran/power_analysis/")
source("./power_analysis_utils.R")

# take a cds

cds <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/ctrl-uninj_24hpf_53k.rds")

# 1. make a ccs
ctrl_ccs = new_cell_count_set(cds,
                              sample_group = "embryo",
                              cell_group = "cell_type_broad")

# 2. get cell type proportions
cell_count_wide = get_cell_wide(ctrl_ccs) # cell_type x embryo
prop_mat = get_prop_mat(cell_count_wide)

# would be nice if insteal of getting cell types
ctrl_penalty_matrix = hooke:::init_penalty_matrix(ctrl_ccs)
rownames(ctrl_penalty_matrix) = paste0("cell_type_",seq(1:dim(cell_count_wide)[1]))
colnames(ctrl_penalty_matrix) = paste0("cell_type_",seq(1:dim(cell_count_wide)[1]))

# 3. fit a drichlet model to the cell type proportions
dfit = fit_drichlet(prop_mat)

cell_type_df = data.frame("cell_type" = colnames(prop_mat))
cell_type_df$rn = paste0("cell_type_", row_number(cell_type_df))
# cell_type_df$rn = row_number(cell_type_df)

# 4. simulate cell type abundances for our fake embryo total count data
sims = seq(1, 10, 1)
emb_num = nrow(prop_mat)

wt_dat = lapply(sims, function(x) {
  dsim <- VGAM::simulate.vlm(dfit, nsim = x)
  aaa <- array(unlist(dsim), c(emb_num*x, ncol(fitted(dfit)), x))
  aaa = aaa[,,1]
})

wt_sim_props = reshape2::melt(wt_dat)
colnames(wt_sim_props) = c("embryo", "cell_type", "cell_prop", "sim")
wt_sim_props$embryo_unq = paste(wt_sim_props$sim, wt_sim_props$embryo, sep='-')

saveRDS(wt_sim_props, "wt_sim_props.rds")

# # plot proportions
# wt_sim_props %>%
#   group_by(cell_type) %>%
#   mutate(ct_means = mean(cell_prop)) %>%
#   distinct(cell_type, ct_means) %>%
#   ggplot(aes(ct_means)) +
#   geom_histogram(bins = 20) +
#   theme_minimal() +
#   ylim("ct proportions")


# edit embryo size

embryo_size = 1000
wt_sim_filt = simulate_seq(cell_count_wide, dfit, genotype = "WT", random_seed=111, size = embryo_size)
mut_sim_filt = simulate_seq(cell_count_wide, dfit, genotype = "MT", random_seed=222, size = embryo_size)


# Test multiple effect sizes with increasing numbers of embryos

emb_interval = 5
embs = sims*emb_interval

v1 = unlist(lapply(embs, function(x) seq(1, x, 1)))
v2 = unlist(lapply(sims, function(x) rep(x, x*emb_interval)))

emb_set = tibble(v2, v1) %>%
  mutate(emb_name = paste(v2, v1, sep = "-")) %>%
  pull(emb_name)

wt_sim = wt_sim_filt %>%
  filter(embryo_unq %in% emb_set)

mut_sim = mut_sim_filt %>%
  filter(embryo_unq %in% emb_set)

cell_types = rownames(cell_count_wide)




res_df = lapply(seq(1, 10), function(which_sim){

  sim_df = run_simulation(which_sim,
                          wt_sim,
                          mut_sim)

  sim_df %>% fwrite(file = paste0(outdir,
                                  "simulation_results/simulation_which_sim_", which_sim ,
                                  "-emb_interval_", emb_interval,
                                  "-embryo_size_", size,
                                  "-num_emb_", num_emb,
                                  ".csv"), sep = ",")

})


# read in results -------------------------------------------------------------


result_files = list.files(path = "~/OneDrive/UW/Trapnell/hooke/examples/simulation_results", pattern = "which_sim", full.names = T, recursive = T)
df = lapply(result_files, function(f){

  fread(f, sep = ",", stringsAsFactors = F, data.table = F, na.strings = "") %>%
    mutate(which_sim = gsub(".csv", "",
                            gsub("/Users/maddyduran/OneDrive/UW/Trapnell/hooke/examples/simulation_results/simulation_which_sim_", "",
                                 f)))

}) %>% bind_rows()
df2 = fread("~/OneDrive/UW/Trapnell/hooke/examples/simulation_results/simulation_effect_1.csv", sep = ",", stringsAsFactors = F, data.table = F, na.strings = "")

df = rbind(df, df2)

hooke_ref <- function(df, wt_sim_props) {

  # add associated mean cell type proportions
  props = wt_sim_props %>%
    mutate(cell_type = paste0("cell_type_", cell_type)) %>%
    group_by(cell_type) %>%
    mutate(ct_means = mean(cell_prop)) %>%
    distinct(cell_type, ct_means)

  df = df %>%
    left_join(props, by = "cell_type") %>%
    mutate(cell_type_prop = paste("prop_", round(ct_means, 3), sep = ""))

  # order cell types by proportion and reset levels for plotting
  ct_order = df %>%
    select(cell_type_prop, ct_means) %>%
    arrange(-ct_means) %>%
    distinct(cell_type_prop) %>%
    pull(cell_type_prop)

  df$cell_type_prop = factor(df$cell_type_prop, levels = ct_order)
  df = df %>%
    mutate(is_sig = case_when(delta_p_value < 0.05 ~ TRUE,
                              delta_p_value > 0.05 ~ FALSE))

  return(df)
}

hooke_res = hooke_ref(df, wt_sim_props) %>%
  left_join(data.frame("which_sim" = as.character(sims), "emb_num" = embs), by = "which_sim") %>%
  group_by(emb_num, effect) %>%
  summarize(pct_pass = 100*sum(delta_p_value < 0.05)/n())

hooke_res%>%
  ggplot(aes(emb_num, pct_pass)) + geom_point() + facet_wrap(~effect) +
  monocle3:::monocle_theme_opts() + ylim(0,100)


ggplot(hooke_res, aes(x = emb_num, y = sig_percent)) +
  geom_point() +
  facet_wrap(~effect_size) +
  theme_minimal() +
  ggtitle("Percent of cell types passing significance across effect sizes")



# BB results ------------------------------------------------------------------

setwd("~/OneDrive/UW/Trapnell/hooke_split_model/benchmarking/")
drive_download("https://drive.google.com/file/d/10_9Pe5ppSSrmZuNuRXi60fbh-vns83kb", overwrite = TRUE)
bb.final.res.bind = fread("cell-abund_sim_res_5emb-groups.csv", sep = ",", data.table = F, stringsAsFactors = F)

wt_sim_props = readRDS("wt_sim_props.rds")

bb_ref <- function(final.res.bind, wt_sim_props) {

  # add associated mean cell type proportions
  props = wt_sim_props %>%
    mutate(cell_type = paste0("cell_type_", cell_type)) %>%
    group_by(cell_type) %>%
    mutate(ct_means = mean(cell_prop)) %>%
    distinct(cell_type, ct_means)

  final.res.bind = final.res.bind %>%
    left_join(props, by = "cell_type") %>%
    mutate(cell_type_prop = paste("prop_", round(ct_means, 3), sep = ""))

  # order cell types by proportion and reset levels for plotting
  ct_order = final.res.bind %>%
    select(cell_type_prop, ct_means) %>%
    arrange(-ct_means) %>%
    distinct(cell_type_prop) %>%
    pull(cell_type_prop)

  final.res.bind$cell_type_prop = factor(final.res.bind$cell_type_prop, levels = ct_order)
  final.res.bind = final.res.bind %>%
    mutate(is_sig = case_when(p.value < 0.05 ~ TRUE,
                              p.value > 0.05 ~ FALSE))

  return(final.res.bind)
}


bb_res_summary = bb_ref(final.res.bind, wt_sim_props) %>%
  group_by(effect_size, emb_num) %>%
  summarize(sig_percent = sum(is_sig)/sum(n()) * 100)

ggplot(bb_res_summary, aes(x = emb_num, y = sig_percent)) +
  geom_point() +
  facet_wrap(~effect_size) +
  theme_minimal() +
  ggtitle("Percent of cell types passing significance across effect sizes")

