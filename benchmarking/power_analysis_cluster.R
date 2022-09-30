suppressPackageStartupMessages({
  library(argparse)
  library(hooke)
  library(tidyverse)
  library(VGAM)
  library(data.table)
})

setwd("/net/trapnell/vol1/home/duran/power_analysis/")

source("./power_analysis_utils.R")

parser <- ArgumentParser()

parser$add_argument("--file_path", type = "character",
                    help = "file path to WT data to simulate from")

parser$add_argument("--batch_seed", type = "integer", default = 42,
                    help = "")

parser$add_argument("--embryo_size", type = "integer", default = 1000,
                    help = "how many cells in an embryo")

# parser$add_argument("--emb_num", type = "integer", default = 8,
#                     help = "how many embryos for each condition")
#
# parser$add_argument("--num_celL_types", type = "integer", default = 85,
#                     help = "the number of cell types")
#
# parser$add_argument("--make_batch_effect", type = "character", default = "no",
#                     help = "should synthetic batch effects be added? (yes/no)")


args <- parser$parse_args()

seed <- args$batch_seed
file_path <- args$file_path
embryo_size <- args$embryo_size

set.seed(seed)

print("Loading dataset...")

cds <- readRDS(file_path)

# where to save data

outdir = "/net/trapnell/vol1/home/duran/hooke_simulations/"


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

  sim_df %>% select(-mutated_sim_ccs) %>%
    fwrite(file = paste0(outdir,
                                  "simulation_results/simulation_which_sim_", which_sim ,
                                  # "-emb_interval_", emb_interval,
                                  "-embryo_size_", embryo_size,
                                  # "-num_emb_", num_emb,
                                  ".csv"), sep = ",")

})







