suppressPackageStartupMessages({
  library(argparse)
  library(hooke)
  library(tidyverse)
  library(VGAM)
  library(data.table)
  library(pbapply)
  library(scuttle)
  # library(speckle)
  library(limma)
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

parser$add_argument("--method", type = "character", default = "hooke",
                    help = "which method to run")

parser$add_argument("--prop_umi", type = "double", default = NULL,
                    help = "downsample UMIs")

parser$add_argument("--prop_cell_type", type = "double", default = NULL,
                    help = "downsample cell types")

parser$add_argument("--cell_group", type = "character", default = "cell_type_broad",
                    help = "cell_type_sub or cell_type_broad")

parser$add_argument("--num_groups_perturbed", type = "integer", default = 1,
                    help = "how many cell types to perturb at once")

parser$add_argument("--num_bootstraps", type = "integer", default = 0,
                    help = "how times to bootstrap vhat")

parser$add_argument("--norm_method", type = "character", default = "size_factors")

parser$add_argument("--run_bb", action = "store_true")
parser$add_argument("--run_hooke",  action = "store_true")
parser$add_argument("--run_hooke_boot",  action = "store_true")
#propeller only runs on R 4.2
parser$add_argument("--run_propeller", action = "store_true")


args <- parser$parse_args()
seed <- args$batch_seed
file_path <- args$file_path
embryo_size <- args$embryo_size
method <- args$method
prop_umi <- args$prop_umi
prop_cell_type <- args$prop_cell_type
cell_group <- args$cell_group
num_groups_perturbed <- args$num_groups_perturbed
num_bootstraps <- args$num_bootstraps
run_propeller <- args$run_propeller
run_hooke <- args$run_hooke
run_bb <- args$run_bb
run_hooke_boot <- args$run_hooke_boot
norm_method <- args$norm_method

set.seed(seed)

print("Loading dataset...")

cds <- readRDS(file_path)


if (is.null(prop_umi) == FALSE) {
 cds = downsample_umis(cds, prop = prop_umi)
}

if (is.null(prop_cell_type) == FALSE) {
  cds = downsample_cell_types(cds, prop = prop_cell_type)
}


# 1. make a ccs
ctrl_ccs = new_cell_count_set(cds,
                              sample_group = "embryo",
                              cell_group = cell_group)

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

wt_sim_filt = simulate_seq(cell_count_wide, dfit, genotype = "WT", random_seed=seed, size = embryo_size)
mut_sim_filt = simulate_seq(cell_count_wide, dfit, genotype = "MT", random_seed=seed*2, size = embryo_size)

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

type_vec = unique(wt_sim$cell_type)

mut_sim = mut_sim_filt %>%
  filter(embryo_unq %in% emb_set)

cell_types = rownames(cell_count_wide)


if (method == "compare") {

  # num_emb = 10

  outdir = "/net/trapnell/vol1/home/duran/hooke_comparison/"
  # outdir = "~/OneDrive/UW/Trapnell/hooke_split_model/benchmarking/analysis/hooke_comparison"  # outdir = "~/OneDrive/UW/Trapnell/hooke_split_model/benchmarking/analysis/hooke_comparison"
  # dir.create(outdir)

  which_sim = 2
  effect_list = c(0.25)
  cell_types = c("cell_type_10", "cell_type_12", "cell_type_20", "cell_type_42")


  # run the hooke version

  if (run_hooke) {
    sim_hooke_df = run_simulation(which_sim,
                                  wt_sim,
                                  mut_sim,
                                  cell_types_to_perturb = cell_types,
                                  effects = effect_list,
                                  norm_method = norm_method,
                                  method = "hooke",
                                  random.seed = seed)

    sim_hooke_df %>%
      fwrite(file = paste0(outdir,
                           "simulation_which_sim_", which_sim,
                           "-embryo_size_", embryo_size,
                           "-method_", "hooke",
                           ifelse(is.null(seed), "",paste0("-seed_", seed)),
                           ifelse(is.null(cell_types), "all_cell_types", ""),
                           ifelse(norm_method == "size_factors", "", norm_method),
                           ".csv"), sep = ",")
  }

  if (run_hooke_boot) {
    sim_hooke_boot_df = run_simulation(which_sim,
                                       wt_sim,
                                       mut_sim,
                                       cell_types_to_perturb = cell_types,
                                       effects = effect_list,
                                       norm_method = norm_method,
                                       method = "hooke",
                                       num_bootstraps = num_bootstraps,
                                       random.seed = seed)
    sim_hooke_boot_df %>%
      fwrite(file = paste0(outdir,
                           "simulation_which_sim_", which_sim,
                           "-embryo_size_", embryo_size,
                           "-method_", "hooke",
                           ifelse(is.null(num_bootstraps), "",paste0("-num-bootstraps_", num_bootstraps)),
                           ifelse(is.null(seed), "",paste0("-seed_", seed)),
                           ifelse(is.null(cell_types), "all_cell_types", ""),
                           ifelse(norm_method == "size_factors", "", norm_method),
                           ".csv"), sep = ",")
  }

  if (run_propeller) {
    # # run the propeller version
    sim_prop_df = run_simulation(which_sim,
                                 wt_sim,
                                 mut_sim,
                                 cell_types_to_perturb = cell_types,
                                 effects = effect_list,
                                 method = "propeller",
                                 random.seed = seed)


    sim_prop_df %>%
      fwrite(file = paste0(outdir,
                           "simulation_which_sim_", which_sim,
                           "-embryo_size_", embryo_size,
                           "-method_", "propeller",
                           ifelse(is.null(seed), "",paste0("-seed_", seed)),
                           ifelse(is.null(cell_types), "all_cell_types", ""),
                           ".csv"), sep = ",")
  }

  if (run_bb) {
    # run bb
    sim_bb_df = run_simulation(which_sim,
                               wt_sim,
                               mut_sim,
                               cell_types_to_perturb = cell_types,
                               effects = effect_list,
                               method = "bb",
                               random.seed = seed)

    sim_bb_df %>%
      fwrite(file = paste0(outdir,
                           "simulation_which_sim_", which_sim,
                           "-embryo_size_", embryo_size,
                           "-method_", "bb",
                           ifelse(is.null(num_bootstraps), "",paste0("-num-bootstraps_", num_bootstraps)),
                           ifelse(is.null(seed), "", paste0("-seed_", seed)),
                           ifelse(is.null(cell_types), "all_cell_types", ""),
                           ".csv"), sep = ",")

  }

}

# try to perturb multiple cell types at once











