#' gets counts
#' @param ccs
#' @return a cell group by embryo count matrix
get_cell_wide <- function(ccs) {
  # cell group x embryo
  cell_group_wide = counts(ccs) %>%
    as.matrix %>%
    as.data.frame

  return(cell_group_wide)
}

#' takes in a ccs and outputs a matrix of
#' cell_type x embryo of proportions
#' @param cell_group_wide a cell group x embryo count matrix
#' @param pseudocount
get_prop_mat <- function(cell_count_wide,
                         pseudocount = 1) {

  prop_mat = cell_count_wide %>%
    rownames_to_column("cell_group") %>%
    pivot_longer(-cell_group, names_to = "embryo", values_to = "count") %>%
    mutate(count = count + pseudocount) %>%
    group_by(embryo) %>%
    mutate(prop = (count)/sum(count)) %>%
    select(-count) %>%
    pivot_wider(names_from = cell_group, values_from = prop) %>%
    column_to_rownames("embryo") %>%
    as.matrix

  return(prop_mat)

}

#' fits a drichlet model based off proportions
#' @param prop_mat
#' @param model_formula_str
#' @param trace
fit_drichlet <- function(prop_mat,
                         model_formula_str = "~ 1",
                         trace = FALSE) {


  model_formula_str = paste("prop_mat", model_formula_str)
  model_formula = as.formula(model_formula_str)

  dfit <- VGAM::vglm(model_formula,
                     data = as.data.frame(prop_mat),
                     family = "dirichlet",
                     trace = trace)

  return(dfit)
}

simulate_seq <- function(cell_count_wide,
                         dfit,
                         genotype,
                         embryo_dim = NULL,
                         random_seed = 111,
                         sim_range = 10,
                         size = 1000) {

  sims = seq(1, sim_range, 1)

  if (is.null(embryo_dim)) {
    embryo_dim = ncol(cell_count_wide)
  }

  all_dat = lapply(sims, function(x){
    dsim <- VGAM::simulate.vlm(dfit, nsim = x)
    aaa <- array(unlist(dsim), c(embryo_dim*x, ncol(fitted(dfit)), x))
    aaa = aaa[,,1]
  })

  sim_props = reshape2::melt(all_dat)

  colnames(sim_props) = c("embryo", "cell_type", "cell_prop", "sim")
  sim_props$embryo_unq = paste(sim_props$sim, sim_props$embryo, sep='-')

  cell_type_names = rownames(cell_count_wide)

  # sim_props = left_join(sim_props %>% mutate(cell_type = as.character(cell_type)),
  #                       data.frame(cell_type_names) %>% rownames_to_column("cell_type"), by = "cell_type") %>%
  #             select(-cell_type)%>% rename("cell_type" = cell_type_names)

  cell_totals = colSums(cell_count_wide)
  my.sd = sd(cell_totals)
  sim.total.cells = round(rnorm(length(unique(sim_props$embryo_unq)), mean = size, sd = my.sd))
  tot_df = data.frame(unique(sim_props$embryo_unq), sim.total.cells, stringsAsFactors = F)
  colnames(tot_df) = c('embryo_unq', 'total_cells')
  sim_count_df = sim_props %>%
    left_join(tot_df, by='embryo_unq') %>%
    mutate(cells = round(total_cells*cell_prop))

  sim_filt = sim_count_df %>%
    left_join(sim_count_df %>%
                select(sim, embryo_unq) %>%
                distinct() %>%
                group_by(sim) %>%
                dplyr::mutate(sim_emb_count = dplyr::row_number()) %>%
                ungroup() %>%
                select(embryo_unq, sim_emb_count), by = "embryo_unq") %>%
    mutate(cell_type = paste("cell_type", cell_type, sep = "_")) %>%
    select(sim_emb_count, cell_type, sim, embryo_unq, cells, total_cells) %>%
    group_by(sim) %>%
    add_count() %>%
    mutate(sim_emb_number = n/max(sim_count_df$cell_type), genotype = genotype) %>%
    select(-n) %>%
    ungroup() %>%
    arrange(embryo_unq)

  return(sim_filt)

}

#' simulate a number of embryos
#' @param cell_group_wide
#' @param dfit
#' @param emb_num
#' @param genotype
simulate <- function(cell_group_wide,
                     dfit,
                     emb_num,
                     genotype,
                     cell_type = NULL,
                     effect_size = NULL,
                     random.seed = 111,
                     size = 1000,
                     num_sim = 10) {

  # simulate wt abundances
  # num_sim = 10
  # sims = seq(1, num_sim, 1)

  set.seed(random.seed)
  embryo_dim = ncol(cell_group_wide)

  # dat = lapply(sims, function(x){
  #   # this returns a (emb_num x cell_group) x nsim
  #   dsim <- VGAM::simulate.vlm(dfit, nsim = x)
  #   # reshape
  #   aaa <- array(unlist(dsim), c(emb_num*x, ncol(fitted(dfit)), x))
  #   aaa = aaa[,,1]
  # })

  # x = 10
  dsim <- VGAM::simulate.vlm(dfit, nsim = num_sim)
  aaa <- array(unlist(dsim), c(embryo_dim, ncol(fitted(dfit)), num_sim))
  aaa = aaa[,,1]
  aaa = aaa[1:emb_num,]

  # this returns a (num of embryo x 10 simulations) x 139

  sim_props = reshape2::melt(aaa)
  colnames(sim_props) = c("embryo", "cell_group_number", "cell_prop") #, "sim")

  cell_group_names = rownames(cell_group_wide)
  df = data.frame("cell_group_number" = 1:length(cell_group_names),
                  "cell_group" = cell_group_names)
  sim_props = sim_props %>% left_join(df, by = "cell_group_number")

  # sim_props = sim_props %>% group_by(embryo) %>% mutate(sim = row_number())

  # sim_props$embryo_unq = paste(sim_props$sim, sim_props$embryo, sep='-')

  # how many total cells does each embryo have
  cell_totals = colSums(cell_group_wide)
  # get standard deviation from embryo totals in our experiment
  my.sd = sd(cell_totals)

  # simulate total cells
  sim.total.cells = round(rnorm(length(unique(sim_props$embryo)), mean = size, sd = my.sd))


  tot_df = data.frame(unique(sim_props$embryo), sim.total.cells, stringsAsFactors = F)
  colnames(tot_df) = c('embryo', 'total_cells')

  # convert to counts
  sim_count_df = sim_props %>%
    # left_join(tot_df, by='embryo_unq') %>%
    left_join(tot_df, by='embryo') %>%
    mutate(cells = round(total_cells*cell_prop)) %>%
    select(-c(total_cells, cell_prop, cell_group_number))

  # mutate these if defined
  # if (is.null(cell_type) & is.null(effect_size)) {
  #
  # }


  sim_count_wide = sim_count_df %>%
    pivot_wider(names_from = "embryo", values_from = cells) %>%
    tibble::column_to_rownames("cell_group") %>%
    as.matrix


  # make a cell_count_cds
  covariates_df = data.frame("embryo" = colnames(sim_count_wide), genotype = genotype)

  cell_count_cds = monocle3::new_cell_data_set(sim_count_wide,
                                               cell_metadata=covariates_df)

  # hack this for now?
  sim_ccs = methods::new("cell_count_set",
                         cell_count_cds,
                         cds=cell_count_cds)

  return(sim_ccs)

  # sim_filt = sim_count_df %>%
  #   left_join(sim_count_df %>%
  #               select(sim, embryo_unq) %>%
  #               distinct() %>%
  #               group_by(sim) %>%
  #               dplyr::mutate(sim_emb_count = dplyr::row_number()) %>%
  #               ungroup() %>%
  #               # select(embryo_unq, sim_emb_count), by = "embryo_unq") %>%
  #               select(embryo_unq, sim_emb_count), by = "embryo_unq") %>%
  #   mutate(cell_type = paste("cell_type", cell_type, sep = "_")) %>%
  #   select(sim_emb_count, cell_type, sim, embryo_unq, cells, total_cells) %>%
  #   group_by(sim) %>%
  #   add_count() %>%
  #   mutate(sim_emb_number = n/max(sim_count_df$cell_type), genotype = genotype) %>%
  #   select(-n) %>%
  #   ungroup() %>%
  #   arrange(embryo_unq)

  # return(sim_filt)
}


#' mutate ccs
#'@param sim_ccs
#'@param cell_type
#'@param effect_size
mutate_ccs = function(sim_ccs,
                      cell_type,
                      effect_size=0) {

  sim_count_wide = get_cell_wide(sim_ccs)
  sim_count_df = sim_count_wide %>%
                 rownames_to_column("cell_group") %>%
                 pivot_longer(-cell_group,
                              names_to = "embryo",
                              values_to = "cells")

  sim_count_df = sim_count_df %>%
                 mutate(cells = ifelse(cell_group == cell_type,
                                       round(cells * (1-effect_size)),
                                       cells))

  sim_count_wide = sim_count_df %>%
                   pivot_wider(names_from = "embryo",
                               values_from = cells) %>%
                   tibble::column_to_rownames("cell_group") %>%
                   as.matrix

  covariates_df = colData(sim_ccs)
  covariates_df$effect = effect_size
  rownames(covariates_df) = seq(1:dim(sim_count_wide)[2])
  colnames(sim_count_wide) = seq(1:dim(sim_count_wide)[2])

  cell_count_cds = monocle3::new_cell_data_set(sim_count_wide,
                                               cell_metadata=covariates_df)

  # now turn it back into a cell count ccs
  mt_ccs = methods::new("cell_count_set",
                        cell_count_cds,
                        cds=cell_count_cds)

  return(mt_ccs)
}

#' plots boxplot of counts between two ccs
#' @param wt_ccs
#' @param mt_ccs
#' @param cell_type
#'
compare_counts <- function(wt_ccs,
                           mt_ccs,
                           cell_type) {

  wt_counts = data.frame("counts" = counts(wt_ccs)[cell_type,])
  mt_counts = data.frame("counts" = counts(mt_ccs)[cell_type,])

  wt_counts$type = "WT"
  mt_counts$type = "perturbed"

  p = rbind(wt_counts, mt_counts) %>%
    ggplot(aes(x = type, y = counts, fill = type)) +
    geom_boxplot() +
    monocle3:::monocle_theme_opts()

  return(p)

}

#'
#'@param sim_df
#'
sim_df_to_ccs <- function(sim_df) {

  df_pivot = sim_df %>% pivot_wider(names_from = cell_type, values_from = cells)
  mat_cols = colnames(df_pivot)[grepl("cell_type", colnames(df_pivot))]
  df_mat = df_pivot[,mat_cols] %>% as.matrix()
  df_cols = colnames(df_pivot)[!grepl("cell_type", colnames(df_pivot))]
  cov_df = df_pivot[,df_cols]

  cell_count_cds = monocle3::new_cell_data_set(t(df_mat),
                                               cell_metadata=cov_df)

  # cell_count_cds <- cds[,Matrix::colSums(exprs(cell_count_cds)) != 0]
  # cell_count_cds <- estimate_size_factors(cell_count_cds)

  # hack this for now?
  sim_ccs = methods::new("cell_count_set",
                         cell_count_cds,
                         cds=cell_count_cds)

  colnames(sim_ccs) = colData(sim_ccs)$embryo_unq


  return(sim_ccs)
}

#' combine the wt and mt ccs into one ccs
#' @param wt_ccs the ctrl cell count set
#' @param mt_ccs the perturbed cell count set
#'
combine_simulation <- function(wt_ccs,
                               mt_ccs) {


  colnames(wt_ccs) = paste0("WT", colnames(wt_ccs))
  colnames(mt_ccs) = paste0("MT", colnames(mt_ccs))

  comb_counts = cbind(counts(wt_ccs),
                      counts(mt_ccs))

  colData(mt_ccs)$genotype = "MT"

  if (("effect" %in% colnames(colData(wt_ccs))) == FALSE) {
    colData(wt_ccs)$effect = 0
  }

  comb_coldata = rbind(colData(wt_ccs), colData(mt_ccs))

  comb_coldata$embryo = paste0(comb_coldata$embryo_unq, comb_coldata$genotype )
  cell_count_cds = monocle3::new_cell_data_set(comb_counts,
                                               cell_metadata=comb_coldata)

  comb_ccs = methods::new("cell_count_set",
                          cell_count_cds,
                          cds=cell_count_cds)

  return(comb_ccs)

}

#' fit a ccm for each perturbation
#' @param wt_ccs
#' @param mt_ccs
#' @param cell_type
#' @param effect_size
#' @param model_formula_str
#'
fit_combine_ccm = function(wt_ccs,
                           mt_ccs,
                           ctrl_penalty_matrix = NULL,
                           model_formula_str = "~ genotype") {


  comb_ccs = combine_simulation(wt_ccs, mt_ccs)

  comb_ccm = new_cell_count_model(comb_ccs,
                                  main_model_formula_str = model_formula_str,
                                  penalty_matrix = ctrl_penalty_matrix)

  return(comb_ccm)

}

#'
#' @param comb_ccm
#' @param cell_type
#' @param q_value_threshold
#'
sim_comb_abund <- function(ccm,
                           cell_type,
                           q_value_threshold = 0.05) {

  # log abund sd is always the same? between both
  wt_abund = estimate_abundances(ccm, tibble(genotype = "WT"))
  mt_abund = estimate_abundances(ccm, tibble(genotype = "MT"))

  comp_abund = compare_abundances(ccm, wt_abund, mt_abund)

  comp_abund = comp_abund %>%
    mutate(delta_log_needed = qnorm(1-q_value_threshold, sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2))) %>%
    mutate(log_diff = abs(delta_log_abund) - delta_log_needed)

  return(comp_abund %>% filter(cell_group == cell_type))

}

#' @param which_sim
#' @param wt_sim
#' @param mt_sim
run_simulation = function(which_sim,
                          wt_sim,
                          mut_sim,
                          effects = c(0.01, 0.1, 0.25, 0.5, 0.75)) {

  wt_df = wt_sim %>% dplyr::filter(sim == which_sim)
  mut_df = mut_sim %>% dplyr::filter(sim == which_sim)

  wt_sim_ccs = sim_df_to_ccs(wt_df)
  mut_sim_ccs = sim_df_to_ccs(mut_df)

  df = expand.grid("cell_type" = unique(wt_df$cell_type),
                   "effect_size" = effects) %>%
    mutate(mutated_sim_ccs = purrr::map2(.f = purrr::possibly(mutate_ccs, NA_character_),
                                         .x = cell_type,
                                         .y = effect_size,
                                         sim_ccs = mut_sim_ccs)) %>%
    mutate(comb_ccm_sim = purrr::map(.f = purrr::possibly(fit_combine_ccm, NA_character_),
                                     .x = mutated_sim_ccs,
                                     wt_ccs = wt_sim_ccs,
                                     ctrl_penalty_matrix = ctrl_penalty_matrix)) %>%
    mutate(comb_abund_sim = purrr::map2(.f = purrr::possibly(sim_comb_abund, NA_character_),
                                        .x = comb_ccm_sim,
                                        .y = cell_type))
  df_unnest = df %>%
    as.data.frame %>%
    select(-comb_ccm_sim) %>%
    tidyr::unnest(c(comb_abund_sim))

  return(df_unnest)

}

