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

  set.seed(random_seed)

  sims = seq(1, sim_range, 1)

  if (is.null(embryo_dim)) {
    embryo_dim = ncol(cell_count_wide)
  }

  all_dat = lapply(sims, function(x){
    dsim <- VGAM::simulate.vlm(dfit, nsim = x, seed = random_seed)
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
  # this was giving negative numbers for really low counts
  # sim.total.cells = round(rnorm(n = length(unique(sim_props$embryo_unq)), mean = size, sd = my.sd))
  sim.total.cells = rpois(length(unique(sim_props$embryo_unq)), size)
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
                      cell_types,
                      effect_size=0) {

  sim_count_wide = get_cell_wide(sim_ccs)
  sim_count_df = sim_count_wide %>%
                 rownames_to_column("cell_group") %>%
                 pivot_longer(-cell_group,
                              names_to = "embryo",
                              values_to = "cells")

  sim_count_df = sim_count_df %>%
                 mutate(cells = ifelse(cell_group %in% cell_types,
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
                        cds = cell_count_cds)

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
                           model_formula_str = "~ genotype",
                           nuisance_model_formula_str = NULL,
                           num_bootstraps = NULL) {


  comb_ccs = combine_simulation(wt_ccs, mt_ccs)

  comb_ccm = new_cell_count_model(comb_ccs,
                                  main_model_formula_str = model_formula_str,
                                  nuisance_model_formula_str = nuisance_model_formula_str,
                                  penalty_matrix = ctrl_penalty_matrix,
                                  num_bootstraps = num_bootstraps)

  return(comb_ccm)

}

#' runs compare abundance between the wt sim and mt sim
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

  # comp_abund = comp_abund %>%
    # mutate(delta_log_needed = qnorm(1-q_value_threshold, sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2))) %>%
    # mutate(log_diff = abs(delta_log_abund) - delta_log_needed)

  return(comp_abund)

}

#' @param num_groups_perturbed
#' @param effect_size
perturb_multiple_cell_types = function(num_groups_perturbed, effect_size, random.seed = 42) {

  set.seed(random.seed)

  cell_types = rownames(ccs)

  cell_types_to_perturb = sample(cell_types, num_groups_perturbed, replace = F)

  ccs = mutate_ccs(ccs, cell_types_to_perturb, effect_size = effect_size)

  return(ccs)
}

#' @param wt_ccs
#' @param mt_ccs
#' @param group_col
fit_propeller = function(wt_ccs,
                         mt_ccs,
                         group_col = "genotype"){

  comb_ccs = combine_simulation(wt_ccs, mt_ccs)

  sample_to_group = colData(comb_ccs) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    select(sample, !!sym(group_col))

  prop_coldata = counts(comb_ccs) %>%
    as.matrix %>%
    as.data.frame()  %>%
    tibble::rownames_to_column("cell_group") %>%
    pivot_longer(-cell_group) %>%
    uncount(value) %>%
    dplyr::rename(sample = name, cluster = cell_group) %>%
    left_join(sample_to_group, by = "sample") %>%
    dplyr::rename(group = genotype)

  prop_res = propeller(clusters = prop_coldata$cluster,
                       sample = prop_coldata$sample,
                       group = prop_coldata$group)

  return(prop_res)

}


fit_betabinomial <- function(wt_df,
                             mt_df) {


}

#' @param which_sim
#' @param wt_sim
#' @param mt_sim
#' @param all_cell_types
#' @param num_bootstraps
run_simulation = function(which_sim,
                          wt_sim,
                          mut_sim,
                          cell_types_to_perturb = NULL,
                          num_bootstraps = NULL,
                          method = c("hooke", "bb", "propeller"),
                          effects = c(0.01, 0.1, 0.25, 0.5, 0.75),
                          random.seed = 42) {

  set.seed(random.seed)

  wt_df = wt_sim %>% dplyr::filter(sim == which_sim)
  mut_df = mut_sim %>% dplyr::filter(sim == which_sim)

  wt_sim_ccs = sim_df_to_ccs(wt_df)
  mut_sim_ccs = sim_df_to_ccs(mut_df)

  if (is.null(cell_types_to_perturb)) {
    cell_types_to_perturb =  unique(wt_df$cell_type)
  }

  df = expand.grid("cell_type" = cell_types_to_perturb,
                   "effect_size" = effects)


  if (method == "hooke") {

    df = df %>%
      mutate(mutated_sim_ccs = purrr::map2(.f = purrr::possibly(mutate_ccs, NA_character_),
                                           .x = cell_type,
                                           .y = effect_size,
                                           sim_ccs = mut_sim_ccs))

    df = df %>%
      mutate(comb_ccm_sim = purrr::map(.f = fit_combine_ccm,
                                      # .f = purrr::possibly(fit_combine_ccm, NA_character_),
                                       .x = mutated_sim_ccs,
                                       wt_ccs = wt_sim_ccs,
                                       model_formula_str = "~ genotype",
                                       ctrl_penalty_matrix = ctrl_penalty_matrix,
                                       num_bootstraps = num_bootstraps))

    df = df %>% mutate(comb_abund_sim = purrr::map2(.f = purrr::possibly(sim_comb_abund, NA_character_),
                                                    .x = comb_ccm_sim,
                                                    .y = cell_type))
    df_unnest = df %>%
      as.data.frame %>%
      select(-c(comb_ccm_sim, mutated_sim_ccs)) %>%
      tidyr::unnest(c(comb_abund_sim)) %>%
      mutate(method = method)

  } else if (method == "bb") {


    df = df %>% mutate(fit_bb = purrr::map2(.f = purrr::possibly(beta_binom_test, NA_character_),
                                            .x = cell_type,
                                            .y = effect_size,
                                            wt_df = wt_df,
                                            mut_df = mut_df))

    df_unnest  = df %>%
      as.data.frame %>%
      select(-cell_type) %>%
      tidyr::unnest(c(fit_bb)) %>%
      mutate(method = method)

  } else if (method == "propeller") {

    df = df %>%
      mutate(mutated_sim_ccs = purrr::map2(.f = purrr::possibly(mutate_ccs, NA_character_),
                                           .x = cell_type,
                                           .y = effect_size,
                                           sim_ccs = mut_sim_ccs))

    df = df %>%
      mutate(prop_res = purrr::map(.f = purrr::possibly(fit_propeller, NA_character_),
                                   .x = mutated_sim_ccs,
                                   wt_ccs = wt_sim_ccs,
                                   group_col = "genotype"))

    df_unnest = df %>%
      as.data.frame() %>%
      select(-mutated_sim_ccs) %>%
      tidyr::unnest(c(prop_res)) %>%
      mutate(cell_group = BaselineProp.clusters) %>%
      mutate(method = method)

  } else if (method == "louvain") {

  }


  # where cell type is the cell that was mutated
  # cell group is the effect size of that cell group

  return(df_unnest)

}




# louvain + glm

# louvain_glm <- function() {
#
#   louvain.model <- model.matrix(~Condition, data=test.meta)
#   louvain.count <- as.matrix(table(sim.clust.merge$Louvain.Clust, sim.clust.merge$Sample))
#   louvain.dge <- DGEList(counts=louvain.count, lib.size=log(colSums(louvain.count)))
#   louvain.dge <- estimateDisp(louvain.dge, louvain.model)
#   louvain.fit <- glmQLFit(louvain.dge, louvain.model, robust=TRUE)
#   louvain.res <- as.data.frame(topTags(glmQLFTest(louvain.fit, coef=2), sort.by='none', n=Inf))
#   table(louvain.res$FDR <= 0.1)
#
# }


beta_binom_test <- function(wt_df, mut_df, which_type, effect_size) {

  # only edit the 1 cell type
  adj_df = mut_df %>% mutate(cells = ifelse(cell_type == which_type,
                                            round(cells * (1-effect_size)),
                                            cells))

  comb_df = rbind(wt_df, adj_df) %>%
    mutate(genotype = forcats::fct_relevel(genotype, c(unique(wt_df$genotype))))

  all_cell_types = unique(comb_df$cell_type)

  run_bb = function(df, x) {
    df = df %>% filter(cell_type == x)
    df$cells = as.integer(df$cells)
    df$total_cells = as.integer(df$total_cells)
    count_df = cbind(df$cells, df$total_cells - df$cells)
    fit = vglm(count_df ~ genotype, data = df, family = betabinomial, trace = F)
    fit_df = tidy.vglm(fit)[3,]
    fit_df
  }

  # cell type is the cell that is perturbed
  test_res =  data.frame("cell_type" = which_type,
                         "cell_group" = all_cell_types) %>%
    mutate(bb.fit = purrr::map(.f = purrr::possibly(run_bb, NA_character_),
                               .x = cell_group,
                               df = comb_df))


  test_res %>% filter(!is.na(bb.fit)) %>% tidyr::unnest(c(bb.fit))

}



# BB

binom_effect = function(wt_df, mut_df, which_type = 1, effect_vec = c(0.01, 0.1, 0.25, 0.5, 0.75)){
  test_res = lapply(effect_vec,
                    FUN = function(x) {

                      # filter for cell type
                      wt_df = wt_df %>% filter(cell_type == which_type)
                      head(wt_df)
                      mut_df = mut_df %>% filter(cell_type == which_type)
                      head(mut_df)

                      # adjust mutant data
                      adj_df = mut_df %>% mutate(cells = round(cells * (1-x)))

                      # combine and run test

                      comb_df = rbind(wt_df, adj_df) %>%
                        mutate(genotype = forcats::fct_relevel(genotype, c(unique(wt_df$genotype))))

                      comb_df$cells = as.integer(comb_df$cells)
                      comb_df$total_cells = as.integer(comb_df$total_cells)

                      count_df = cbind(comb_df$cells, comb_df$total_cells - comb_df$cells)
                      fit = vglm(count_df ~ genotype, data = comb_df, family = betabinomial, trace = F)
                      fit_df = tidy.vglm(fit)[3,]})
  test_res
}




binom_type = function(wt_df,
                      mut_df,
                      type_vec = type_vec,
                      effect_vec = effect_vec,
                      which_sim = 1) {

  wt_df = wt_df %>% dplyr::filter(sim == which_sim)
  mut_df = mut_df %>% dplyr::filter(sim == which_sim)
  celltype_res = lapply(type_vec, function(x){
    print(x)
    res = binom_effect(wt_df = wt_df, mut_df = mut_df, which_type = x, effect_vec = effect_vec)
    names(res) = effect_vec
    res.bind = do.call(rbind, res)
    res.bind$effect_size = row.names(res.bind)
    res.bind$cell_type = x
    res.bind
  })
  celltype_res
}

tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
  co <- as.data.frame(coef(summary(x)))
  names(co) <- c("estimate","std.error","statistic","p.value")
  if (conf.int) {
    qq <- qnorm((1+conf.level)/2)
    co <- transform(co,
                    conf.low=estimate-qq*std.error,
                    conf.high=estimate+qq*std.error)
  }
  co <- data.frame(term=rownames(co),co)
  rownames(co) <- NULL
  return(co)
}

compare_abundance = function(cell_df, wt_df){
  comb_df = rbind(as.data.frame(cell_df), as.data.frame(wt_df)) %>%
    mutate(genotype = fct_relevel(genotype, c(unique(wt_df$genotype))))

  cell.types = unique(comb_df$cell_type)

  test_res = sapply(cell.types,
                    FUN = function(x) {
                      type_df = comb_df %>% filter(cell_type == x)
                      count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
                      fit =  VGAM::vglm(count_df ~ genotype, data = type_df, family = "betabinomial", trace = F)
                      fit_df = tidy.vglm(fit)[3,]}, USE.NAMES = T, simplify = F)

  test_res = do.call(rbind, test_res)
  test_res = test_res %>% tibble::rownames_to_column(var = "cell_group")
  test_res %>% arrange(desc(estimate))
}


downsample_umis = function(cds, prop, bycol = FALSE) {
  counts(cds) = downsampleMatrix(counts(cds), prop = prop, bycol=bycol)
  cds = estimate_size_factors(cds)
  return(cds)
}

total_umi = function(cds) {
  colData(cds)$totalUMI <- Matrix::colSums(exprs(cds))
  return(cds)
}

downsample_cell_types = function(cds, prop, colname = "cell_type_broad", random.seed = 42) {
  set.seed(random.seed)
  cell_types = colData(cds)[[colname]] %>% unique()
  num_cell_types = cell_types %>% length()
  sampled_cell_types = sample(cell_types, size = round(num_cell_types*prop), replace = F)
  cds = cds[,colData(cds)[[colname]] %in% sampled_cell_types]
  return(cds)
}

downsample_counts <- function(ccs, prop=0.8) {
  counts(ccs) = scuttle::downsampleMatrix(counts(ccs), prop=prop)
  return(ccs)
}

read_result_files <- function(result_files, result_dir) {

  res_df = lapply(result_files, function(f){

    f = gsub(f, pattern = ".csv", replacement = "")
    string_split = stringr::str_split(f, pattern = "_|-") %>% unlist()
    which_sim = string_split[4]
    embryo_size = string_split[7]
    cell_grouping = paste0("cell_type_", string_split[10])
    prop_ct = string_split[13]

    data.table::fread(paste0(result_dir, f, ".csv"), sep = ",", stringsAsFactors = F, data.table = F, na.strings = "") %>%
      mutate("which_sim" = which_sim, "embryo_size" = embryo_size, "prop_ct" = prop_ct, "cell_grouping" = cell_grouping)

  }) %>% bind_rows()

  return(res_df)
}


prepare_res_df <- function(res_df, wt_sim_props) {
  sims = seq(1, 10, 1)
  emb_interval = 5
  embs = sims*emb_interval

  props = wt_sim_props %>%
    mutate(cell_type = paste0("cell_type_", cell_type)) %>%
    group_by(cell_type) %>%
    mutate(ct_means = mean(cell_prop)) %>%
    distinct(cell_type, ct_means)

  res_df = res_df %>%
    left_join(props, by = "cell_type") %>%
    mutate(cell_type_prop = paste("prop_", round(ct_means, 3), sep = ""))

  # order cell types by proportion and reset levels for plotting
  ct_order = res_df %>%
    select(cell_type_prop, ct_means) %>%
    arrange(-ct_means) %>%
    distinct(cell_type_prop) %>%
    pull(cell_type_prop)

  res_df = res_df %>%
    left_join(data.frame("which_sim" = as.character(sims), "emb_num" = embs), by = "which_sim")

  res_df$cell_type_prop = factor(res_df$cell_type_prop, levels = ct_order)
  res_df = res_df %>%
    mutate(is_sig = case_when(delta_p_value < 0.05 ~ TRUE,
                              delta_p_value > 0.05 ~ FALSE))

  num_cell_counts = res_df %>%
    group_by(emb_num, effect_size, embryo_size, prop_ct, which_sim, cell_group) %>%
    tally() %>%
    group_by(prop_ct, embryo_size) %>% summarise("num_cell_types" = mean(n))

  res_df = left_join(res_df, num_cell_counts, by = c("prop_ct", "embryo_size"))

  res_df$embryo_size = factor(res_df$embryo_size,
                              levels = unique(res_df$embryo_size) %>% as.numeric %>% sort() %>% as.character())
  return(res_df)
}


calculate_TPR <- function(res_df) {

  sum_df = res_df %>%
    filter(!is.na(cell_group)) %>% # idk why this is here
    mutate(type = case_when(
      (cell_type == cell_group) & (delta_p_value < q_value_threshold) ~ "TP",
      (cell_type == cell_group) & (delta_p_value >= q_value_threshold) ~ "FN",
      (cell_type != cell_group) & (delta_p_value < q_value_threshold) ~ "FP",
      (cell_type != cell_group) & (delta_p_value >= q_value_threshold) ~ "TN",
    )) %>%
    # group_by(cell_type, type) %>%
    group_by(cell_type, effect_size, embryo_size, emb_num, which_sim, prop_ct, type) %>%
    tally() %>%
    pivot_wider(names_from = type, values_from = n) %>%
    dplyr::union_all(dplyr::tibble(TP = integer(), FP = integer(),
                                   TN = integer(), FN = integer())) %>%
    replace(is.na(.), 0) %>%
    mutate(TPR = TP / (TP + FN),
           FPR = FP / (FP + TN))

  return(sum_df)

}


calc_tpr = function(df, q_value_threshold) {

  # missing_values =
  # df %>%
  #   group_by(cell_type, embryo_size, rep, method) %>%
  #   tally() %>%
  #   pivot_wider(values_from = n, names_from = method) %>%
  #   replace(is.na(.), 0) #%>%
  #   mutate_at(c("hooke", "hooke_bootstrap", "propeller", "bb"), ~.-87) %>%
  #   ungroup %>%
  #   group_by(cell_type, embryo_size) %>%
  #   summarise_at(c("hooke", "hooke_bootstrap", "propeller", "bb"), sum) %>%
  #   pivot_longer(-c(cell_type, embryo_size), values_to = "missing", names_to = "method")

  df  %>%
    mutate(type = case_when(
      (cell_type == as.character(cell_group)) & (p.value < q_value_threshold) ~ "TP",
      (cell_type == as.character(cell_group)) & (p.value >= q_value_threshold) ~ "FN",
      (cell_type != as.character(cell_group)) & (p.value < q_value_threshold) ~ "FP",
      (cell_type != as.character(cell_group)) & (p.value >= q_value_threshold) ~ "TN",
      (cell_type == as.character(cell_group)) & (is.na(p.value)) ~ "FN",
      (cell_type != as.character(cell_group)) & (is.na(p.value)) ~ "TN"
    )) %>%
    group_by(cell_type, type, effect_size, embryo_size, method) %>%
    tally() %>%
    pivot_wider(names_from = type, values_from = n) %>%
    dplyr::union_all(dplyr::tibble(TP = integer(), FP = integer(),
                                   TN = integer(), FN = integer())) %>%
    replace(is.na(.), 0) %>%
    # left_join(missing_values, by = c("cell_type", "embryo_size", "method")) %>%
    # mutate(FN = FN - missing) %>% select(-missing) %>%
    mutate(TPR = TP / (TP + FN),
           FPR = FP / (FP + TN)) %>%
    replace(is.na(.), 0) %>% as.data.frame()

}
