# this is for all benchmarking functions


load_lineage_tree <- function() {
  # fix to load from gdrive or github?
  G <- readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/lineage_map_igraph_210824.rds")
  # edges = igraph::as_data_frame(G, what = "edges")
  # nodes = igraph::as_data_frame(G, what = "vertices") #%>%
  # #   rename("cell_type_broad"="name") %>%
  # #   left_join(anno_df, by = "cell_type_broad") %>%
  # #   filter(cell_type_broad != "neuron (cranial ganglion)")
  #
  # # no_edge = setdiff(nodes$cell_type_broad, union(edges$from, edges$to))
  # # edges = rbind(edges, data.frame("from" = no_edge, "to" = no_edge))
  #
  # n = network::network(edges, directed = T,loops = T)
  # nodes = nodes[match(network.vertex.names(n), nodes$cell_type_broad),]
  # n %v% "id" = network.vertex.names(n)
  # n %v% "tissue" = nodes[["tissue"]]
  # n %v% "cell_type_broad" = nodes[["cell_type_broad"]]
  # n %v% "cell_type_sub" = nodes[["cell_type_sub"]]

  return(G)
}

# lit_tree <- load_lineage_tree()

# a help function to tidy the vgam model output - used in compare_abundance function
tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
  co <- as.data.frame(coef(VGAM::summary(x)))
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

# this shouldn't change size factors
subset_ccs = function(ccs, ctrl_ids, gene_target, time) {
  gene_ids = c(ctrl_ids, gene_target)
  ccs = ccs[, colData(ccs)$timepoint == time]
  ccs = ccs[, colData(ccs)$gene_target %in% gene_ids]
}


bb_compare_abundance <- function(ccs,
                                 comp_col,
                                 model_formula,
                                 ctrl_ids,
                                 nuisance_cols = NULL,
                                 nuisance_formula = NULL,
                                 ...) {

  # make ccs
  colData(ccs)[[comp_col]] = as.data.frame(colData(ccs)) %>%
    mutate(knockout = ifelse(!!sym(comp_col) %in% ctrl_ids, "ctrl", !!sym(comp_col))) %>%
    pull(knockout)

  ccs_coldata = colData(ccs) %>% as.data.frame %>% select(sample, !!sym(comp_col),
                                                          Size_Factor, all_of(nuisance_cols))

  count_df = counts(ccs) %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_group") %>%
    tidyr::pivot_longer(-"cell_group", names_to = "sample", values_to = "cells") %>%
    group_by(sample) %>%
    left_join(ccs_coldata, by="sample") %>%
    mutate(cells = round(cells/Size_Factor)) %>%
    mutate(total_cells = sum(cells))

  fc_df = count_df %>%
    group_by_at(vars(comp_col, cell_group, nuisance_cols)) %>%
    summarize(cell_mean = mean(cells)) %>%
    tidyr::pivot_wider(names_from = comp_col, values_from = "cell_mean") %>%
    tidyr::pivot_longer(-c(cell_group, ctrl, nuisance_cols), names_to = comp_col, values_to = "cells") %>%
    mutate(abund_log2fc = log2((cells + 1)/(ctrl+1))) %>%
    arrange(!!sym(comp_col))


  # to fix new col name
  count_df[[comp_col]] = forcats::fct_relevel(count_df[[comp_col]], "ctrl")

  cell.groups = unique(count_df$cell_group)

  # a help function to tidy the vgam model output - used in compare_abundance function
  tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
    co <- as.data.frame(coef(VGAM::summary(x)))
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

  fit_beta_binomial = function(cg, count_df, model_formula, trace = TRUE, ...) {
    type_df = count_df %>% filter(cell_group == cg)
    count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
    fit =  VGAM::vglm(as.formula(model_formula), "betabinomial", data = type_df, trace = trace, ...)
    fit_df = tidy.vglm(fit)

  }

  model_formula = stringr::str_replace_all(model_formula, "~", "")
  model_formula_str = paste0("count_df ~", model_formula)

  if (is.null(nuisance_formula) == FALSE) {
    nuisance_formula = stringr::str_replace_all(nuisance_formula, "~", "")
    model_formula_str = paste0(model_formula_str, "+", nuisance_formula)
  }

  bb_res = data.frame("cell_group" = cell.groups) %>%
    mutate("beta_binomial" = purrr::map(.f  = purrr::possibly(fit_beta_binomial, NA_character_),
                                        .x = cell_group,
                                        count_df = count_df,
                                        model_formula = model_formula_str)) %>%
    filter(!is.na(beta_binomial)) %>%
    tidyr::unnest(c(beta_binomial)) %>%
    arrange(desc(estimate)) %>%
    filter(!grepl("Intercept", term))

  bb_res = left_join(bb_res,
                     fc_df %>% select(cell_group, abund_log2fc, "ctrl_mean" = "ctrl"),
                     by = "cell_group")

  return(bb_res)
}


# it's lost because its parent is lost
# it's lost because some other cell type it depends on is lost (i.e. signaling)
# it's lost because it directly depends on the gene


parent_child_foldchanges <- function(state_graph, cond_b_v_a_tbl, genes = list(), qval_threshold = 0.05) {

  cond_b_v_a_tbl = cond_b_v_a_tbl %>%
    mutate(delta_log_abund = ifelse(delta_q_value < qval_threshold, delta_log_abund, 0))

  if (length(genes) > 0) {
    gene_expr = aggregated_expr_data(ccm@ccs@cds, group_cells_by = ccm@ccs@info$cell_group)
    sub_gene_expr = gene_expr %>%
      filter(gene_short_name %in% genes)
    cond_b_v_a_tbl = cond_b_v_a_tbl %>% left_join(sub_gene_expr, by = c("cell_group"))

  }

  state_graph_edges = igraph::as_data_frame(state_graph, what = "edges") %>%
    dplyr::rename("parent"="from", "child"="to") %>%
    left_join(cond_b_v_a_tbl, by = c("parent" = "cell_group")) %>%
    left_join(cond_b_v_a_tbl, by = c("child" = "cell_group"), suffix = c(".parent", ".child")) %>%
    mutate(fold_changes = case_when(
      delta_log_abund.parent < 0 & delta_log_abund.child < 0 ~ "parent decrease, descendents decrease",
      delta_log_abund.parent < 0 & delta_log_abund.child == 0 ~ "parent decrease, descendents no change",
      delta_log_abund.parent < 0 & delta_log_abund.child > 0 ~ "parent decrease, descendents decrease",
      delta_log_abund.parent > 0 & delta_log_abund.child > 0 ~ "parent increase, descendents increase",
      delta_log_abund.parent > 0 & delta_log_abund.child == 0 ~ "parent increase, descendents no change",
      delta_log_abund.parent > 0 & delta_log_abund.child < 0 ~ "parent increase, descendents decrease",
      TRUE ~ "no change in either"
    ))

  if (length(genes) > 0) {
  state_graph_edges = state_graph_edges %>%
    mutate(gene_expr = case_when(
      fraction_expressing.parent >= 0.1 & fraction_expressing.child >= 0.1 ~ "both express",
      fraction_expressing.parent >= 0.1 & fraction_expressing.child < 0.1 ~ "parent expr, child doesn't expr",
      fraction_expressing.parent < 0.1 & fraction_expressing.child >= 0.1 ~ "parent doesn't expr, child expr",
      fraction_expressing.parent < 0.1 & fraction_expressing.child < 0.1 ~ "neither expresses",
    ))
  }

  state_graph_edges %>%
    select(parent, child, fold_changes) %>%
    group_by(parent, fold_changes) %>%
    tally()

}








