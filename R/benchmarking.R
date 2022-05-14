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

# run the beta binomial test
# to compare the values
# return output similar to compare_abundances
run_beta_binomial <- function(ccs,
                          model_formula = "count_df ~ genotype",
                          covariates = c("experiment", "gene_target", "timepoint", "knockout"),
                          ...) {


  count_df = counts(ccs) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("cell_group") %>%
    tidyr::pivot_longer(-"cell_group", names_to = "sample", values_to = "cells")

  coldata = ccs@colData %>%
    as.data.frame %>%
    select(all_of(covariates), sample) %>%
    distinct()
  rownames(coldata) = NULL

  comb_df = left_join(count_df,
            coldata,
            by = "sample") %>%
    group_by(sample) %>%
    mutate(total_cells = sum(cells)) %>%
    ungroup() %>%
    # i think you have to remove 0s?
    filter(cells > 0 )

  wt_df = comb_df %>% filter(knockout == FALSE)
  comb_df = comb_df %>%
    mutate(genotype = forcats::fct_relevel(comb_df$gene_target, c(unique(wt_df$gene_target))))

  cell.groups = unique(comb_df$cell_group)

  test_res = sapply(cell.groups,
                    FUN = function(x) {
                      type_df = comb_df %>% filter(cell_group == x)
                      count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
                      fit =  VGAM::vglm(as.formula(model_formula), "betabinomial", data = type_df, trace = TRUE, ...)
                      fit_df = tidy.vglm(fit)}, USE.NAMES = T, simplify = F)

  test_res = do.call(rbind, test_res)
  test_res = test_res %>% tibble::rownames_to_column(var = "cell_group")
  test_res %>% arrange(desc(estimate))

  return(test_res)
}


fit_beta_binomial = function(ccs,
                             genotype,
                             timepoint,
                             ctrl_ids = c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa", "ctrl-tbx16", "ctrl-met"),
                             ...) {

  # subset the ccs to gene target and timepoint
  subset_ccs = ccs[,colData(ccs)$gene_target == genotype | colData(ccs)$gene_target %in% ctrl_ids]
  subset_ccs = subset_ccs[,colData(subset_ccs)$timepoint == timepoint]
  colData(subset_ccs)$knockout = (colData(subset_ccs)$gene_target == genotype)

  test_res = run_beta_binomial(subset_ccs)
  test_res = test_res %>%
    filter(!(grepl("intercept", term, ignore.case = T))) %>%
    dplyr::mutate(qval = p.adjust(p.value, method = "BH")) %>%
    mutate(geno.v.wt = stringr::str_sub(term, 9)) %>%
    tidyr::separate(cell_group, into = c("cell_group", NULL), sep = "\\.") %>%
    select(geno.v.wt, everything(), -term) %>%
    arrange(geno.v.wt) %>%
    mutate(timepoint = timepoint)

  # make it match column names for ease of plotting
  test_res %>% select(cell_group, "delta_log_abund" = estimate, "delta_q_value" = qval)

}

# want to run this

# genotype_df %>% mutate(timepoint = 18) %>%
#   dplyr::mutate(bb_res = purrr::map2(.f = purrr::possibly(fit_beta_binomial, NA_real_),
#                                     .x = gene_target,
#                                     .y = timepoint,
#                                     ccs = ccs,
#                                     ctrl_ids=c("wt", "ctrl-inj", "ctrl-noto", "ctrl-mafba", "ctrl-hgfa")))


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











