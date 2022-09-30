collect_pln_gene_edges <- function(ccm,
                                   cond_b_vs_a_tbl,
                                   log_abundance_thresh = 1-5) {

  pln_model = model(ccm, "reduced")
  abundance_corr_tbl = correlate_abundance_changes(pln_model, cond_b_vs_a_tbl)

  abundance_corr_tbl = abundance_corr_tbl %>% dplyr::filter(
    (to_log_abund_x > log_abundance_thresh | to_log_abund_y > log_abundance_thresh) &  # Keep if the "to" node is above abundance thresh in at least one condition
      (from_log_abund_x > log_abundance_thresh | from_log_abund_y > log_abundance_thresh) # Keep if the "to" node is above abundance thresh in at least one condition
  )

  corr_edge_delta_abund = abundance_corr_tbl

  corr_edge_delta_abund = corr_edge_delta_abund %>% dplyr::mutate(scaled_weight = -pcor)

  # corr_edge_delta_abund = corr_edge_delta_abund %>%
  #   dplyr::mutate(edge_type = dplyr::case_when(
  #     # pos pcor, same fold change
  #     from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor > 0 ~ "activator",
  #     from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
  #     # pos pcor, reciprocal fold change
  #     from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
  #     from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
  #     # neg pcor, same fold change
  #     from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
  #     from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor < 0 ~ "undirected",
  #     # neg pcor, reciprocal fold change
  #     from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor < 0 ~ "repressor",
  #     from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
  #     TRUE ~ "hidden"
  #   ))

  return (corr_edge_delta_abund)
}


comp_abund <- function(to, from, ccm) {
  to_abund = estimate_abundances(ccm, tibble("cell_state" = to))
  from_abund = estimate_abundances(ccm, tibble("cell_state" = from))
  return(compare_abundances(ccm, to_abund, from_abund))
}


id_to_shortname = function(cds, gene_ids) {
  rowData(cds) %>%
    as.data.frame() %>%
    select(gene_short_name, id) %>%
    filter(id %in% gene_ids) %>%
    pull(gene_short_name)
}

shortname_to_id = function(cds, short_names) {
  rowData(cds) %>%
    as.data.frame() %>%
    select(gene_short_name, id) %>%
    filter(gene_short_name %in% short_names) %>%
    pull(id)
}


# directed edges
# fc direction
# pcor act/rep

find_reg_states <- function(path,
                            from_state,
                            to_state,
                            upreg = T,
                            state_1 = T,
                            log_abund_threshold = 0,
                            delta_p_value_threshold = 1.0) {


  path_edges = path %>%
    filter(from == from_state, to == to_state) %>%
    select(filter_edges) %>%
    tidyr::unnest(filter_edges) %>%
    filter(to_delta_p_value < delta_p_value_threshold |
             from_delta_p_value < delta_p_value_threshold)

  # either get upreg or downreg states
  if (upreg) {
    path_edges = path_edges %>%
      filter(to_delta_log_abund > log_abund_threshold |
               from_delta_log_abund > log_abund_threshold)
  } else {
    path_edges = path_edges %>%
      filter(to_delta_log_abund < log_abund_threshold |
             from_delta_log_abund < log_abund_threshold)
  }

  # is this state the first state (from) or the second state (to)
  if (state_1) {
    path_edges = path_edges %>%
      select_if(grepl("from", colnames(path_edges)) | colnames(path_edges) %in% c("from", "to", "pcor"))
  } else {
    path_edges = path_edges %>%
      select_if(grepl("to", colnames(path_edges)) | colnames(path_edges) %in% c("from", "to", "pcor"))
  }

  return(path_edges)

}

find_all_edges = function(path, sub_states, ...) {

  # both upregulated from 1 to 2 and 2 to 3
  upreg_1_to_2 = find_reg_states(path, from = sub_states[1], to = sub_states[2], upreg = T, state_1 = T, ...)
  upreg_2_to_3 = find_reg_states(path, from = sub_states[2], to = sub_states[3], upreg = T, state_1 = F, ...)
  downreg_1_to_2 = find_reg_states(path, from = sub_states[1], to = sub_states[2], upreg = F, state_1 = T, ...)
  downreg_2_to_3 = find_reg_states(path, from = sub_states[2], to = sub_states[3], upreg = F, state_1 = F, ...)

  join_df = function(df1, df2) {
    joint_df = inner_join(df1, df2,
               by = c("from", "to"),
               suffix = c(".from", ".to")) %>%
      rename("pcor" = "pcor.from") %>%
      select(-c(pcor.to))
    return(joint_df)
  }

  # both upreg
  both_upreg_edges = join_df(upreg_1_to_2, upreg_2_to_3)
  # both downreg
  both_downreg_edges = join_df(downreg_1_to_2, downreg_2_to_3)
  # reciprocal changes
  upreg_downreg_edges = join_df(upreg_1_to_2, downreg_2_to_3)
  downreg_upreg_edges = join_df(downreg_1_to_2, upreg_2_to_3)


  # stick them together
  all_edges = rbind(both_upreg_edges,
                    both_downreg_edges,
                    upreg_downreg_edges,
                    downreg_upreg_edges) %>%
    dplyr::mutate(scaled_weight  = abs(pcor) / max(abs(pcor)))

  all_edges = all_edges %>%
    dplyr::mutate(edge_type = dplyr::case_when(
        # pos pcor, same fold change
        from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor > 0 ~ "activator",
        from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
        # pos pcor, reciprocal fold change
        from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor > 0 ~ "undirected",
        from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor > 0 ~ "undirected",
        # neg pcor, same fold change
        from_delta_log_abund > 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
        from_delta_log_abund < 0 & to_delta_log_abund < 0 & pcor < 0 ~ "undirected",
        # neg pcor, reciprocal fold change
        from_delta_log_abund > 0 & to_delta_log_abund < 0 & pcor < 0 ~ "repressor",
        from_delta_log_abund < 0 & to_delta_log_abund > 0 & pcor < 0 ~ "undirected",
        TRUE ~ "hidden"
      ))

  return(all_edges)
}
