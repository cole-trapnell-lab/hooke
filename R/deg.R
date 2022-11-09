

# Calculate the probability vector
makeprobsvec <- function(p) {
  phat <- p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

# Calculate the probability matrix for a relative abundance matrix
makeprobs <- function(a) {
  colSums<-apply(a,2,sum)
  b <- Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# Calculate the Shannon entropy based on the probability vector
# shannon.entropy <- function(p) {
#   if (min(p) < 0 || (p) <=0)
#     return(Inf)
#   p.norm <- p[p>0]/sum(p)
#   -sum(log2(p.norm)*p.norm)
# }

shannon_entropy <- function(p) {
  #if (Matrix::rowMin(p) < 0 || (p) <=0)
  #  return(Inf)
  p.norm <- p / Matrix::rowSums(p)
  lg_pnorm = log2(p.norm) * p.norm
  lg_pnorm[p.norm == 0] = 0
  SE = -Matrix::rowSums(lg_pnorm)
  return (SE)
}

# Calculate the Jensen-Shannon distance for two probability distribution
js_dist_to_pattern <- function (x, pattern)
{
  stopifnot(ncol(x) == length(pattern))
  avg_x_pattern = sweep(x, 2, pattern, "+") / 2

  JSdiv = shannon_entropy(avg_x_pattern) -
    (shannon_entropy(x) + shannon_entropy(matrix(pattern, nrow=1))) * 0.5
  JSdiv[is.infinite(JSdiv)] = 1
  JSdiv[JSdiv < 0] = 0
  JSdist <- sqrt(JSdiv)
  pattern_match_score = 1 - JSdist
  return(pattern_match_score)
}

score_genes_for_expression_pattern <- function(cell_state, gene_patterns, state_graph, estimate_matrix, state_term="cell_group", cores=1){

  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))

  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))

  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))

  #states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()

  expr_df = tibble(gene_id=row.names(estimate_matrix))



  self_and_parent = exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents), drop=FALSE])
  self_parent_sibs = exp(estimate_matrix[gene_patterns$gene_id, c(cell_state, parents, siblings), drop=FALSE])

  num_parents = length(parents)
  num_siblings = length(siblings)
  gene_patterns_and_scores = gene_patterns %>%
    tidyr::unnest(interpretation) %>%
    #group_by(interpretation) %>%
    mutate(pattern_match_score =
             case_when(interpretation %in% c("Absent") ~ 0,
                       interpretation %in% c("Maintained") ~ js_dist_to_pattern(self_and_parent,c(1, rep(1, num_parents))),
                       interpretation %in% c("Selectively maintained", "Specifically maintained") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(1, num_parents), rep(0, num_siblings))),
                       interpretation %in% c("Upregulated", "Activated") ~ js_dist_to_pattern(self_and_parent, c(1, rep(0, num_parents))),
                       interpretation %in% c("Selectively upregulated", "Specifically upregulated", "Selectively activated", "Specifically activated") ~ js_dist_to_pattern(self_parent_sibs, c(1, rep(0, num_parents), rep(0, num_siblings))),
                       interpretation %in% c("Downregulated", "Dectivated") ~ js_dist_to_pattern(self_and_parent, c(0, rep(1, num_parents))),
                       interpretation %in% c("Selectively downregulated", "Specifically downregulated", "Selectively deactivated", "Specifically deactivated") ~ js_dist_to_pattern(self_parent_sibs, c(0, rep(1, num_parents), rep(0, num_siblings))),
                       TRUE ~ 0)
           )
  return(gene_patterns_and_scores)
}


classify_genes_in_cell_state <- function(cell_state, state_graph, estimate_matrix, stderr_matrix, state_term="cell_group", log_fc_thresh=1, abs_expr_thresh = 1e-3, sig_thresh=0.05, cores=1){
  #expr_self = expr_mat[,cell_state]

  parents = get_parents(state_graph, cell_state) #igraph::neighbors(state_graph, cell_state, mode="in")
  parents = intersect(parents, colnames(estimate_matrix))

  children = get_children(state_graph, cell_state)#igraph::neighbors(state_graph, cell_state, mode="out")
  children = intersect(children, colnames(estimate_matrix))

  siblings = get_siblings(state_graph, cell_state)#igraph::neighbors(state_graph, parents, mode="out")
  siblings = intersect(siblings, colnames(estimate_matrix))

  states_in_contrast = c(cell_state, parents, children, siblings) %>% unique()

  expr_df = tibble(gene_id=row.names(estimate_matrix))

  message("      examining coeffficients ", cell_state)

  expr_df$expr_self = pnorm(estimate_matrix[,cell_state] - log(abs_expr_thresh), sd = stderr_matrix[,cell_state], lower.tail=FALSE)
  expr_df$expr_self = p.adjust(expr_df$expr_self, method="BH") < sig_thresh

  expr_df$expressed_in_parents = NA
  expr_df$expressed_in_siblings = NA
  expr_df$higher_than_parents = NA
  expr_df$lower_than_parents = NA
  expr_df$higher_than_all_siblings = NA
  expr_df$lower_than_all_siblings = NA
  expr_df$higher_than_siblings = NA
  expr_df$lower_than_siblings = NA

  if (length(parents) > 0){
    expressed_in_parents_mat = pnorm(estimate_matrix[,parents, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,parents, drop=F], lower.tail=FALSE)
    expressed_in_parents_mat = apply(expressed_in_parents_mat, 2, p.adjust, method="BH")

    expressed_in_parents_mat = expressed_in_parents_mat < sig_thresh
    expr_df$expressed_in_parents = Matrix::rowSums(expressed_in_parents_mat) > 0

    higher_than_parents_stat = -t(sweep(t(estimate_matrix[,parents, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_parents_pval = pnorm(higher_than_parents_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_parents_pval = apply(higher_than_parents_pval, 2, p.adjust, method="BH")

    higher_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & higher_than_parents_pval < sig_thresh
    expr_df$higher_than_parents = Matrix::rowSums(higher_than_parents_mat) > 0

    lower_than_parents_pval = pnorm(-higher_than_parents_stat,
                                    sd = sqrt(sweep(t(stderr_matrix[,parents, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_parents_pval = apply(lower_than_parents_pval, 2, p.adjust, method="BH")

    lower_than_parents_mat = abs(higher_than_parents_stat) > log_fc_thresh & lower_than_parents_pval < sig_thresh
    expr_df$lower_than_parents = Matrix::rowSums(lower_than_parents_mat) > 0
  }else{
    expr_df$expressed_in_parents = NA
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_parents = NA
    expr_df$lower_than_parents = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA

  }

  if (length(siblings) > 0){
    expressed_in_siblings_mat = pnorm(estimate_matrix[,siblings, drop=F] - log(abs_expr_thresh), sd = stderr_matrix[,siblings, drop=F], lower.tail=FALSE)
    expressed_in_siblings_mat = apply(expressed_in_siblings_mat, 2, p.adjust, method="BH")

    expressed_in_siblings_mat = expressed_in_siblings_mat < sig_thresh
    expr_df$expressed_in_siblings = Matrix::rowSums(expressed_in_siblings_mat) > 0

    higher_than_siblings_stat = -t(sweep(t(estimate_matrix[,siblings, drop=F]), 2, as.numeric(estimate_matrix[,cell_state]) , `-`))
    higher_than_siblings_pval = pnorm(higher_than_siblings_stat,
                                      sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    higher_than_siblings_pval = apply(higher_than_siblings_pval, 2, p.adjust, method="BH")

    higher_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & higher_than_siblings_pval < sig_thresh
    expr_df$higher_than_all_siblings = Matrix::rowSums(higher_than_siblings_mat) == ncol(higher_than_siblings_pval)
    expr_df$higher_than_siblings = Matrix::rowSums(higher_than_siblings_mat) > 0

    lower_than_siblings_pval = pnorm(-higher_than_siblings_stat,
                                     sd = sqrt(sweep(t(stderr_matrix[,siblings, drop=F]^2), 2, as.numeric(stderr_matrix[,cell_state, drop=F]^2), `+`)), lower.tail=FALSE)
    lower_than_siblings_pval = apply(lower_than_siblings_pval, 2, p.adjust, method="BH")

    lower_than_siblings_mat = abs(higher_than_siblings_stat) > log_fc_thresh & lower_than_siblings_pval < sig_thresh
    expr_df$lower_than_all_siblings = Matrix::rowSums(lower_than_siblings_mat) == ncol(lower_than_siblings_mat)
    expr_df$lower_than_siblings = Matrix::rowSums(lower_than_siblings_mat) > 0


  }else{
    expr_df$expressed_in_siblings = NA
    expr_df$higher_than_all_siblings = NA
    expr_df$lower_than_all_siblings = NA
    expr_df$higher_than_siblings = NA
    expr_df$lower_than_siblings = NA
  }

  expr_df = expr_df %>% tidyr::nest(data = !gene_id)

  message("      interpreting patterns")
  interpret_expression_pattern = function(pat_df){
    if (pat_df$expr_self){
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Maintained")
      }else if (pat_df$expressed_in_parents){
        # Expressed in self and parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self and parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Upregulated")
          else if(pat_df$lower_than_parents)
            return("Downregulated")
          else
            return("Maintained")
        } else {
          # Expressed in self and parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than siblings
              return("Specifically upregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than siblings
              return("Selectively upregulated")
            }
            else if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than siblings
              return("Upregulated")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
          else if(pat_df$lower_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than siblings
              return("Selectively downregulated")
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Selectively downregulated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than  all siblings
              return("Specifically downregulated")
            }
            else { # same as parent, same as siblings
              return("Downregulated")
            }
          }
          else { # same as parent
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Specifically maintained")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Selectively maintained")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than siblings
              return("Maintained")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
        }

      }else{
        # Expressed in self but not in parent
        if (is.na(pat_df$expressed_in_siblings)){
          # Expressed in self, not in parent and there are no siblings
          if (pat_df$higher_than_parents)
            return("Activated")
          else if(pat_df$lower_than_parents)
            return("Downregulated") # shouldn't happen
          else
            return("Maintained") # shouldn't happen
        } else {
          # Expressed in self and not in parent and there are siblings
          if (pat_df$higher_than_parents){
            if (pat_df$expressed_in_siblings == FALSE | pat_df$higher_than_all_siblings){
              # Higher than parent, and higher than all siblings
              return("Specifically activated")
            } else if (pat_df$higher_than_siblings){
              # Higher than parent, and higher than some siblings
              return("Selectively activated")
            }
            else if(pat_df$lower_than_all_siblings){
              # Higher than parent, but lower than all siblings
              return("Activated")
            }
            if(pat_df$lower_than_siblings){
              # Higher than parent, but lower than some siblings
              return("Activated")
            }
            else { # same as parent, same as siblings
              return("Maintained")
            }
          }
          else if(pat_df$lower_than_parents){
            # if the gene is lower in the parent, which is off, just mark the gene absent
            if (pat_df$higher_than_all_siblings){
              # Lower than parent, and higher than all siblings
              return("Absent") # shouldn't happen
            }
            else if (pat_df$higher_than_siblings){
              # Lower than parent, and higher than some siblings
              return("Absent") # shouldn't happen
            }
            else if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
          else { # same as parent (which is off)
            if (pat_df$higher_than_all_siblings){
              # Same as parent, and higher than all siblings
              return("Absent")
            }
            else if (pat_df$higher_than_siblings){
              # Same as parent, and higher than some siblings
              return("Absent")
            }
            else if(pat_df$lower_than_all_siblings){
              # Same as parent, but lower than all siblings
              return("Absent")
            }
            else if(pat_df$lower_than_siblings){
              # Same as parent, but lower than some siblings
              return("Absent")
            }
            else { # same as parent, same as siblings
              return("Absent")
            }
          }
        }
      }
      return ("Expressed")
    }else{
      # Not expressed in self
      if (is.na(pat_df$expressed_in_parents)){
        # no parents, therefore no siblings
        return ("Absent")
      }else if (pat_df$expressed_in_parents){
        # Not expressed in self, but expressed in parents
        if (is.na(pat_df$expressed_in_siblings)){
          # Not expressed in self, expressed parent and there are no siblings
          if(pat_df$lower_than_parents)
            return("Deactivated")
          else
            return("Absent") # shouldn't happen
        } else {
          # Not expressed in self, expressed in parent and there are siblings
          if(pat_df$lower_than_parents){
            # Lower than parent
            if(pat_df$lower_than_all_siblings){
              # Lower than parent and  lower than siblings
              return("Specifically deactivated")
            }
            else if(pat_df$lower_than_siblings){
              # Lower than parent and  lower than siblings
              return("Selectively deactivated")
            }
            return("Deactivated")
          }
          else {
            #Not expressed in self, not lower than parent
            return ("Absent")
          }
        }
      }else{
        # Not expressed in self or parents
        return ("Absent")
      }
      return ("Absent")
    }
    return ("Absent")
    #match_row = match(data.frame(t(pat_df)), data.frame(t(interp_table)))
    #interpetation[match_row]
  }
  #debug(interpret_expression_pattern)
  expr_df = expr_df %>% mutate(interpretation = purrr::map(.f = purrr::possibly(
    interpret_expression_pattern, NA_character_), .x = data))
  message("      completed ", cell_state)
  return(expr_df)
}
#debug(classify_genes_in_cell_state)


#' Classify each gene's pattern of expresison in each state in a state transition graph
#' @export
classify_genes_over_graph <- function(ccm,
                                      state_graph,
                                      gene_ids = NULL,
                                      group_nodes_by=NULL,
                                      log_fc_thresh=1,
                                      abs_expr_thresh = 1e-3,
                                      sig_thresh=0.05,
                                      min_samples_detected = 2,
                                      min_cells_per_pseudobulk = 3,
                                      cores=1,
                                      ...){
  if (is.null(group_nodes_by)){
    pb_cds = pseudobulk_cds_for_states(ccm)
    state_term = "cell_group"
  }else{
    pb_cds = pseudobulk_cds_for_states(ccm, state_col = group_nodes_by)
    state_term = group_nodes_by
  }

  if (!is(state_graph, "igraph")){
    state_graph = state_graph %>% igraph::graph_from_data_frame()
  }


  #cds_to_test = pb_cds[,as.character(colData(pb_cds)[,state_term]) %in% states_in_model]

  #colData(cds_to_test)[,state_term] = factor(as.character(colData(cds_to_test)[,state_term]), levels=states_in_model) # set the "self" state as the reference level

  #norm_expr_mat = normalized_counts(pb_cds, "size_only", pseudocount = 0)

  if (is.null(gene_ids) == FALSE){
    pb_cds = pb_cds[gene_ids,]
  }

  # expr_over_thresh = threshold_expression_matrix(normalized_counts(pb_cds, "size_only", pseudocount = 0), ...)
  expr_over_thresh = normalized_counts(pb_cds, "size_only", pseudocount = 0)
  genes_to_test = which(Matrix::rowSums(expr_over_thresh) >= min_samples_detected)
  pb_cds = pb_cds[genes_to_test,]

  pseudobulks_to_test = which(colData(pb_cds)$num_cells_in_group >= min_cells_per_pseudobulk)

  message("fitting regression models")
  pb_cds = pb_cds[,pseudobulks_to_test]
  pb_group_models = fit_models(pb_cds,
                               model_formula_str=paste("~ 0 + ", state_term),
                               weights=colData(pb_cds)$num_cells_in_group,
                               cores=cores) %>% dplyr::select(gene_short_name, id, model, model_summary)

  message("      collecting coefficients")
  pb_group_models = coefficient_table(pb_group_models) %>%
    dplyr::select(gene_short_name, id, term, estimate, std_err) %>%
    mutate(term = stringr::str_replace_all(term, state_term, ""))
  estimate_matrix = pb_group_models %>% dplyr::select(id, term, estimate)
  estimate_matrix = estimate_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  estimate_matrix = estimate_matrix %>% tidyr::pivot_wider(names_from=term, values_from=estimate, values_fill=0)

  gene_ids = estimate_matrix$id
  estimate_matrix$id = NULL
  estimate_matrix = as.matrix(estimate_matrix)
  row.names(estimate_matrix) = gene_ids
  colnames(estimate_matrix) = as.character(colnames(estimate_matrix))

  stderr_matrix = pb_group_models %>% dplyr::select(id, term, std_err)
  stderr_matrix = stderr_matrix %>% mutate(term = factor(term, levels=unique(colData(pb_cds)[,state_term])))
  stderr_matrix = stderr_matrix %>% tidyr::pivot_wider(names_from=term, values_from=std_err, values_fill=0)

  gene_ids = stderr_matrix$id
  stderr_matrix$id = NULL
  stderr_matrix = as.matrix(stderr_matrix)
  row.names(stderr_matrix) = gene_ids
  colnames(stderr_matrix) = as.character(colnames(stderr_matrix))

  #p_val_matrix = pnorm(estimate_matrix - log(abs_expr_thresh), sd = stderr_matrix, lower.tail=FALSE)

  #expr_thresh_mat = p_val_matrix < sig_thresh

  #cell_states = tibble(cell_state = unlist(igraph::V(state_graph)$name))
  states_to_assess = intersect(as.character(unique(colData(pb_cds)[,state_term])), unlist(igraph::V(state_graph)$name))
  cell_states = tibble(cell_state = states_to_assess)

  cell_states = cell_states %>%
    dplyr::mutate(gene_classes = purrr::map(.f = purrr::possibly(
      classify_genes_in_cell_state, NA_real_), .x = cell_state,
      state_graph, estimate_matrix, stderr_matrix, state_term,
      log_fc_thresh=log_fc_thresh,
      abs_expr_thresh = abs_expr_thresh,
      sig_thresh=sig_thresh,
      cores=cores))

  cell_states = cell_states %>%
    filter(is.na(gene_classes) == FALSE) %>%
    dplyr::mutate(gene_class_scores = purrr::map2(.f = purrr::possibly(
      score_genes_for_expression_pattern, NA_real_),
      .x = cell_state,
      .y = gene_classes,
      state_graph,
      estimate_matrix))
  return(cell_states)
}

#' get the parent(s) of a state in a state transition graph
get_parents = function(state_graph, cell_state){
  parents = igraph::neighbors(state_graph, cell_state, mode="in")
  if (length(parents) > 0)
    return (parents$name)
  else
    return (c())
}

#' get the children of a state in a state transition graph
get_children = function(state_graph, cell_state){
  children = igraph::neighbors(state_graph, cell_state, mode="out")
  if (length(children) > 0)
    return (children$name)
  else
    return (c())
}

#' get the siblings of a state in a state transition graph
get_siblings = function(state_graph, cell_state){
  parents = get_parents(state_graph, cell_state)
  siblings = igraph::neighbors(state_graph, parents, mode="out")
  siblings = setdiff(siblings$name, cell_state) #exclude self
  return(siblings)
}

#' Compute a pseudobulk expression matrix for a model
pseudobulk_cds_for_states <- function(ccm, state_col=NULL, collapse_samples=FALSE){

  if (is.null(state_col)){
    cell_group_df = tibble::rownames_to_column(ccm@ccs@metadata[["cell_group_assignments"]])
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = cell_group)
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, "cell_group", sep="_")) %>% dplyr::select(rowname, pseudobulk_id, cell_group)
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, cell_group) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
    #%>% select(rowname, cell_group)
  }else{
    cell_group_df = tibble::rownames_to_column(ccm@ccs@metadata[["cell_group_assignments"]])
    cds_group_df = colData(ccm@ccs@cds) %>%
      as.data.frame %>% tibble::rownames_to_column() %>% dplyr::select(rowname, !!sym(state_col))
    cell_group_df = left_join(cell_group_df, cds_group_df, by=c("rowname"))
    if (collapse_samples)
      cell_group_df = cell_group_df %>% mutate(group_id = !!sym(state_col))
    cell_group_df = cell_group_df %>%
      dplyr::mutate(pseudobulk_id = paste(group_id, !!sym(state_col), sep="_")) %>% dplyr::select(rowname, pseudobulk_id, !!sym(state_col))
    agg_coldata = cell_group_df %>%
      dplyr::group_by(pseudobulk_id, !!sym(state_col)) %>%
      dplyr::summarize(num_cells_in_group = n()) %>%
      as.data.frame
  }

  agg_expr_mat = monocle3::aggregate_gene_expression(ccm@ccs@cds,
                                                     cell_group_df=cell_group_df,
                                                     norm_method="size_only",
                                                     scale_agg_values = FALSE,
                                                     pseudocount=0,
                                                     cell_agg_fun="mean")

  agg_expr_mat = agg_expr_mat[,agg_coldata$pseudobulk_id]

  row.names(agg_coldata) = agg_coldata$pseudobulk_id
  agg_coldata = agg_coldata[colnames(agg_expr_mat),]

  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(ccm@ccs@cds) %>% as.data.frame)
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  return(pseudobulk_cds)
}
