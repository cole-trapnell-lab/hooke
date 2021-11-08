setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))
setOldClass(c("PLNnetworkfamily"), prototype=structure(list(), class="PLNnetworkfamily"))
setOldClass(c("PLNnetworkfit"), prototype=structure(list(), class="PLNnetworkfit"))


#' The cell_count_set class
#'
#' The main class used by Hooke to hold cell abundances data.
#' cell_count_set extends the Monocle's cell_data_set class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field cds cell_data_set, the Monocle cell data set object that this class models.
#' @name cell_count_set
#' @rdname cell_count_set
#' @aliases cell_count_set-class
#' @exportClass cell_count_set
setClass("cell_count_set",
         contains = "cell_data_set",
         slots = c(cds = "cell_data_set")
)


#' The cell_count_model class
#'
#' The main class used by Monocle3 to hold single-cell expression data.
#' cell_count_model extends the Bioconductor SingleCellExperiment class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field ccs cell_count_set the underlying object of cell counts this class models.
#' @field model_formula_str character string, specifies the model for the cell abundances.
#' @field model_family PLNnetworkfamily a family of PLNnetwork models for the cell count data.
#' @field best_model PLNnetworkfit the current best PLN network model for the cell count data.
#' @field model_aux SimpleList auxiliary information from from PLN model construction
#' @name cell_count_model
#' @rdname cell_count_model
#' @aliases cell_count_model-class
#' @import PLNmodels
#' @exportClass cell_count_model
setClass("cell_count_model",
         slots = c(ccs = "cell_count_set",
                   model_formula_str = "character",
                   # FIXME: for some reason I can't include these as slots. Maybe because they are R6? dunno. For now, they live in model_aux
                   model_family = "PLNnetworkfamily",
                   best_model = "PLNnetworkfit",
                   model_aux = "SimpleList")
)



#' Create a new cell_data_set object.
#'
#' @param cds A Monocle cell data set object.
#' @param sample_group A column in colData(cds) that specifes how cells are grouped into samples
#' @param cell_group A column in colData(cds) that specifies how cells are grouped into types or states (e.g. cluster)
#' @param sample_metadata data frame containing attributes of individual samples, where
#' @param cell_metadata data frame containing attributes of individual cell groups, where
#'   \code{row.names(cell_metadata)} are entries in \code{cell_group}
#' @return a new cell_data_set object
#' @importFrom dplyr summarize
#' @importFrom dplyr %>%
#' @export
new_cell_count_set <- function(cds,
                               sample_group,
                               cell_group,
                               sample_metadata = NULL,
                               cell_metadata = NULL) {

  coldata_df = colData(cds) %>% tibble::as_tibble()
  coldata_df$cluster = clusters(cds)
  coldata_df$partition = partitions(cds)

  coldata_df = coldata_df %>% dplyr::rename_("sample" = sample_group, "cell_group" = cell_group)

  coldata_df$group_id = coldata_df %>%
    dplyr::group_indices_("sample", "cell_group") %>% as.character

  cds_summary = coldata_df %>%
    dplyr::group_by_("sample", "cell_group") %>%
    dplyr::summarize(cells = dplyr::n())

  cds_covariates_df = coldata_df %>%
    dplyr::group_by_("sample") %>%
    dplyr::summarize(across(where(is.numeric), mean),
                     across(where(is.factor), function(x) { tail(names(sort(table(x))), 1) }),
                     across(where(is.character), function(x) { tail(names(sort(table(x))), 1) } ))

  if (is.null(sample_metadata) == FALSE){
    cds_covariates_df = left_join(cds_covariates_df, sample_metadata, by=c("sample"="sample"))
  }

  cds_covariates_df = cds_covariates_df %>% as.data.frame(cds_covariates_df)
  row.names(cds_covariates_df) = cds_covariates_df %>% dplyr::pull(sample)

  cell_counts_wide = tidyr::spread(cds_summary, sample, cells, fill=0)
  cell_states = as.character(cell_counts_wide %>% dplyr::pull("cell_group"))
  cell_counts_wide = as.matrix(cell_counts_wide[,3:ncol(cell_counts_wide)])
  row.names(cell_counts_wide) = cell_states
  #cell_counts_wide = t(cell_counts_wide)

  cds_covariates_df = cds_covariates_df[colnames(cell_counts_wide),]

  # This is super confusing because of the way the arguments are named in new_cell_data_set.
  # We are making a matrix of dimension MxN, where M are cell types and N are samples (e.g. embryos, replicates, etc).
  # The "gene" metadata monocle normally expects will actually be used to hold
  ccs = methods::new("cell_count_set",
               monocle3::new_cell_data_set(cell_counts_wide,
                                           cell_metadata=cds_covariates_df,
                                           gene_metadata=cell_metadata),
               cds=cds)


  # assertthat::assert_that(class(expression_data) == "matrix" ||
  #                           is_sparse_matrix(expression_data),
  #                         msg = paste("Argument expression_data must be a",
  #                                     "matrix - either sparse from the",
  #                                     "Matrix package or dense"))
  # if (!is.null(cell_metadata)) {
  #   assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data),
  #                           msg = paste("cell_metadata must be NULL or have",
  #                                       "the same number of rows as columns",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(cell_metadata)) &
  #                             all(row.names(cell_metadata) == colnames(expression_data)),
  #                           msg = paste("row.names of cell_metadata must be equal to colnames of",
  #                                       "expression_data"))
  # }
  #
  # if (!is.null(gene_metadata)) {
  #   assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data),
  #                           msg = paste("gene_metadata must be NULL or have",
  #                                       "the same number of rows as rows",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(
  #     row.names(gene_metadata) == row.names(expression_data)),
  #     msg = paste("row.names of gene_metadata must be equal to row.names of",
  #                 "expression_data"))
  # }
  #
  # if (is.null(cell_metadata)) {
  #   cell_metadata <- data.frame(cell = colnames(expression_data),
  #                               row.names = colnames(expression_data))
  # }
  #
  # if(!('gene_short_name' %in% colnames(gene_metadata))) {
  #   warning(paste("Warning: gene_metadata must contain a column verbatim",
  #                 "named 'gene_short_name' for certain functions."))
  # }
  #
  # sce <- SingleCellExperiment(list(counts=methods::as(expression_data, "dgCMatrix")),
  #                             rowData = gene_metadata,
  #                             colData = cell_metadata)
  #
  # cds <- methods::new("cell_data_set",
  #                     assays = SummarizedExperiment::Assays(
  #                       list(counts=methods::as(expression_data, "dgCMatrix"))),
  #                     colData = colData(sce),
  #                     int_elementMetadata =int_elementMetadata(sce),
  #                     int_colData = int_colData(sce),
  #                     int_metadata = int_metadata(sce),
  #                     metadata = metadata(sce),
  #                     NAMES = NULL,
  #                     elementMetadata = elementMetadata(sce)[,0],
  #                     rowRanges = rowRanges(sce))
  #
  # metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  # clusters <- stats::setNames(SimpleList(), character(0))
  # cds <- estimate_size_factors(cds)
  # cds

  ccs@metadata[["cell_group_assignments"]] = coldata_df %>% dplyr::select(group_id, sample, cell_group) %>% as.data.frame
  row.names(ccs@metadata[["cell_group_assignments"]]) = colnames(cds)

  return (ccs)
}


#' Create a new cell_count_model object.
#'
#' Fits a PLNnetwork according to a formula. Accepts a matrix of penalties as a
#' way of encoding a graph prior. Automatically selects sparsity parameter, but
#' allows user to update it.
#'
#' @param ccs A Hooke cell_count_set object.
#' @param model_formula_str A character string specifying the model of cell abundances across samples.
#' where terms refer to columns in\code{colData(ccs)}
#' @param penalty_matrix A numeric NxN symmetric matrix specifying penalties for
#'   the PLN model, where N is the number of cell types. Entries must be
#'   positive. Use to specify an undirected graph prior for the PLN model.
#' @return a new cell_count_model object
#' @importFrom PLNmodels prepare_data
#' @importFrom PLNmodels PLNnetwork
#' @importFrom PLNmodels getBestModel
#' @importFrom PLNmodels getModel
#' @export
new_cell_count_model <- function(ccs, model_formula_str, penalty_matrix = NULL, whitelist=NULL, blacklist=NULL, ...) {


  pln_data <- PLNmodels::prepare_data(counts = counts(ccs),
                                      covariates = colData(ccs) %>% as.data.frame,
                                      offset = "none")#size_factors(ccs))

  model_formula_str = paste("Abundance", model_formula_str)
  model_formula = as.formula(model_formula_str)
  #pln_data <- as.name(deparse(substitute(pln_data)))

  initial_penalties = init_penalty_matrix(ccs, whitelist=whitelist, blacklist=blacklist, ...)
  initial_penalties = initial_penalties[colnames(pln_data$Abundance), colnames(pln_data$Abundance)]
  # INSANE R BULLSHIT ALERT: for reasons I do not understand,
  # calling the fit via do.call prevents a weird error with formula
  # created with as.formula (e.g. after pasting).
  pln_model <- do.call(PLNmodels::PLNnetwork, args=list(model_formula,
                                     data=pln_data,
                                     control_init=list(min.ratio=0.001),
                                     control_main=list(penalty_weights=initial_penalties)))

  # Choose a model that isn't very aggressively sparsified
  best_model <- PLNmodels::getBestModel(pln_model, "EBIC")
  best_model <- PLNmodels::getModel(pln_model, var=best_model$penalty/10)

  # assertthat::assert_that(class(expression_data) == "matrix" ||
  #                           is_sparse_matrix(expression_data),
  #                         msg = paste("Argument expression_data must be a",
  #                                     "matrix - either sparse from the",
  #                                     "Matrix package or dense"))
  # if (!is.null(cell_metadata)) {
  #   assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data),
  #                           msg = paste("cell_metadata must be NULL or have",
  #                                       "the same number of rows as columns",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(cell_metadata)) &
  #                             all(row.names(cell_metadata) == colnames(expression_data)),
  #                           msg = paste("row.names of cell_metadata must be equal to colnames of",
  #                                       "expression_data"))
  # }
  #
  # if (!is.null(gene_metadata)) {
  #   assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data),
  #                           msg = paste("gene_metadata must be NULL or have",
  #                                       "the same number of rows as rows",
  #                                       "in expression_data"))
  #   assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(
  #     row.names(gene_metadata) == row.names(expression_data)),
  #     msg = paste("row.names of gene_metadata must be equal to row.names of",
  #                 "expression_data"))
  # }
  #
  # if (is.null(cell_metadata)) {
  #   cell_metadata <- data.frame(cell = colnames(expression_data),
  #                               row.names = colnames(expression_data))
  # }
  #
  # if(!('gene_short_name' %in% colnames(gene_metadata))) {
  #   warning(paste("Warning: gene_metadata must contain a column verbatim",
  #                 "named 'gene_short_name' for certain functions."))
  # }
  #
  ccm <- methods::new("cell_count_model",
                      ccs = ccs,
                      model_formula_str = model_formula_str,
                      best_model = best_model,
                      model_family = pln_model,
                      model_aux = SimpleList()
                      )
  #
  # metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  # clusters <- stats::setNames(SimpleList(), character(0))
  # cds <- estimate_size_factors(cds)
  # cds
  #ccm@model_aux[["best_model"]] = best_model
  #ccm@model_aux[["model_family"]] = pln_model


  ccm
}

#' Select the model a cell_count_model should use
#'
#' @param ccm A cell_count_model object.
#' @param criterion a character string specifying the PLNmodels criterion to use. Must be one of "BIC", "EBIC" or "StARS".
#' @return an updated cell_count_model object
#' @importFrom PLNmodels getBestModel
#' @importFrom PLNmodels getModel
#' @export
select_model <- function(ccm, criterion = "EBIC", sparsity_factor=1.0)
{
  best_model <- PLNmodels::getBestModel(pln_model, criterion)
  best_model <- PLNmodels::getModel(pln_model, var=best_model$penalty * sparsity_factor)
  ccm@best_model = best_model
  return(ccm)
}


#' Initialize the PLN network penalty matrix, accepting optional whitelists and
#' blacklists of edges that are "free" or "off limits" between cell groups
#' @param ccs A cell_count_set of aggregated cell counts
#' @param whitelist a data frame with two columns corresponding to (undirected) edges that should receive no penalty
#' @param blacklist a data frame with two columns corresponding to (undirected) edges that should receive very high penalty
#' @param dist_fun A function that returns a penalty based given a distance between two clusters
init_penalty_matrix = function(ccs, whitelist=NULL, blacklist=NULL, base_penalty = 1, min_penalty=0.01, max_penalty=1e6){
  cell_group_centroids = centroids(ccs)
  dist_matrix = as.matrix(dist(cell_group_centroids, method = "euclidean", upper=T, diag = T))

  row.names(dist_matrix) <- row.names(cell_group_centroids)
  colnames(dist_matrix) <- row.names(cell_group_centroids)

  # TODO: do I need this? Probably the caller can and should do this.
  #dist_matrix = dist_matrix[colnames(data$Abundance), colnames(data$Abundance)]

  get_rho_mat <- function(DM, distance_parameter = 1, s=1, xmin = NULL) {
    if (is.null(xmin)){
      xmin = min(DM[DM > 0]) / 2
    }
    #out <- (1-(xmin/DM)^s) * distance_parameter
    out = min_penalty + (DM / max(DM))^2
    #out =  1 + DM^s
    # penalties have to be > 0
    out[!is.finite(out)] <- min_penalty
    out[out < 0] <- min_penalty
    return(out)
  }

  penalty_matrix = base_penalty * (1 + get_rho_mat(dist_matrix, distance_parameter=1, s=2))

  # FIXME: This might only actually work when grouping cells by clusters and cluster names are
  # integers. We should make sure this generalizes when making white/black lists of cell groups
  # by type or other groupings
  if (is.null(whitelist) == FALSE){
    penalty_matrix[as.matrix(whitelist[,c(1,2)])] = min_penalty
    penalty_matrix[as.matrix(whitelist[,c(2,1)])] = min_penalty

  }

  if (is.null(blacklist) == FALSE){
    penalty_matrix[as.matrix(blacklist[,c(1,2)])] = max_penalty
    penalty_matrix[as.matrix(blacklist[,c(2,1)])] = max_penalty
  }

  # TODO: add support for whitelisting and blacklistint
  #qplot(as.numeric(dist_matrix), as.numeric(out))
  return(penalty_matrix)
}

