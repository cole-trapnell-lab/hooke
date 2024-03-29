% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.R
\name{pseudobulk_ccs_for_states}
\alias{pseudobulk_ccs_for_states}
\title{Compute a pseudobulk expression matrix for a ccs}
\usage{
pseudobulk_ccs_for_states(ccs, state_col = NULL, collapse_samples = FALSE)
}
\arguments{
\item{ccs}{a cell count set object}

\item{state_col}{column to aggregate expression. Defaults using the cell_group used in ccs construction.}

\item{collapse_samples}{boolean Whether to collapse sample groups into one.}
}
\description{
Compute a pseudobulk expression matrix for a ccs
}
