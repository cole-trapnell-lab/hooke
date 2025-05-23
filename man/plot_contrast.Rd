% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_contrast}
\alias{plot_contrast}
\title{Plot a UMAP colored by how cells shift in a given contrast}
\usage{
plot_contrast(
  ccm,
  cond_b_vs_a_tbl,
  log_abundance_thresh = -5,
  scale_shifts_by = c("receiver", "sender", "none"),
  edge_size = 2,
  cell_size = 1,
  q_value_thresh = 1,
  group_label_size = 2,
  plot_labels = c("significant", "all", "none"),
  fc_limits = c(-3, 3),
  sender_cell_groups = NULL,
  receiver_cell_groups = NULL,
  plot_edges = c("none", "all", "directed", "undirected"),
  edge_significance = c("both", "one-sided"),
  keep_colors = FALSE,
  label_cell_groups = list(),
  repel_labels = TRUE,
  model_for_pcors = c("reduced", "full"),
  switch_label = NULL,
  sub_cds = NULL,
  alpha = 1,
  mask = list(),
  downsample = NULL,
  x = 1,
  y = 2
)
}
\arguments{
\item{ccm}{A cell_count_model object.}

\item{cond_b_vs_a_tbl}{data.frame A data frame from compare_abundances.}

\item{log_abundance_thresh}{numeric Select cell groups by log abundance.}

\item{scale_shifts_by}{string A scale directed graph edges by "sender",
"receiver", or "none".}

\item{edge_size}{numeric The size of edges in the plot.}

\item{cell_size}{numeric The size of cells in the plot.}

\item{q_value_thresh}{numeric Remove contrasts whose change in
q-value exceeds q_value_thresh.}

\item{group_label_size}{numeric The size of group labels in the plot.}

\item{plot_labels}{string Choose cell groups to label.}

\item{fc_limits}{vector The range of cell abundance changes to
include in the plot.}

\item{sender_cell_groups}{list Sender cell groups of directed graph.}

\item{receiver_cell_groups}{list Receiver cell groups of directed graph.}

\item{plot_edges}{string Type of edges to plot.}

\item{label_cell_groups}{list The cell_group labels to include in the plot.}

\item{repel_labels}{logical Repel overlapping plot labels.}

\item{model_for_pcors}{string The model to use for orienting graph
edges from the PLNnetwork.}

\item{switch_label}{string The name of the cell_data_set column with cell_group identifiers.}

\item{sub_cds}{string A cell_data_set.}

\item{alpha}{numeric A the ggplot opacity. A value between 0 and 1.}

\item{downsample}{how much to downsample the plots}

\item{x}{numeric The column number for the UMAP x coordinate.}

\item{y}{numeric The column number for the UMAP y coordinate.}
}
\value{
A ggplot2 plot object.
}
\description{
Plot a UMAP colored by how cells shift in a given contrast
}
