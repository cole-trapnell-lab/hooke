Plot a UMAP colored by how cells shift in a given contrast — plot\_contrast • hooke

Toggle navigation [hooke](../index.html) 0.0.1

*   [Reference](../reference/index.html)
*   [Articles](#)
    *   [Hooke Tutorial](../articles/hooke_tutorial.html)

Plot a UMAP colored by how cells shift in a given contrast
==========================================================

`plot_contrast.Rd`

Plot a UMAP colored by how cells shift in a given contrast

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
      x = 1,
      y = 2
    )

Arguments
---------

ccm

A cell\_count\_model object.

cond\_b\_vs\_a\_tbl

data.frame A data frame from compare\_abundances.

log\_abundance\_thresh

numeric Select cell groups by log abundance.

scale\_shifts\_by

string A scale directed graph edges by "sender", "receiver", or "none".

edge\_size

numeric The size of edges in the plot.

cell\_size

numeric The size of cells in the plot.

q\_value\_thresh

numeric Remove contrasts whose change in q-value exceeds q\_value\_thresh.

group\_label\_size

numeric The size of group labels in the plot.

plot\_labels

string Choose cell groups to label.

fc\_limits

vector The range of cell abundance changes to include in the plot.

sender\_cell\_groups

list Sender cell groups of directed graph.

receiver\_cell\_groups

list Receiver cell groups of directed graph.

plot\_edges

string Type of edges to plot.

label\_cell\_groups

list The cell\_group labels to include in the plot.

repel\_labels

logical Repel overlapping plot labels.

model\_for\_pcors

string The model to use for orienting graph edges from the PLNnetwork.

switch\_label

string The name of the cell\_data\_set column with cell\_group identifiers.

sub\_cds

string A cell\_data\_set.

alpha

numeric A the ggplot opacity. A value between 0 and 1.

x

numeric The column number for the UMAP x coordinate.

y

numeric The column number for the UMAP y coordinate.

Value
-----

A ggplot2 plot object.

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.