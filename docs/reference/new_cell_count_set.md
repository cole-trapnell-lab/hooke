Create a new cell\_data\_set object.

    new_cell_count_set(
      cds,
      sample_group,
      cell_group,
      sample_metadata = NULL,
      cell_metadata = NULL,
      lower_threshold = NULL,
      upper_threshold = NULL,
      keep_cds = TRUE,
      norm_method = c("size_factors", "TSS", "CSS", "RLE", "GMPR", "Wrench", "none"),
      size_factors = NULL,
      pseudocount = 0
    )

Arguments
---------

cds

A Monocle cell data set object.

sample\_group

A column in colData(cds) that specifes how cells are grouped into samples.

cell\_group

A column in colData(cds) that specifies how cells are grouped into types or states (e.g. cluster).

sample\_metadata

Data frame containing attributes of individual samples, where the column named 'sample' has entries in `sample_group`.

cell\_metadata

Data frame containing attributes of individual cell groups, where `row.names(cell_metadata)` are entries in `cell_group`

lower\_threshold

numeric Minimum number of cells in retained cell\_groups.

upper\_threshold

numeric Maximum number of cells in retained cell\_groups.

norm\_method

string Normalization method used to compute scaling factors used as offset during PLN inference.

size\_factors

numeric vector or matrix User supplied vector or matrix of offsets passed the PLNmodels::prepare\_data() method.

Value
-----

a new cell\_data\_set object

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.