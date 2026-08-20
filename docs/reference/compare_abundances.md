Compare two estimates of cell abundances from a Hooke model.

    compare_abundances(
      ccm,
      cond_x,
      cond_y,
      by = "cell_group",
      method = c("BH", "bonferroni", "hochberg", "hommel", "BY"),
      alpha = 0.05,
      power = 0.8,
      convert_scale = FALSE,
      adjust_q_values = FALSE,
      log_scale = c("log", "log10", "log2")
    )

Arguments
---------

ccm

A cell\_count\_model.

cond\_x

tibble A cell type abundance estimate from estimate\_abundances().

cond\_y

tibble A cell type abundance estimate from estimate\_abundances().

by

string The column name used to join the two estimates.

method

string A method for correcting P-value multiple comparisons. This can be "BH" (Benjamini & Hochberg), "bonferroni" (Bonferroni), "hochberg" (Hochberg), "hommel", (Hommel), or "BYH" (Benjamini & Yekutieli).

alpha

Desired significance level

power

Desired power level for calculating minimum detectable fold change

convert\_scale

Whether to convert to log2 scale.

adjust\_q\_values

logical Whether to recompute `delta_q_value` using only rows whose power meets the `power` threshold.

log\_scale

Log scale used for estimate\_abundances.

Value
-----

tibble A table contrasting cond\_x and cond\_y (interpret as Y/X).
