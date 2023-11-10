Compare two estimates of cell abundances from a Hooke model.

    compare_abundances(
      ccm,
      cond_x,
      cond_y,
      method = c("BH", "bonferroni", "hochberg", "hommel", "BY")
    )

Arguments
---------

ccm

A cell\_count\_model.

cond\_x

tibble A cell type abundance estimate from estimate\_abundances().

cond\_y

tibble A cell type abundance estimate from estimate from estimate\_abundances().

method

string A method for correcting P-value multiple comparisons. This can be "BH" (Benjamini & Hochberg), "bonferroni" (Bonferroni), "hochberg" (Hochberg), "hommel", (Hommel), or "BYH" (Benjamini & Yekutieli).

Value
-----

tibble A table contrasting cond\_x and cond\_y (interpret as Y/X).

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.