Predict cell type abundances given a PLN model and a set of inputs for its covariates — estimate\_abundances • hooke

Toggle navigation [hooke](../index.html) 0.0.1

*   [Reference](../reference/index.html)
*   [Articles](#)
    *   [Hooke Tutorial](../articles/hooke_tutorial.html)

Predict cell type abundances given a PLN model and a set of inputs for its covariates
=====================================================================================

`estimate_abundances.Rd`

Predict cell type abundances given a PLN model and a set of inputs for its covariates

    estimate_abundances(ccm, newdata, min_log_abund = -5)

Arguments
---------

ccm

A cell\_count\_model.

newdata

tibble A tibble of variables used for the prediction.

min\_log\_abund

numeric Minimum log abundance value.

Value
-----

A tibble of cell abundance predictions.

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.