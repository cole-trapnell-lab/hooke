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
