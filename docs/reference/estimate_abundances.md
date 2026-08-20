Predict cell type abundances given a PLN model and a set of inputs for its covariates

    estimate_abundances(
      ccm,
      newdata,
      min_log_abund = -5,
      cell_group = "cell_group",
      scale = c("log", "log10", "log2", "per_1000")
    )

Arguments
---------

ccm

A cell\_count\_model.

newdata

tibble A tibble of variables used for the prediction.

min\_log\_abund

numeric Minimum log abundance value.

cell\_group

string The name of the groups that are being estimated.

scale

Desired log scale for the output. Default is natural log.

Value
-----

A tibble of cell abundance predictions.
