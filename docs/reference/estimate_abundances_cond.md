Predict cell type abundances given a PLN model and a set of inputs for its covariates and observed counts

    estimate_abundances_cond(
      ccm,
      newdata,
      cond_responses,
      min_log_abund = -5,
      cell_group = "cell_group",
      type = c("link", "response"),
      pln_model = c("full", "reduced")
    )

Arguments
---------

ccm

A cell\_count\_model, fit with `covariance_type = "full"`.

newdata

tibble A tibble of variables used for the prediction. Must either be a single row or a tibble with one row per sample of the cell count set for ccm.

cond\_responses

A data frame containing the observed counts of the conditioning cell type(s); column names must match cell types in the model.

min\_log\_abund

numeric Minimum log abundance value.

cell\_group

string The name of the groups that are being estimated.

type

string Scale of the predictions. Either "link" (default, latent scale) or "response" (mean response scale).

pln\_model

string Which PLN model form to use for prediction. Either "full" (default) or "reduced".

Value
-----

A tibble of cell abundance predictions for the unconditioned cell types.
