# Conditional prediction

A `cell_count_model` doesn't just describe how each cell type's abundance depends on your covariates independently — because it's fit as a full multivariate PLN model, it also captures how cell types' abundances covary with each other. [`estimate_abundances()`](https://cole-trapnell-lab.github.io/hooke/about/#estimating-abundances) uses that model to answer "what abundance does the model predict for each cell type, given these covariate values?" [`estimate_abundances_cond()`](https://cole-trapnell-lab.github.io/hooke/reference/estimate_abundances_cond) answers a different question: "given that I've *observed* these counts for some cell types in this sample, what does the model now predict for the others?"

That second question only has an answer if the model actually has a notion of cell types covarying in the first place. Fit `new_cell_count_model()` with the default `covariance_type` (`"spherical"`) and there's nothing to condition on — every cell type's variance is independent. You need `covariance_type = "full"`:

```
ccm_full = new_cell_count_model(ccs,
                                main_model_formula_str = "~ exposed",
                                covariance_type = "full")
```

_See the [`new_cell_count_model()` reference page](https://cole-trapnell-lab.github.io/hooke/reference/new_cell_count_model) for the other covariance options — `"full"` is more expensive to fit than `"spherical"`/`"diagonal"`, which is why it isn't the default, but it's required for anything in this page._

## Predicting one cell type's abundance from another's observed count

[`estimate_abundances_cond()`](https://cole-trapnell-lab.github.io/hooke/reference/estimate_abundances_cond) takes the same kind of `newdata` as `estimate_abundances()`, plus a `cond_responses` table of the observed counts you want to condition on:

* `ccm` - a `cell_count_model` fit with `covariance_type = "full"`
* `newdata` - a tibble of covariate values, same as for `estimate_abundances()`
* `cond_responses` - a data frame of observed counts for the cell type(s) you're conditioning on; column names must match cell types in the model
* `min_log_abund` - floor for predicted log-abundance, default `-5`
* `cell_group` - the name of the groups being estimated
* `type` - `"link"` (latent scale, default) or `"response"` (mean count scale)
* `pln_model` - `"full"` (default) or `"reduced"`, which PLN model form to predict from

Suppose we've observed an unusually high Aerocyte count in a new exposed sample, and want to know what that implies for Endothelial Capillary abundance — a cell type Aerocytes are known to associate with:

```
observed_counts = data.frame(Aerocytes = 40)

capillary_pred_cond = estimate_abundances_cond(ccm_full,
                                               newdata = tibble::tibble(exposed = "exposed"),
                                               cond_responses = observed_counts,
                                               type = "response")

capillary_pred_cond %>% head()
```

| exposed  | cell_group            | log_abund |
|----------|------------------------|-----------|
| exposed  | Endothelial Capillary  | 5.62      |
| exposed  | Endothelial Arterial   | 3.10      |
| exposed  | Endothelial Venous     | 3.24      |
| exposed  | Fibroblasts            | 5.31      |

Every cell type other than the one(s) in `cond_responses` gets a row — not just the one you're interested in — so filter down to the cell type(s) you actually care about.

## Conditional vs. unconditional predictions

Comparing this against the unconditional prediction for the same covariates — `estimate_abundances(ccm_full, tibble::tibble(exposed = "exposed"))` — shows how much (if at all) knowing the Aerocyte count actually moved the model's expectation for Endothelial Capillary abundance. If the two predictions are nearly identical, the model isn't finding much covariance between those two cell types once `exposed` is accounted for; a large, consistent gap across many samples is evidence of a real coupling between them.

_Cross-validating that gap across genotypes or conditions — fitting on one condition, conditioning on observed counts from another, and scoring the prediction error — is a real, useful pattern built on top of [`estimate_abundances_cond()`](https://cole-trapnell-lab.github.io/hooke/reference/estimate_abundances_cond), but it isn't packaged as a dedicated helper yet. You'd assemble it yourself the way the example above does, just repeated across held-out samples._
