# hooke 0.0.3

### Breaking changes

* `compare_abundances()` no longer returns `mdfc`. `calculate_mdfc()` has been
  fixed rather than frozen, and the corrected quantity is published as
  `mdfc80`. The two are **not interchangeable**: the value differs on
  essentially every row, and any threshold chosen against `mdfc` needs
  rechoosing. Code reading `mdfc` now fails loudly rather than silently getting
  a number that means something else.
  * The old value was `exp((z_alpha + z_observed) * SE)`, because a
    `dplyr::mutate()` bound `power = power` to the column created on the line
    above instead of to the formal argument, feeding the observed power in
    where a requested power level belonged. It was understated on most rows,
    smallest exactly where the data are weakest, and `Inf` on the *most*
    significant rows (`qnorm(1) = Inf`) rather than on degenerate ones.
  * `calculate_mdfc()` now takes `df` and `base`, and its `Inf` guard tests the
    standard errors alone. The old guard also tested `power == 0`, so it fired
    only as a side effect of the contamination; with a constant power it would
    have reported a 1-fold detection limit for a degenerate fit. Such rows are
    now `NA`, with the reason in `contrast_note`.
* `calculate_power()` is renamed `calculate_observed_power()`, and the column is
  published as `observed_power`. `power` remains as a **deprecated alias** for
  one release and is bit-identical to 0.0.2.

### Changes

* `compare_abundances()` gains the columns needed to tell "no phenotype" apart
  from "not powered". Every retained column, `power` included, is bit-identical
  to 0.0.2.
  * `mdfc80` -- minimum detectable fold change at the requested `power`: the
    smallest change the contrast could have caught.
  * `power_at_margin` -- power to detect a change of `margin`. The dual of
    `mdfc80`: it fixes the effect size and reports the power, where `mdfc80`
    fixes the power and reports the effect size. Computed exactly, from the
    noncentral t.
  * `margin_fold_change` -- the `margin` used, so the table records what it was
    powered for rather than leaving a reader to guess.
  * `delta_log_abund_lo` / `delta_log_abund_hi` -- Wald confidence interval on
    `delta_log_abund` at level `1 - alpha`.
  * `df_resid` -- residual degrees of freedom. Previously computed and
    discarded, so nothing downstream could rebuild an interval or an
    equivalence test from the published table.
  * `contrast_note` -- why `mdfc80` and `power_at_margin` are `NA`:
    `"degenerate_fit"` or `"insufficient_df"`; `NA` when the row is fine.
* New `margin` argument, the effect size to power against, given as a **fold
  change** (default `2`) so it means the same thing regardless of `log_scale`
  and `convert_scale`. It is converted to the result's log scale internally.
* `mdfc80` and `power_at_margin` are **effect-independent**: functions of
  `delta_log_abund_se`, `alpha`, `df_resid` and a declared constant only, never
  of the observed effect. That is what makes either a fair answer to "were we
  powered to see a phenotype here?", and what makes them safe to filter on.
  Both use t quantiles, matching `delta_p_value`, rather than the normal
  approximation; at 8 residual df the normal approximation understates the
  detectable effect by roughly 15%. Both exponentiate in the base implied by
  `log_scale`/`convert_scale`, where the old `calculate_mdfc()` hardcoded
  `exp()`.
* `observed_power` is retained but should not be used to decide whether a null
  result is meaningful. It substitutes the observed Wald statistic for the true
  effect, which makes it strictly increasing in `|Z|` and so a deterministic
  restatement of `delta_p_value`: `>= 0.8` is exactly `p <= ~0.008`. It does not
  measure precision -- a cell type with a huge standard error and a fluke
  estimate scores high, while a tightly measured genuine null scores the floor.
* `compare_abundances(adjust_q_values = TRUE)` now restricts the correction
  using `power_at_margin` rather than `power`, and warns. Filtering on `power`
  selected rows by their p-value and then adjusted those same p-values, which is
  anti-conservative; `power_at_margin` is effect-independent, so the filter is
  defensible. It still changes what the FDR guarantee covers -- the adequately
  powered subset, not every row tested. The default is unchanged (`FALSE`), and
  it is not enabled in production.

# hooke 0.0.2

### Changes

* Fix per-sample aggregation of logical covariates in `new_cell_count_set()`. A
  logical `colData` column was collapsed with `sum(x) == 1`, so a flag that is
  constant within a sample (for example `knockout`) became `FALSE` for every
  sample and dropped out of downstream models. Logical columns now take a
  majority vote, matching how factor and character columns are collapsed.
  Models fit on data with logical covariates will give different results.
* `subset_ccs()` has a help page again. Its roxygen block was malformed, so the
  exported function had no documentation.

# hooke 0.0.1

Release notes start here. Hooke has not yet had a version bump, so there is no
prior history to document.

When you bump `Version:` in `DESCRIPTION`, add a section above this one:

```
# hooke <new version>

### Changes

* One bullet per user-visible change.
```

Document what a caller would notice — new or removed exported functions, changed
defaults, changed return shapes, bug fixes that alter results. Internal
refactors that leave the API and the numbers alone do not need an entry.
