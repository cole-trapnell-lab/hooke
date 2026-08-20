Predict cell type abundances given a PLN model over a range of time or other interval

    estimate_abundances_over_interval(
      ccm,
      interval_start,
      interval_stop,
      interval_col = "timepoint",
      interval_step = 2,
      min_log_abund = -5,
      newdata = tibble()
    )

Arguments
---------

ccm

A cell\_count\_model.

interval\_start

numeric Interval start value.

interval\_stop

numeric Interval stop value.

interval\_col

character Interval values are taken from the interval\_var data. Default is "timepoint".

interval\_step

numeric Interval size. Default is 2.

min\_log\_abund

numeric Minimum log abundance value.

newdata

tibble Additional covariate values to combine with the interval sequence (cross-joined with each interval value). If it contains the `interval_col` column, that column is dropped before joining. Defaults to an empty tibble.

Value
-----

A tibble of cell abundance predictions.
