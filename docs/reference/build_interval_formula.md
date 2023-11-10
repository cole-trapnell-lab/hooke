Builds a model formula for time series models based on the range of the data. This is a utility function that puts the knots in reasonable positions based on the range of the data.

    build_interval_formula(
      ccs,
      num_breaks,
      interval_var = "timepoint",
      interval_start = NULL,
      interval_stop = NULL
    )

Arguments
---------

character

interval\_var

numeric

interval\_stop Interval stop value.

Value
-----

An interval model formula.
