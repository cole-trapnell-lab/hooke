Helper function to plot kinetics — compare\_ko\_to\_wt\_at\_timepoint • hooke

Toggle navigation [hooke](../index.html) 0.0.1

*   [Reference](../reference/index.html)
*   [Articles](#)
    *   [Hooke Tutorial](../articles/hooke_tutorial.html)

Helper function to plot kinetics
================================

`compare_ko_to_wt_at_timepoint.Rd`

Helper function to plot kinetics

    compare_ko_to_wt_at_timepoint(
      tp,
      perturbation_ccm,
      wt_pred_df,
      ko_pred_df,
      interval_col
    )

Arguments
---------

tp

timepoint

perturbation\_ccm

a cell count model with a perturbation

wt\_pred\_df

control output from estimate\_abundances\_over\_interval()

ko\_pred\_df

perturbation output from estimate\_abundances\_over\_interval()

interval\_col

column that matches the timepoint information

Contents
--------

Developed by Cole Trapnell.

Site built with [pkgdown](https://pkgdown.r-lib.org/) 2.0.7.