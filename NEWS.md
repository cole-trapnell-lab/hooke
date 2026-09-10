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
