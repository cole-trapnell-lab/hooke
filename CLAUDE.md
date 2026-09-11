# hooke

R package: differential analysis of **cell counts / abundances** in single-cell experiments.
Upstream of `platt`; consumed by sulston, mcclintock, lewis, zscape_portal.

## Commands (verified in Makefile)
```bash
make fast     # == make test
make test     # Rscript -e 'testthat::test_local(".", reporter="summary", stop_on_failure=TRUE)'
make check    # R CMD check --no-manual .
```
3 test files: `test-cell_count_model.R`, `test-cell_count_set.R`, `test-contrasts.R`.

## Gotchas
- `HOOKE_SKIP_DEP_CHECK=1` **skips the entire test suite**, not just a dependency probe.
  A green `make fast` under that flag means nothing ran. (The root README describes it
  as bypassing a "dependency check" — that wording is wrong.)
- `Depends: PLNmodels`, which is not on CRAN in a usable version — CI installs
  `PLN-team/PLNmodels@master`. Locally: `Rscript scripts/install_hooke_deps.R` from the stack root.
- `Remotes: bioc::Rgraphviz`. `LinkingTo: Rcpp` + `src/` — changing `src/` needs a recompile.
- 0.0.3 **removed the `mdfc` column**; `calculate_mdfc()` was fixed in place and now emits
  `mdfc80`. The two are not interchangeable — the value differs on essentially every row,
  so thresholds picked against `mdfc` need rechoosing. `calculate_power()` is now
  `calculate_observed_power()`, published as `observed_power`, with `power` a deprecated
  alias for one release. `observed_power` is `delta_p_value` restated (strictly increasing
  in `|Z|`), so filtering on it before `p.adjust()` is anti-conservative. Filter on
  `mdfc80` or `power_at_margin` — the only effect-independent columns.

## Release notes
`NEWS.md` exists (seeded at 0.0.1, no prior history). Bump `Version:` in `DESCRIPTION` and add the
matching `# hooke <version>` / `### Changes` section in the same commit. See the root `CLAUDE.md`
convention and `repos/monocle3/NEWS.md` for the reference format.

## Architecture
`R/` is the whole surface; exported API in `NAMESPACE`, roxygen markdown enabled (7.3.2).

## CI
`.github/workflows/check_on_push.yml`, `on: [push]`: `R CMD build` + `R CMD check` inside
`ghcr.io/cole-trapnell-lab/monocle3_depend:v1.3.0_1`, with `TRAVIS=true` and
`_R_CHECK_TESTS_NLINES_=0`.
