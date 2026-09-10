# Deprecated functions moved from utils.R.

#' @noRd
#' @note Deprecated.
collect_genotype_effects <- function(ccm, timepoint = 24, expt = "GAP16") {
  .Deprecated("collect_genotype_effects", msg = "collect_genotype_effects is deprecated.")
  control_abund <- estimate_abundances(ccm, tibble(knockout = FALSE, timepoint = timepoint, expt = expt))
  knockout_abund <- estimate_abundances(ccm, tibble(knockout = TRUE, timepoint = timepoint, expt = expt))
  genotype_comparison_tbl <- compare_abundances(ccm, control_abund, knockout_abund)
}
