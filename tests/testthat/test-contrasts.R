context('test-contrasts')


test_that('estimate_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('estimate_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('estimate_abundances problems', {

  # On Github Actions
  # skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances problems', {

  # Not on Github Actions
  # skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})



# ---------------------------------------------------------------------------
# power / mdfc replacement: mdfc80 and power_at_margin
# ---------------------------------------------------------------------------

test_that('calculate_mdfc no longer depends on the observed effect', {

  se_x <- c(0.1, 0.4, 1.2)
  se_y <- c(0.2, 0.4, 0.3)

  mdfc80 <- calculate_mdfc(se_x, se_y, alpha = 0.05, power = 0.8, df = 40)

  # A strictly increasing function of the contrast SE and nothing else. The
  # column it replaces was essentially uncorrelated with precision.
  se <- sqrt(se_x^2 + se_y^2)
  expect_equal(order(mdfc80), order(se))
  expect_true(all(diff(mdfc80[order(se)]) > 0))
  expect_equal(cor(mdfc80, se, method = "spearman"), 1)

  # Closed form, with t quantiles to match delta_p_value.
  expect_equal(mdfc80, exp((qt(0.975, 40) + qt(0.8, 40)) * se))

  # Requesting more power demands a larger detectable effect.
  expect_true(all(
    calculate_mdfc(se_x, se_y, power = 0.9, df = 40) > mdfc80
  ))
})


test_that('calculate_mdfc is base aware', {

  # SEs on a log2 scale must exponentiate in base 2, not base e. The removed
  # calculate_mdfc() hardcoded exp(), inflating log2 results.
  se <- sqrt(0.5^2 + 0.5^2)
  expect_equal(
    calculate_mdfc(0.5, 0.5, power = 0.8, df = 40, base = 2),
    2^((qt(0.975, 40) + qt(0.8, 40)) * se)
  )
  expect_equal(
    calculate_mdfc(0.5, 0.5, power = 0.8, df = 40, base = 10),
    10^((qt(0.975, 40) + qt(0.8, 40)) * se)
  )
  expect_gt(
    calculate_mdfc(0.5, 0.5, power = 0.8, df = 40, base = exp(1)),
    calculate_mdfc(0.5, 0.5, power = 0.8, df = 40, base = 2)
  )

  expect_equal(log_base_value("log"), exp(1))
  expect_equal(log_base_value("log2"), 2)
  expect_equal(log_base_value("log10"), 10)
  expect_error(log_base_value("ln"), "Unrecognized log scale")
})


test_that('calculate_power_at_margin does not depend on the observed effect', {

  se_x <- c(0.1, 0.4, 1.2)
  se_y <- c(0.2, 0.4, 0.3)
  se <- sqrt(se_x^2 + se_y^2)
  margin_log <- log(2)

  pwr <- calculate_power_at_margin(se_x, se_y, alpha = 0.05,
                                   margin_log = margin_log, df = 40)

  # Strictly DEcreasing in SE: the noisier the contrast, the less power.
  expect_equal(order(pwr), rev(order(se)))
  expect_equal(cor(pwr, se, method = "spearman"), -1)

  # Closed form: two-sided noncentral t.
  tcrit <- qt(0.975, 40)
  expect_equal(
    pwr,
    (1 - pt(tcrit, 40, ncp = margin_log / se)) + pt(-tcrit, 40, ncp = margin_log / se)
  )

  # A bigger change is easier to detect; a stricter alpha is harder.
  expect_true(all(
    calculate_power_at_margin(se_x, se_y, margin_log = log(4), df = 40) > pwr
  ))
  expect_true(all(
    calculate_power_at_margin(se_x, se_y, alpha = 0.01,
                              margin_log = margin_log, df = 40) < pwr
  ))

  # Sign of the margin is irrelevant -- it is a magnitude.
  expect_equal(
    calculate_power_at_margin(se_x, se_y, margin_log = -margin_log, df = 40), pwr
  )

  # Bounded by alpha from below (a zero margin is just the type I error rate).
  expect_equal(
    calculate_power_at_margin(0.3, 0.3, margin_log = 0, df = 40), 0.05
  )
})


test_that('mdfc80 and power_at_margin are duals of each other', {

  # The two columns answer the same question from opposite ends: fixing the
  # power and reporting the effect, or fixing the effect and reporting the
  # power. Powering against a row's own mdfc80 must return the power that
  # mdfc80 was computed at.
  se_x <- c(0.1, 0.3, 0.8)
  se_y <- c(0.2, 0.3, 0.5)

  for (target in c(0.8, 0.9)) {
    mdfc <- calculate_mdfc(se_x, se_y, power = target, df = 40)
    back <- calculate_power_at_margin(se_x, se_y, margin_log = log(mdfc), df = 40)
    # Not exact: mdfc80 uses the (tcrit + t_power) normal-theory form rather
    # than inverting the noncentral t, which is very slightly conservative.
    expect_equal(back, rep(target, length(se_x)), tolerance = 0.02)
    expect_true(all(back >= target))
  }
})


test_that('the new columns return NA where no answer is defined', {

  # A degenerate fit has no detection limit. It must not report a 1-fold limit,
  # which is what the removed calculate_mdfc() did once the observed-power
  # contamination it accidentally depended on was taken away.
  expect_true(is.na(calculate_mdfc(0, 0, power = 0.8, df = 40)))
  expect_true(is.na(calculate_power_at_margin(0, 0, margin_log = log(2), df = 40)))

  expect_true(is.na(calculate_mdfc(NA_real_, 0.3, power = 0.8, df = 40)))
  expect_true(is.na(calculate_power_at_margin(NA_real_, 0.3, margin_log = log(2), df = 40)))

  # Invalid residual df yields no answer either, rather than falling back to a
  # normal approximation that hides the problem.
  for (bad_df in list(0, -1, NULL, NA_real_)) {
    expect_true(is.na(calculate_mdfc(0.3, 0.3, power = 0.8, df = bad_df)))
    expect_true(is.na(calculate_power_at_margin(0.3, 0.3, margin_log = log(2), df = bad_df)))
  }

  # Valid rows are untouched by the guard, and keep their position.
  mdfc80 <- calculate_mdfc(c(0, 0.3), c(0, 0.3), power = 0.8, df = 40)
  pwr <- calculate_power_at_margin(c(0, 0.3), c(0, 0.3), margin_log = log(2), df = 40)
  expect_true(is.na(mdfc80[1]) && is.finite(mdfc80[2]))
  expect_true(is.na(pwr[1]) && is.finite(pwr[2]))
})


test_that('t based detection limits are not the normal approximation', {

  # delta_p_value has always been t based while power/mdfc were z based. With
  # few samples per arm that mismatch is material, which is why mdfc80 uses t.
  se <- sqrt(0.3^2 + 0.3^2)
  z_based <- exp((qnorm(0.975) + qnorm(0.8)) * se)

  expect_gt(calculate_mdfc(0.3, 0.3, power = 0.8, df = 8) / z_based, 1.1)
  expect_equal(calculate_mdfc(0.3, 0.3, power = 0.8, df = 1e6), z_based,
               tolerance = 1e-4)
})


test_that('observed_power reproduces the 0.0.2 power column', {

  # `power` survives as a deprecated alias of `observed_power`, so its values
  # must not move. This pins them against a from-scratch restatement of the
  # 0.0.2 formula.
  set.seed(2016)
  n <- 500
  se_x <- c(0, abs(rnorm(n - 1, 0.3, 0.1)))
  se_y <- c(0, abs(rnorm(n - 1, 0.3, 0.1)))
  beta_x <- rnorm(n)
  beta_y <- rnorm(n)

  legacy_power <- local({
    Z <- (beta_x - beta_y) / sqrt(se_x^2 + se_y^2)
    Z_alpha <- qnorm(1 - 0.05 / 2)
    p <- 1 - pnorm(Z_alpha - Z) + pnorm(-Z_alpha - Z)
    ifelse(se_x == 0 & se_y == 0, 0, p)
  })
  expect_identical(
    calculate_observed_power(beta_x, se_x, beta_y, se_y, alpha = 0.05),
    legacy_power
  )
})


test_that('mdfc80 is not the old mdfc, and differs in the documented direction', {

  # calculate_mdfc() is fixed rather than frozen, so the published quantity
  # moves. This records how, so the change stays deliberate.
  set.seed(2016)
  n <- 2000
  se_x <- abs(rnorm(n, 0.3, 0.1))
  se_y <- abs(rnorm(n, 0.3, 0.1))
  beta_x <- rnorm(n)
  beta_y <- rnorm(n)
  se <- sqrt(se_x^2 + se_y^2)

  obs_power <- calculate_observed_power(beta_x, se_x, beta_y, se_y, alpha = 0.05)

  # The 0.0.2 value: the observed power vector fed in where a requested power
  # level belonged, and exp() regardless of scale.
  legacy_mdfc <- exp((qnorm(0.975) + qnorm(obs_power)) * se)
  mdfc80 <- calculate_mdfc(se_x, se_y, alpha = 0.05, power = 0.8, df = 40)

  # It moves nearly everywhere, and understates on the clear majority of rows.
  moved <- mean(abs(legacy_mdfc - mdfc80) > 1e-9, na.rm = TRUE)
  expect_gt(moved, 0.99)
  finite <- is.finite(legacy_mdfc) & is.finite(mdfc80)
  expect_gt(mean(legacy_mdfc[finite] < mdfc80[finite]), 0.5)

  # The old value tracked the observed effect; the new one tracks precision.
  expect_lt(abs(cor(legacy_mdfc[finite], se[finite], method = "spearman")), 0.5)
  expect_equal(cor(mdfc80, se, method = "spearman"), 1)

  # The old Inf rows were the MOST significant ones (qnorm(1) = Inf), not the
  # degenerate ones. mdfc80 is finite there and NA only for a degenerate fit.
  inf_rows <- !is.finite(legacy_mdfc)
  if (any(inf_rows)) {
    expect_true(all(obs_power[inf_rows] > 0.99))
    expect_true(all(is.finite(mdfc80[inf_rows])))
  }
  expect_true(is.na(calculate_mdfc(0, 0, alpha = 0.05, power = 0.8, df = 40)))
})
