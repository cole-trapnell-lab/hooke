context('test-lrt')

# Likelihood-ratio test and the other opt-in calling options of
# compare_abundances(). These test the internals on simulated counts; none of
# them fits a PLN model.

sim_design <- function(n_per = 8, stages = c(24, 48, 72)) {
  d <- expand.grid(rep = seq_len(n_per), knockout = c(FALSE, TRUE), timepoint = stages)
  X <- model.matrix(~ factor(timepoint) * knockout, d)
  list(d = d, X = X)
}

sim_counts <- function(X, beta, sigma, off, seed) {
  set.seed(seed)
  rpois(nrow(X), exp(off + as.vector(X %*% beta) + rnorm(nrow(X), 0, sigma)))
}

test_that('the Gauss-Hermite rule integrates against a standard normal', {
  gh <- gauss_hermite_normal(20)
  w <- exp(gh$logw)
  expect_equal(sum(w), 1)
  expect_equal(sum(w * gh$nodes), 0)
  expect_equal(sum(w * gh$nodes^2), 1)
  expect_equal(sum(w * gh$nodes^4), 3)
  expect_equal(sum(w * exp(gh$nodes)), exp(0.5), tolerance = 1e-10)
})

test_that('the analytic gradient matches finite differences', {
  s <- sim_design()
  off <- rep(0.1, nrow(s$X))
  beta <- c(2, 0.3, -0.2, -0.5, 0.1, 0.2)
  y <- sim_counts(s$X, beta, 0.4, off, seed = 1)
  gh <- gauss_hermite_normal(20)
  f <- function(p) pln1_loglik_grad(p[1:6], p[7], y, lgamma(y + 1), s$X, off, gh)$ll
  p0 <- c(beta, log(0.4))
  num <- sapply(seq_along(p0), function(i) {
    h <- 1e-5
    e <- replace(numeric(7), i, h)
    (f(p0 + e) - f(p0 - e)) / (2 * h)
  })
  g <- pln1_loglik_grad(p0[1:6], p0[7], y, lgamma(y + 1), s$X, off, gh)
  expect_equal(c(g$grad_beta, g$grad_log_sigma), num, tolerance = 1e-5)
})

test_that('with almost no extra variance the test reduces to the Poisson GLM test', {
  s <- sim_design()
  off <- rep(0, nrow(s$X))
  y <- sim_counts(s$X, c(3, 0.2, -0.1, -0.4, 0, 0), 0, off, seed = 2)
  Y <- matrix(y, ncol = 1, dimnames = list(NULL, "a"))
  # knockout effect at the first stage
  contrast <- c(0, 0, 0, 1, 0, 0)
  r <- pln_contrast_lrt(Y, s$X, off, contrast, sigma_fixed = c(a = 0.01), reference = "chisq")
  full <- glm(y ~ s$X - 1, family = poisson())
  red <- glm(y ~ s$X[, -4] - 1, family = poisson())
  expect_equal(r$lrt_statistic, deviance(red) - deviance(full), tolerance = 0.02)
})

test_that('a near-total loss is called and an unchanged cell group is not', {
  s <- sim_design()
  off <- rep(0, nrow(s$X))
  contrast <- c(0, 0, 1, 1, 0, 1)   # knockout minus control at the last stage
  set.seed(3)
  unchanged <- sim_counts(s$X, c(3, 0.2, 0.1, 0, 0, 0), 0.3, off, seed = 4)
  lost <- sim_counts(s$X, c(3, 0.2, 0.1, 0, 0, 0), 0.3, off, seed = 5)
  lost[s$d$knockout & s$d$timepoint == 72] <- 0
  Y <- cbind(unchanged = unchanged, lost = lost, empty = 0)
  contrast <- c(0, 0, 0, 1, 0, 1)
  for (sf in list(NULL, c(unchanged = 0.3, lost = 0.3, empty = 0.3))) {
    r <- pln_contrast_lrt(Y, s$X, off, contrast, sigma_fixed = sf)
    expect_equal(r$cell_group, c("unchanged", "lost", "empty"))
    expect_true(all(r$lrt_statistic[1:2] >= 0))
    expect_lt(r$lrt_p_value[2], 1e-6)
    expect_gt(r$lrt_p_value[1], 0.01)
    expect_equal(r$lrt_note, c(NA, NA, "no_counts"))
    expect_true(is.na(r$lrt_p_value[3]))
  }
})

test_that('offsets may be one per sample or a samples x cell groups matrix', {
  s <- sim_design()
  off <- log(runif(nrow(s$X), 0.5, 2))
  Y <- sapply(1:3, function(i) sim_counts(s$X, c(3, 0.2, 0.1, -0.6, 0, 0), 0.3, off, seed = 20 + i))
  colnames(Y) <- paste0("c", 1:3)
  contrast <- c(0, 0, 0, 1, 0, 0)
  v <- pln_contrast_lrt(Y, s$X, off, contrast)
  m <- pln_contrast_lrt(Y, s$X, matrix(off, nrow(Y), ncol(Y)), contrast)
  expect_equal(v, m)
  expect_true(all(is.na(v$lrt_note)))
  expect_error(pln_contrast_lrt(Y, s$X, rep(off, 2), contrast))
})

test_that('the F reference is never less conservative than chi-square', {
  s <- sim_design(n_per = 4)
  off <- rep(0, nrow(s$X))
  Y <- sapply(1:5, function(i) sim_counts(s$X, c(2, 0.2, 0.1, 0.5, 0, 0), 0.4, off, seed = 10 + i))
  colnames(Y) <- paste0("c", 1:5)
  contrast <- c(0, 0, 0, 1, 0, 0)
  f <- pln_contrast_lrt(Y, s$X, off, contrast, reference = "F")
  ch <- pln_contrast_lrt(Y, s$X, off, contrast, reference = "chisq")
  expect_equal(f$lrt_statistic, ch$lrt_statistic)
  expect_true(all(f$lrt_p_value >= ch$lrt_p_value))
})

test_that('a contrast the design cannot estimate is reported, not guessed', {
  s <- sim_design()
  keep <- !(s$d$knockout & s$d$timepoint == 72)     # no knockout samples at the last stage
  X <- s$X[keep, ]
  off <- rep(0, nrow(X))
  Y <- matrix(sim_counts(X, c(3, 0, 0, 0, 0, 0), 0.3, off, seed = 6), ncol = 1, dimnames = list(NULL, "a"))
  r <- pln_contrast_lrt(Y, X, off, c(0, 0, 0, 1, 0, 1))
  expect_equal(r$lrt_note, "not_estimable")
  expect_true(is.na(r$lrt_p_value))
  # the same design still supports the contrast at a stage that has both arms
  r <- pln_contrast_lrt(Y, X, off, c(0, 0, 0, 1, 0, 0))
  expect_true(is.na(r$lrt_note))
  expect_true(is.finite(r$lrt_p_value))
})

test_that('the covariance standard error cancels what the two estimates share', {
  # two cell groups, coefficients (intercept, batch, knockout)
  v1 <- matrix(c(0.04, 0, 0, 0, 1.00, 0, 0, 0, 0.09), 3)   # a very uncertain batch term
  v2 <- diag(c(0.01, 0.25, 0.16))
  v_hat <- as.matrix(Matrix::bdiag(v1, v2))
  x_ctrl <- c(1, 1, 0)
  x_ko <- c(1, 1, 1)
  se_cov <- contrast_se_blocks(v_hat, x_ko - x_ctrl, p = 2)
  expect_equal(se_cov, c(0.3, 0.4))
  se_indep <- c(
    sqrt(x_ctrl %*% v1 %*% x_ctrl + x_ko %*% v1 %*% x_ko),
    sqrt(x_ctrl %*% v2 %*% x_ctrl + x_ko %*% v2 %*% x_ko)
  )
  expect_true(all(se_indep > se_cov))
  expect_error(contrast_se_blocks(v_hat, c(1, 0), p = 2), "does not match")
})

# The vhat a bootstrap leaves when one resample had no samples of one condition
# (hooke#65): the knockout coefficient's row and column are NA in every block.
na_knockout_vhat <- function() {
  v <- diag(c(0.04, 0.09))                      # coefficients (intercept, knockout)
  v[2, ] <- NA
  v[, 2] <- NA
  methods::as(Matrix::bdiag(v, v), "CsparseMatrix")   # two cell groups, sparse like a ccm's vhat
}

test_that('a contrast whose coefficient has no variance gets an NA standard error, not 0', {
  v_hat <- na_knockout_vhat()
  expect_equal(contrast_se_blocks(v_hat, c(0, 1), p = 2), c(NA_real_, NA_real_))
  # a contrast that does not use the NA coefficient keeps its standard error
  expect_equal(contrast_se_blocks(v_hat, c(1, 0), p = 2), c(0.2, 0.2))
})

test_that('only the predictions that use an NA coefficient are flagged', {
  v_hat <- na_knockout_vhat()
  expect_equal(uses_na_coefficient(v_hat, c(1, 0), p = 2), c(FALSE, FALSE))   # control arm
  expect_equal(uses_na_coefficient(v_hat, c(1, 1), p = 2), c(TRUE, TRUE))     # perturbed arm
  expect_equal(uses_na_coefficient(v_hat, c(0, 0), p = 2), c(FALSE, FALSE))
  # an NA in one cell group's block flags that group only
  v <- as.matrix(Matrix::bdiag(diag(c(0.04, 0.09)), diag(c(0.01, 0.16))))
  v[4, 4] <- NA
  expect_equal(uses_na_coefficient(v, c(1, 1), p = 2), c(FALSE, TRUE))
  expect_equal(contrast_se_blocks(v, c(0, 1), p = 2), c(0.3, NA_real_))
  expect_error(uses_na_coefficient(v, c(1, 0, 0), p = 2), "does not match")
})

test_that('a covariance without NA gives the same standard errors as before', {
  v1 <- matrix(c(0.04, 0, 0, 0, 1.00, 0, 0, 0, 0.09), 3)
  v_hat <- as.matrix(Matrix::bdiag(v1, diag(c(0.01, 0.25, 0.16))))
  expect_equal(uses_na_coefficient(v_hat, c(1, 1, 1), p = 2), c(FALSE, FALSE))
  expect_equal(contrast_se_blocks(v_hat, c(0, 0, 1), p = 2), c(0.3, 0.4))
})

test_that('the new compare_abundances() arguments default to the old behaviour', {
  f <- formals(compare_abundances)
  expect_equal(eval(f$contrast_se)[1], "independent")
  expect_equal(eval(f$test)[1], "wald")
  expect_null(f$support_by)
  expect_null(f$q_value_groups)
})
