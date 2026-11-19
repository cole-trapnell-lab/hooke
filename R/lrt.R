# Likelihood-ratio test for an abundance contrast, one cell group at a time.
#
# With a spherical or diagonal covariance, a PLN model is one Poisson
# log-normal regression per cell group:
#
#   y_i ~ Poisson(exp(offset_i + x_i' beta + sigma * z_i)),  z_i ~ N(0, 1)
#
# and with a full covariance each cell group's marginal distribution still has
# this form. Its likelihood is a one-dimensional integral over z, computed here
# by Gauss-Hermite quadrature, so the test below is an ordinary likelihood-ratio
# test on that marginal likelihood. It does not use the variational bound that
# PLNmodels maximizes.
#
# Everything in this file is internal. compare_abundances(test = "lrt_*") is the
# entry point.

# Gauss-Hermite rule for an expectation over a standard normal:
#   E[g(z)] ~= sum(exp(logw) * g(nodes)).
# Golub-Welsch: nodes are the eigenvalues of the Jacobi matrix of the Hermite
# recurrence, weights the squared first components of its eigenvectors.
gauss_hermite_normal <- function(n = 20) {
  i <- seq_len(n - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1)] <- sqrt(i / 2)
  jacobi[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(jacobi, symmetric = TRUE)
  o <- order(e$values)
  list(nodes = sqrt(2) * e$values[o], logw = log(e$vectors[1, o]^2))
}

# Log-likelihood of one cell group and its gradient in (beta, log sigma).
pln1_loglik_grad <- function(beta, log_sigma, y, lgamma_y, X, off, gh) {
  s <- exp(log_sigma)
  eta <- off + as.vector(X %*% beta)
  e <- outer(eta, s * gh$nodes, "+")
  lam <- exp(e)
  lp <- y * e - lam - lgamma_y
  lp <- sweep(lp, 2, gh$logw, "+")
  m <- lp[cbind(seq_along(y), max.col(lp, ties.method = "first"))]
  w <- exp(lp - m)
  rs <- rowSums(w)
  r <- w * (y - lam) / rs
  list(
    ll = sum(m + log(rs)),
    grad_beta = as.vector(crossprod(X, rowSums(r))),
    grad_log_sigma = s * sum(r %*% gh$nodes)
  )
}

# Bounds on the optimizer. Coefficients are on design columns rescaled to
# max |x| = 1, so 25 is a 7e10-fold change; sigma is the log-scale SD between
# samples beyond Poisson noise.
.pln1_beta_bound <- 25
.pln1_sigma_range <- c(0.01, 5)

# Maximum-likelihood fit of one cell group. `sigma_fixed = NULL` estimates
# sigma; otherwise sigma is held at that value.
pln1_fit <- function(y, lgamma_y, X, off, gh, beta0, sigma0 = 0.3, sigma_fixed = NULL) {
  k <- ncol(X)
  beta0[!is.finite(beta0)] <- 0
  beta0 <- pmin(pmax(beta0, -.pln1_beta_bound), .pln1_beta_bound)

  if (!is.null(sigma_fixed)) {
    ls <- log(sigma_fixed)
    cache <- new.env()
    cache$p <- NULL
    ev <- function(p) {
      if (is.null(cache$p) || any(p != cache$p)) {
        cache$p <- p
        cache$r <- pln1_loglik_grad(p, ls, y, lgamma_y, X, off, gh)
      }
      cache$r
    }
    o <- stats::optim(beta0, function(p) -ev(p)$ll, function(p) -ev(p)$grad_beta,
      method = "L-BFGS-B",
      lower = rep(-.pln1_beta_bound, k), upper = rep(.pln1_beta_bound, k)
    )
    return(list(ll = -o$value, beta = o$par, sigma = sigma_fixed, convergence = o$convergence))
  }

  if (!is.finite(sigma0)) sigma0 <- 0.3
  ls0 <- log(min(max(sigma0, .pln1_sigma_range[1] * 1.01), .pln1_sigma_range[2] * 0.998))
  cache <- new.env()
  cache$p <- NULL
  ev <- function(p) {
    if (is.null(cache$p) || any(p != cache$p)) {
      cache$p <- p
      cache$r <- pln1_loglik_grad(p[seq_len(k)], p[k + 1], y, lgamma_y, X, off, gh)
    }
    cache$r
  }
  o <- stats::optim(c(beta0, ls0), function(p) -ev(p)$ll,
    function(p) {
      r <- ev(p)
      -c(r$grad_beta, r$grad_log_sigma)
    },
    method = "L-BFGS-B",
    lower = c(rep(-.pln1_beta_bound, k), log(.pln1_sigma_range[1])),
    upper = c(rep(.pln1_beta_bound, k), log(.pln1_sigma_range[2]))
  )
  list(ll = -o$value, beta = o$par[seq_len(k)], sigma = exp(o$par[k + 1]), convergence = o$convergence)
}

# Poisson GLM coefficients as a starting point.
pln1_start <- function(y, X, off) {
  tryCatch(
    {
      b <- suppressWarnings(stats::glm.fit(X, y, offset = off, family = stats::poisson())$coefficients)
      pmin(pmax(b, -20), 20)
    },
    error = function(e) rep(0, ncol(X))
  )
}

# Likelihood-ratio statistic for H0: contrast' beta = 0 in one cell group.
#
# The reduced model is the full design times a basis of the contrast's null
# space, so it is nested in the full model. Each fit is restarted from the
# other's optimum (the full fit from the reduced optimum, which it contains;
# the reduced fit from the projected full optimum), so the statistic cannot be
# negative because one optimizer stopped early.
pln1_lrt <- function(y, X, off, null_basis, gh, sigma_fixed = NULL) {
  lgamma_y <- lgamma(y + 1)
  best_of <- function(a, b) if (b$ll > a$ll) b else a
  fit <- function(Xm, beta0, sigma0 = 0.3) {
    pln1_fit(y, lgamma_y, Xm, off, gh, beta0, sigma0 = sigma0, sigma_fixed = sigma_fixed)
  }
  Xr <- X %*% null_basis
  full <- fit(X, pln1_start(y, X, off))
  red <- fit(Xr, pln1_start(y, Xr, off))
  full <- best_of(full, fit(X, as.vector(null_basis %*% red$beta), red$sigma))
  red <- best_of(red, fit(Xr, stats::lm.fit(Xr, X %*% full$beta)$coefficients, full$sigma))
  full <- best_of(full, fit(X, as.vector(null_basis %*% red$beta), red$sigma))
  list(statistic = max(0, 2 * (full$ll - red$ll)), sigma = full$sigma)
}

# Likelihood-ratio test of one contrast in every cell group.
#
# Y: samples x cell groups counts. X: samples x coefficients design matrix.
# off: log offsets, one per sample, or a samples x cell groups matrix (the
#   layout PLNmodels::prepare_data() gives). contrast: length ncol(X).
# sigma_fixed: NULL (sigma estimated in both fits, per cell group) or a named
#   vector of the sigma to hold fixed for each cell group.
# reference: "F" compares the statistic to F(1, df2), where df2 is the residual
#   degrees of freedom (samples minus coefficients, minus one more when sigma
#   is estimated); "chisq" compares it to chi-square(1).
#
# Returns a data frame with one row per cell group: lrt_statistic, lrt_p_value
# and lrt_note (NA when the test ran; otherwise why there is no p-value).
pln_contrast_lrt <- function(Y, X, off, contrast, sigma_fixed = NULL,
                             reference = c("F", "chisq"), cores = 1, n_nodes = 20) {
  reference <- match.arg(reference)
  groups <- colnames(Y)
  out <- data.frame(
    cell_group = groups, lrt_statistic = NA_real_, lrt_p_value = NA_real_,
    lrt_note = NA_character_, stringsAsFactors = FALSE
  )

  # Rescale design columns to max |x| = 1 so a bounded coefficient cannot
  # overflow exp(); the contrast is rescaled to match.
  d <- 1 / pmax(apply(abs(X), 2, max), 1e-12)
  Xs <- sweep(X, 2, d, "*")
  cs <- contrast * d

  # Aliased design columns are dropped. A contrast outside the row space of the
  # design has no estimate (for example a condition with no samples).
  q <- qr(Xs)
  if (qr(rbind(Xs, cs))$rank > q$rank) {
    out$lrt_note <- "not_estimable"
    return(out)
  }
  keep <- sort(q$pivot[seq_len(q$rank)])
  Xs <- Xs[, keep, drop = FALSE]
  cs <- cs[keep]
  if (all(abs(cs) < 1e-12)) {
    out$lrt_note <- "not_estimable"
    return(out)
  }

  n <- nrow(Xs)
  k <- ncol(Xs)
  df2 <- n - k - if (is.null(sigma_fixed)) 1 else 0
  if (reference == "F" && df2 < 1) {
    out$lrt_note <- "insufficient_df"
    return(out)
  }

  if (is.matrix(off)) {
    stopifnot(nrow(off) == n, ncol(off) == length(groups))
  } else {
    stopifnot(length(off) == n)
  }
  null_basis <- qr.Q(qr(matrix(cs, ncol = 1)), complete = TRUE)[, -1, drop = FALSE]
  gh <- gauss_hermite_normal(n_nodes)

  one <- function(j) {
    y <- Y[, j]
    if (sum(y) <= 0) {
      return(list(note = "no_counts"))
    }
    s <- if (is.null(sigma_fixed)) NULL else sigma_fixed[[groups[j]]]
    if (!is.null(s) && (!is.finite(s) || s <= 0)) {
      return(list(note = "fit_failed"))
    }
    off_j <- if (is.matrix(off)) off[, j] else off
    tryCatch(pln1_lrt(y, Xs, off_j, null_basis, gh, sigma_fixed = s),
      error = function(e) list(note = "fit_failed")
    )
  }
  res <- if (cores > 1) {
    parallel::mclapply(seq_along(groups), one, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_along(groups), one)
  }

  for (j in seq_along(groups)) {
    r <- res[[j]]
    if (!is.list(r) || (is.null(r$statistic) && is.null(r$note))) {
      out$lrt_note[j] <- "fit_failed"
    } else if (!is.null(r$note)) {
      out$lrt_note[j] <- r$note
    } else {
      out$lrt_statistic[j] <- r$statistic
    }
  }
  out$lrt_p_value <- if (reference == "F") {
    stats::pf(out$lrt_statistic, 1, df2, lower.tail = FALSE)
  } else {
    stats::pchisq(out$lrt_statistic, 1, lower.tail = FALSE)
  }
  out
}

# The counts, design matrix and log offsets the full model was fitted to,
# rebuilt from the model's cell count set the way new_cell_count_model() built
# them.
ccm_full_model_data <- function(ccm) {
  ccs <- ccm@ccs
  if (ncol(ccs) == 0) {
    stop(
      "This option needs the samples the model was fitted to, but the model's ",
      "cell count set is empty. Fit the model with keep_ccs = TRUE."
    )
  }
  pln_data <- PLNmodels::prepare_data(
    counts = counts(ccs),
    covariates = colData(ccs) %>% as.data.frame(),
    offset = monocle3::size_factors(ccs)
  )
  full_model_formula <- stats::as.formula(ccm@full_model_formula)
  frame <- stats::model.frame(full_model_formula[-2], pln_data)
  X <- stats::model.matrix(stats::terms(frame), frame)
  # prepare_data() stores the offset as a samples x cell groups matrix, so the
  # model frame's offset is one too.
  off <- stats::model.offset(frame)
  stored <- as.numeric(ccm@model_aux[["full_model_offsets"]])
  if (length(stored) != length(off) || !isTRUE(all.equal(stored, as.numeric(off)))) {
    stop(
      "The model's cell count set does not match the samples the model was ",
      "fitted to (offsets differ), so the contrast cannot be refitted."
    )
  }
  Y <- as.matrix(pln_data$Abundance)
  if (is.matrix(off)) {
    if (!identical(dim(off), dim(Y))) {
      stop("The model's offsets do not match its counts (", paste(dim(off), collapse = " x "), " vs ", paste(dim(Y), collapse = " x "), ").")
    }
  } else {
    off <- as.numeric(off)
  }
  list(Y = Y, X = X, off = off)
}

# The single design-matrix row behind a table from estimate_abundances().
condition_design_row <- function(ccm, cond, which_cond = "cond") {
  model_terms <- ccm@model_aux[["full_model_terms"]]
  xlev <- ccm@model_aux[["full_model_xlevels"]]
  vars <- setdiff(all.vars(model_terms), "Offset")
  missing_vars <- setdiff(vars, colnames(cond))
  if (length(missing_vars) > 0) {
    stop(
      which_cond, " is missing model variable(s): ", paste(missing_vars, collapse = ", "),
      ". Pass the table returned by estimate_abundances()."
    )
  }
  newdata <- unique(as.data.frame(cond)[, vars, drop = FALSE])
  if (nrow(newdata) != 1) {
    stop(which_cond, " must describe exactly one condition; it has ", nrow(newdata), ".")
  }
  newdata$Offset <- 1
  newdata <- coerce_newdata_to_xlevels(newdata, xlev)
  X <- stats::model.matrix(model_terms, newdata, xlev = xlev)
  stats::setNames(as.numeric(X[1, ]), colnames(X))
}

# Standard error of contrast' beta for every cell group, from the model's
# coefficient covariance (natural-log scale). Unlike
# sqrt(se_x^2 + se_y^2), this keeps the covariance between the two predictions,
# so uncertainty they share (a batch term, a shared time curve) cancels.
contrast_se_from_vhat <- function(ccm, contrast) {
  se <- contrast_se_blocks(ccm@vhat, contrast, model(ccm)$p)
  stats::setNames(se, colnames(model(ccm)$model_par$Sigma))
}

# v_hat: (groups x coefficients)-square covariance, one block of coefficients
# per cell group, in cell-group order (the layout of a cell_count_model's vhat).
contrast_se_blocks <- function(v_hat, contrast, p) {
  v_hat[is.na(v_hat)] <- 0
  k <- length(contrast)
  if (nrow(v_hat) != p * k) {
    stop("The model's coefficient covariance does not match its design (", nrow(v_hat), " vs ", p, " x ", k, ").")
  }
  C <- Matrix::bdiag(rep.int(list(Matrix::Matrix(matrix(contrast, nrow = 1), sparse = TRUE)), p))
  v <- Matrix::diag(as.matrix(C %*% v_hat %*% Matrix::t(C)))
  sqrt(pmax(v, 0))
}

# Samples and cells behind one condition: the samples whose values in the
# `support_by` columns equal the condition's.
condition_support <- function(ccm, cond, support_by, which_cond = "cond") {
  ccs <- ccm@ccs
  if (ncol(ccs) == 0) {
    stop("support_by needs the model's cell count set. Fit the model with keep_ccs = TRUE.")
  }
  cd <- as.data.frame(colData(ccs))
  missing_cols <- setdiff(support_by, intersect(colnames(cd), colnames(cond)))
  if (length(missing_cols) > 0) {
    stop("support_by column(s) not in ", which_cond, " and the cell count set: ", paste(missing_cols, collapse = ", "))
  }
  vals <- unique(as.data.frame(cond)[, support_by, drop = FALSE])
  if (nrow(vals) != 1) {
    stop(which_cond, " must describe exactly one condition; it has ", nrow(vals), ".")
  }
  sel <- rep(TRUE, nrow(cd))
  for (nm in support_by) {
    sel <- sel & !is.na(cd[[nm]]) & as.character(cd[[nm]]) == as.character(vals[[nm]])
  }
  cnt <- counts(ccs)
  cells <- if (any(sel)) Matrix::rowSums(cnt[, sel, drop = FALSE]) else rep(0, nrow(cnt))
  list(n_samples = sum(sel), cells = stats::setNames(as.numeric(cells), rownames(cnt)))
}
