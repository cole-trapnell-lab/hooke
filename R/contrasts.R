# NOTE: This ignores offsets. DO NOT USE TO RECOVER FITTED VALUES FROM THE ORIGINAL DATA
#' @noRd
my_plnnetwork_predict <- function(ccm, newdata, type = c("link", "response"), envir = parent.frame()) {
  type <- match.arg(type)
  n_new <- nrow(newdata)
  X <- model.matrix(ccm@model_aux[["full_model_terms"]], newdata,
    xlev = ccm@model_aux[["full_model_xlevels"]]
  )

  # #O <- model.offset(ccm@model_aux[["full_model_frame"]])
  EZ <- tcrossprod(X, t(model(ccm)$model_par$B))
  # #if (!is.null(O))
  # #  EZ <- EZ + O
  # EZ <- sweep(EZ, 2, 0.5 * Matrix::diag(model(ccm)$model_par$Sigma), "+")

  colnames(EZ) <- colnames(model(ccm)$model_par$Sigma)

  M <- matrix(0, nrow = n_new, ncol = model(ccm)$p)
  S <- matrix(Matrix::diag(model(ccm)$model_par$Sigma), nrow = n_new, ncol = model(ccm)$p, byrow = TRUE)

  # results <- switch(type, link = EZ, response = exp(EZ))
  results <- switch(type,
    link = EZ + M,
    response = exp(EZ + M + 0.5 * S)
  )
  attr(results, "type") <- type
  results
}



#' Predict Conditional Means from PLN Model
#'
#' Computes conditional predictions from a PLN (Poisson Log-Normal) model, given new data and conditional responses.
#'
#' @param ccm An object containing the fitted PLN model and associated data structures.
#' @param newdata A data frame or matrix of new covariate values for prediction.
#' @param cond_responses A matrix or data frame of observed responses for the conditional species (columns must be a subset of the species in the model).
#' @param type Character string specifying the scale of the predictions. Either \code{"link"} (default, latent scale) or \code{"response"} (mean response scale).
#' @param var_par Logical; if \code{TRUE}, attaches the conditional mean (\code{M}) and variance (\code{S}) as attributes to the result.
#' @param envir Environment in which to evaluate model terms (default is parent frame).
#' @param pln_model Character string specifying which PLN model to use for prediction. Either \code{"full"} (default) or \code{"reduced"}.
#'
#' @return A matrix of predicted values for the unconditioned species. The type of prediction is determined by the \code{type} argument. If \code{var_par = TRUE}, the result has attributes \code{"M"} (conditional mean) and \code{"S"} (conditional variance).
#'
#' @details
#' This function computes the conditional mean (and optionally variance) of the latent variables or response for a subset of species, given observed values for another subset (the conditional responses), under a fitted PLN model. It handles both the full and reduced model forms, and can return predictions on the latent or response scale.
my_pln_predict_cond <- function(ccm,
                                newdata,
                                cond_responses,
                                type = c("link", "response"),
                                var_par = FALSE,
                                envir = parent.frame(),
                                pln_model = c("full", "reduced")) {
  type <- match.arg(type)
  pln_model <- match.arg(pln_model)

  # Checks
  Yc <- as.matrix(cond_responses)
  sp_names <- colnames(model(ccm, model_to_return = pln_model)$model_par$B)
  if (!any(colnames(cond_responses) %in% sp_names)) {
    stop("Yc must be a subset of the species in responses")
  }
  if (!nrow(Yc) == nrow(newdata)) {
    stop("The number of rows of Yc must match the number of rows in newdata")
  }

  # Dimensions and subsets
  n_new <- nrow(Yc)
  cond <- sp_names %in% colnames(Yc)
  cond_2 <- colnames(Yc) %in% sp_names
  ## Extract the model matrices from the new data set with initial formula
  # X <- model.matrix(formula(private$formula)[-2], newdata, xlev = attr(private$formula, "xlevels"))
  # X <- model.matrix(terms(ccm@model_aux[["model_frame"]]), newdata,
  # xlev = ccm@model_aux[["xlevels"]])

  if (pln_model == "full") {
    # X <- model.matrix(terms(ccm@model_aux[["full_model_frame"]]), newdata,
    #                   xlev = ccm@model_aux[["full_model_xlevels"]])
    X <- model.matrix(ccm@model_aux[["full_model_terms"]], newdata,
      xlev = ccm@model_aux[["full_model_xlevels"]]
    )
  } else if (pln_model == "reduced") {
    # X <- model.matrix(terms(ccm@model_aux[["reduced_model_frame"]]), newdata,
    #                   xlev = ccm@model_aux[["reduced_model_xlevels"]])
    X <- model.matrix(ccm@model_aux[["reduced_model_terms"]], newdata,
      xlev = ccm@model_aux[["reduced_model_xlevels"]]
    )
  }

  # O <- model.offset(model.frame(formula(ccm@full_model_formula)[-2], newdata))
  O <- NULL
  # O <- model.offset(model.frame(formula(private$formula)[-2], newdata))
  if (is.null(O)) {
    O <- matrix(0, n_new, model(ccm, model_to_return = pln_model)$p)
  }

  # Compute parameters of the law
  Sigma <- model(ccm, model_to_return = "reduced")$model_par$Sigma
  vcov11 <- Sigma[cond, cond, drop = FALSE]
  vcov22 <- Sigma[!cond, !cond, drop = FALSE]
  vcov12 <- Sigma[cond, !cond, drop = FALSE]
  prec11 <- solve(vcov11)

  # A <- crossprod(vcov12, prec11)
  A <- crossprod(as.matrix(Sigma[cond, , drop = FALSE]), prec11)

  # Sigma21 <- vcov22 - A %*% vcov12
  Sigma21 <- as.matrix(Sigma[, , drop = FALSE]) - A %*% as.matrix(Sigma[cond, , drop = FALSE])

  VE <- model(ccm, model_to_return = pln_model)$optimize_vestep(
    covariates = X,
    offsets = O[, cond, drop = FALSE],
    responses = Yc[, cond_2, drop = FALSE],
    weights = rep(1, n_new),
    B = model(ccm, model_to_return = pln_model)$model_par$B[, cond, drop = FALSE],
    Omega = prec11
  )

  M <- tcrossprod(VE$M, A)

  S <- map(1:n_new, ~ crossprod(VE$S[., ] * t(A)) + Sigma21) %>% simplify2array()

  ## mean latent positions in the parameter space

  EZ <- tcrossprod(X, t(model(ccm, model_to_return = pln_model)$model_par$B[, , drop = FALSE])) + M + O[, , drop = FALSE]
  # EZ <- sweep(EZ, 2, 0.5 * Matrix::diag(model(ccm, model_to_return = pln_model)$model_par$Sigma[, , drop = FALSE]), "+")
  colnames(EZ) <- colnames(model(ccm, model_to_return = pln_model)$model_par$Sigma[, , drop = FALSE])

  # EZ <- X %*% model(ccm)$model_par$B[, !cond, drop = FALSE] + M + O[, !cond, drop = FALSE]
  # colnames(EZ) <- setdiff(sp_names, colnames(Yc))

  # ! We should only add the .5*diag(S2) term only if we want the type="response"
  if (type == "response") {
    if (ncol(EZ) == 1) {
      EZ <- EZ + .5 * S
    } else {
      EZ <- EZ + .5 * t(apply(S, 3, diag))
    }
  }
  results <- switch(type,
    link = EZ,
    response = exp(EZ)
  )
  attr(results, "type") <- type
  if (var_par) {
    attr(results, "M") <- M
    attr(results, "S") <- S
  }
  results
}


#' Predict cell type abundances given a PLN model and a set of inputs for its covariates
#'
#' @param ccm A cell_count_model.
#' @param newdata tibble A tibble of variables used for the prediction.
#' @param min_log_abund numeric Minimum log abundance value.
#' @param cell_group string The name of the groups that are being estimated.
#' @param log_scale Desired log scale for the output. Default is natural log.
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances <- function(ccm,
                                newdata,
                                min_log_abund = -5,
                                cell_group = "cell_group",
                                scale = c("log", "log10", "log2", "per_1000")) {
  if (!tibble::is_tibble(newdata)) {
    newdata <- newdata %>% as_tibble()
  }

  scale <- match.arg(scale)
  if (scale == "per_1000") {
    type <- "response"
  } else {
    type <- "link"
  }

  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(tibble::is_tibble(newdata))
  assertthat::assert_that(is.numeric(min_log_abund))
  assertthat::assert_that(is.character(cell_group))

  # check that all terms in new data have been specified
  # missing_terms = setdiff(names(ccm@model_aux$xlevels), names(newdata))
  #
  # if (length(missing_terms) >= 1) {
  #
  #   default_df = lapply(missing_terms, function(term){
  #     df = data.frame(t = levels(factor(colData(ccm@ccs)[[term]]))[1])
  #     names(df) = term
  #     df
  #   }) %>% bind_cols()
  #
  #   newdata = cbind(newdata, tibble(default_df))
  #
  #   print( paste0(paste(missing_terms,collapse = ", "),
  #                 " missing from specified newdata columns. Assuming default values: ",
  #                 paste(default_df[1,],collapse = ", ")))
  #
  #
  # }

  newdata <- fill_missing_terms_with_default_values(ccm, newdata, pln_model = "full")

  estimate_abundance_row <- function(ccm, model_terms, newdata, min_log_abund, type) {
    # assertthat::assert_that(
    #   tryCatch(expr = length(missing_terms) == 0,
    #            error = function(e) FALSE),
    #   msg = paste0(missing_terms, " missing from newdata columns"))

    # stopifnot(nrow(newdata) == 1)
    newdata$Offset <- 1

    # model_terms = terms(ccm@model_aux[["full_model_frame"]])
    model_terms <- ccm@model_aux[["full_model_terms"]]
    base_X <- Matrix::sparse.model.matrix(model_terms, newdata,
      xlev = ccm@model_aux[["full_model_xlevels"]]
    )

    # base_X <- model.matrix(formula(ccm@model_formula_str)[-2], newdata,
    #                  xlev = ccm@model_aux[["xlevels"]])
    X <- Matrix::bdiag(rep.int(list(base_X), model(ccm)$p))

    # if it doesn't exist, use orig computation
    # if (is.na(ccm@bootstrapped_vhat)[[1]]) {
    #   v_hat = compute_vhat(ccm)
    # } else {
    #   v_hat = ccm@bootstrapped_vhat
    # }

    # vhat_coef <- coef(model(ccm), type="main")

    # vcov_type <- grep('vcov', names(attributes(vhat_coef)), value=TRUE)
    v_hat <- ccm@vhat
    v_hat[is.na(v_hat)] <- 0
    v_hat_method <- ccm@vhat_method

    se_fit <- sqrt(Matrix::diag(as.matrix(X %*% v_hat %*% Matrix::t(X))))

    # if (v_hat_method == "wald") {
    #   se_fit = sqrt(Matrix::diag(as.matrix(X %*% v_hat %*% Matrix::t(X)))) / sqrt(model(ccm)$n)
    # } else {
    #   se_fit = sqrt(Matrix::diag(as.matrix(X %*% v_hat %*% Matrix::t(X))))
    # }

    pred_out <- my_plnnetwork_predict(ccm, newdata = newdata, type = type)
    # pred_out = max(pred_out, -5)
    # log_abund = pred_out[1,]
    log_abund <- as.numeric(pred_out)

    log_abund_sd <- sqrt(Matrix::diag(coef(model(ccm), type = "covariance")))
    names(log_abund_sd) <- colnames(coef(model(ccm), type = "covariance"))
    log_abund_se <- se_fit

    below_thresh <- log_abund < min_log_abund
    log_abund[below_thresh] <- min_log_abund
    log_abund_se[below_thresh] <- 0

    # max_log_abundances = log(matrixStats::colMaxs(pln_model$fitted))
    # min_log_abundances = log(matrixStats::colMins(pln_model$fitted))
    # percent_max = 100 * (exp(log_abund)/exp(max_log_abundances))
    # percent_range = 100 * (exp(log_abund) - exp(min_log_abundances)) / (exp(max_log_abundances) - exp(min_log_abundances))
    pred_out_tbl <- tibble::tibble(
      !!cell_group := colnames(pred_out),
      log_abund,
      log_abund_se
    )
    pred_out_tbl <- left_join(pred_out_tbl, tibble::tibble(!!cell_group := names(log_abund_sd), log_abund_sd), by = cell_group)
    # max_log_abundances,
    # min_log_abundances,
    # percent_max,
    # percent_range)
    newdata$Offset <- NULL
    pred_out_tbl <- cbind(newdata, pred_out_tbl)
    pred_out_tbl <- tibble::tibble(pred_out_tbl)
  }

  pred_out_tbl <- newdata %>%
    as.data.frame() %>%
    group_split(row_number(), .keep = FALSE) %>%
    purrr::map_df(tidyr::nest) %>%
    mutate(timepoint_abund = purrr::map(
      .f = estimate_abundance_row,
      .x = data,
      ccm = ccm,
      model_terms = model_terms,
      min_log_abund = min_log_abund,
      type = type
    )) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  if (scale == "log2") {
    pred_out_tbl <- pred_out_tbl %>%
      mutate(
        log2_abund = log2(exp(log_abund)),
        log_abund_se = log2(exp(log_abund_se)),
        log_abund_se = log2(exp(log_abund_sd))
      ) %>%
      select(-c(log_abund, log_abund_se, log_abund_sd))
  } else if (scale == "log10") {
    pred_out_tbl <- pred_out_tbl %>%
      mutate(
        log_abund = log10(exp(log_abund)),
        log_abund_se = log10(exp(log_abund_se)),
        log_abund_sd = log10(exp(log_abund_sd))
      )
  } else if (scale == "per_1000") {
    pred_out_tbl <- pred_out_tbl %>%
      ## FIXME -- the output of plnnet_predict is already on the response scale. A rename may be more appropriate if it will not break any dependencies
      mutate(
        abund_per_1000 = log_abund,
        abund_per_1000_se = log_abund_se,
        abund_per_1000_sd = log_abund_sd
      )
  }

  pred_out_tbl <- tibble::tibble(pred_out_tbl)
  return(pred_out_tbl)
}

#' Predict cell type abundances given a PLN model and a set of inputs for its covariates
#' and observed counts
#'
#' @param ccm A cell_count_model.
#' @param newdata tibble A tibble of variables used for the prediction.
#  Must either be a single row or a tibble with one row per sample of the cell
#' count set for ccm.
#' @param cond_responses a data frame containing the counts of the observed variables
#' @param min_log_abund numeric Minimum log abundance value.
#' @param cell_group string The name of the groups that are being estimated.
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances_cond <- function(ccm,
                                     newdata,
                                     cond_responses,
                                     min_log_abund = -5,
                                     cell_group = "cell_group",
                                     type = c("link", "response"),
                                     pln_model = c("full", "reduced")) {
  if (!tibble::is_tibble(newdata)) {
    newdata <- newdata %>% as_tibble()
  }

  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(tibble::is_tibble(newdata))
  assertthat::assert_that(is.numeric(min_log_abund) | is.null(min_log_abund))
  assertthat::assert_that(is.character(cell_group))

  type <- match.arg(type)
  pln_model <- match.arg(pln_model)

  newdata <- fill_missing_terms_with_default_values(ccm, newdata, pln_model)
  newdata$Offset <- 1

  cond_responses <- t(as.matrix(cond_responses))
  if (nrow(cond_responses) != nrow(newdata) & nrow(newdata) == 1) {
    newdata <- newdata %>% slice(rep(1:n(), each = nrow(cond_responses)))
  } else if (nrow(cond_responses) == nrow(newdata)) {
    newdata <- newdata
  } else {
    stop("The number of rows of cond_responses must match the number of rows in newdata or the newdata must have one row.")
  }


  estimate_abundance_cond_row <- function(ccm, newdata, cond_responses,
                                          type = type,
                                          pln_model = pln_model,
                                          min_log_abund = -5) {
    newdata$Offset <- 1

    pred_out <- my_pln_predict_cond(ccm, newdata,
      cond_responses,
      type = type,
      pln_model = pln_model
    )
    # log_abund = as.numeric(pred_out)
    log_abund <- as.numeric(t(pred_out))
    newdata$Offset <- NULL

    if (is.null(min_log_abund) == FALSE) {
      below_thresh <- log_abund < min_log_abund
      log_abund[below_thresh] <- min_log_abund
      # log_abund_se[below_thresh] = 0
    }

    # pred_out_tbl = tibble::tibble(cell_group=colnames(pred_out), log_abund)
    # pred_out_tbl = cbind(newdata, pred_out_tbl)
    # pred_out_tbl <- tibble::tibble(pred_out_tbl)
    pred_out_tbl <- cbind(
      cell_group = rep(colnames(pred_out), each = length(log_abund) / ncol(pred_out)),
      log_abund,
      # log_abund_se,
      do.call("rbind", replicate(length(log_abund) / nrow(newdata), newdata, simplify = FALSE))
    )

    # pred_out_tbl = left_join(pred_out_tbl,
    #                          tibble(log_abund = log_abund,
    #                                 cell_group=rep(colnames(pred_out), times=length(log_abund)/ncol(pred_out)),
    #                                 sample=rep(newdata$sample, each=length(log_abund)/nrow(newdata))),
    #                          by=c("cell_group", "sample"))

    pred_out_tbl <- tibble::tibble(pred_out_tbl)
    return(pred_out_tbl)
  }

  pred_out_tbl <- cbind(newdata, cond_responses) %>%
    group_split(row_number(), .keep = FALSE) %>%
    purrr::map_df(tidyr::nest,
      data = colnames(newdata),
      cond_response = colnames(cond_responses)
    ) %>%
    mutate(timepoint_abund = purrr::map2(
      .f = estimate_abundance_cond_row,
      .x = data,
      .y = cond_response,
      ccm = ccm,
      type = type,
      pln_model = pln_model,
      min_log_abund = min_log_abund
    )) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  pred_out_tbl$rn <- NULL


  return(pred_out_tbl)
}


# To do : need better error message for when you are missing a column that needs to be specified for the model
#' Predict cell type abundances given a PLN model over a range of time or other interval
#' @param ccm A cell_count_model.
#' @param interval_start numeric Interval start value.
#' @param interval_stop numeric Interval stop value.
#' @param interval_col character Interval values are taken from the interval_var data. Default is "timepoint".
#' @param interval_step numeric Interval size. Default is 2.
#' @inheritParams estimate_abundances
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances_over_interval <- function(ccm,
                                              interval_start,
                                              interval_stop,
                                              interval_col = "timepoint",
                                              interval_step = 2,
                                              min_log_abund = -5,
                                              newdata = tibble(),
                                              scale = c("log", "log10", "log2", "per_1000")) {
  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(is.numeric(interval_start))
  assertthat::assert_that(is.numeric(interval_stop))
  assertthat::assert_that(interval_stop >= interval_start)
  assertthat::assert_that(is.numeric(interval_step))
  scale <- match.arg(scale)

  # assertthat::assert_that(interval_col %in% attr(terms(ccm@model_aux[['model_frame']]), 'term.labels'))

  # make it so that if new data has interval col in it, override

  timepoint_pred_df <- tibble(IV = seq(interval_start, interval_stop, interval_step))
  colnames(timepoint_pred_df)[1] <- interval_col

  if (interval_col %in% colnames(newdata)) {
    newdata[[interval_col]] <- NULL
  }

  if (nrow(newdata) > 0) {
    timepoint_pred_df <- cross_join(timepoint_pred_df, newdata)
  }

  timepoint_pred_df <- timepoint_pred_df %>%
    group_split(row_number(), .keep = FALSE) %>%
    purrr::map_df(tidyr::nest) %>%
    mutate(timepoint_abund = purrr::map(
      .f = estimate_abundances,
      .x = data,
      ccm = ccm,
      min_log_abund = min_log_abund,
      scale = scale
    )) %>%
    select(timepoint_abund) %>%
    tidyr::unnest(c(timepoint_abund))

  # time_interval_pred_helper = function(tp, ...){
  #   tp_tbl = tibble(IV=tp, ...)
  #   colnames(tp_tbl)[1] = interval_col
  #   estimate_abundances(ccm, tp_tbl, min_log_abund = min_log_abund)
  # }
  #
  # cross_join(timepoint_pred_df, tibble(expt = "GAP16"))
  #
  # timepoint_pred_df = timepoint_pred_df %>%
  #   dplyr::mutate(timepoint_abund = purrr::map(.f = purrr::possibly(
  #     .f = time_interval_pred_helper, NA_real_),
  #     .x = !!sym(interval_col),
  #     ...)) %>%
  #   select(timepoint_abund) %>%
  #   tidyr::unnest(c(timepoint_abund))

  return(timepoint_pred_df)
}


# Post-hoc ("observed") power for a contrast.
#
# Renamed from calculate_power() in 0.0.3 to say what it computes. The body is
# unchanged, so the deprecated `power` column keeps its published values.
#
# WARNING: this is not a power calculation in the design sense. It substitutes
# the OBSERVED Wald statistic for the true effect, which makes it strictly
# increasing in |Z| and therefore a deterministic restatement of
# `delta_p_value` -- `power >= 0.8` is exactly `p <= ~0.008`. It floors at
# `alpha` when Z = 0 and saturates at 1 for |Z| >= ~10.26.
#
# It does not measure precision: a cell type with a huge SE and a fluke
# estimate scores high, while a tightly measured genuine null scores the floor.
# Filtering on it before p.adjust() selects on the statistic being adjusted and
# is anti-conservative.
#
# For "were we powered to see a change here?" use calculate_power_at_margin()
# or calculate_mdfc(), both of which are effect-independent.
calculate_observed_power <- function(beta_x, SE_x, beta_y, SE_y, alpha = 0.05) {
  # Wald test statistic for comparing two groups
  Z <- (beta_x - beta_y) / sqrt(SE_x^2 + SE_y^2)

  # Critical value for a two-tailed test at alpha
  Z_alpha <- qnorm(1 - alpha / 2)

  # Compute power: 1 - Type II error probability
  power <- 1 - pnorm(Z_alpha - Z) + pnorm(-Z_alpha - Z)

  # return 0 if divide by 0
  power <- ifelse(SE_x == 0 & SE_y == 0, 0, power)

  return(power)
}

# Linear-scale base implied by a log scale name.
log_base_value <- function(log_scale) {
  switch(log_scale,
    log = exp(1),
    log2 = 2,
    log10 = 10,
    stop("Unrecognized log scale: ", log_scale)
  )
}

# Minimum detectable fold change at a REQUESTED power level.
#
# Fixed in 0.0.3. Before that this function was called with the OBSERVED power
# vector rather than a requested power level, because a dplyr::mutate() bound
# `power = power` to the column created on the line above instead of to the
# formal argument. The stored value was therefore
# exp((z_alpha + z_observed) * SE): understated on most rows, smallest exactly
# where the data are weakest, and Inf on the most significant rows
# (qnorm(1) = Inf). It also hardcoded exp(), inflating log2 SEs, and its Inf
# guard tested `power == 0` alongside the SEs, so it only ever fired as a side
# effect of that same contamination.
#
# It now depends only on the standard errors, `alpha`, the residual degrees of
# freedom and `power` -- never on the observed effect. That effect-independence
# is what makes it eligible as a filter, or as the basis of a "we would have
# caught a change this large" claim.
#
# Uses t quantiles so that the detection limit and `delta_p_value` refer to the
# same distribution; with 6-13 samples per arm the normal approximation is
# optimistic. This is the normal-theory MDE formula with t quantiles
# substituted, not an exact noncentral-t solution -- close enough at these df,
# and slightly conservative.
#
# `base` is the linear-scale base matching the log scale of the SEs, so that a
# log2 contrast yields a log2-consistent fold change rather than an exp() one.
#
# Returns NA (never Inf, never 1) where no detection limit is defined: a zero
# or non-finite SE is a degenerate fit, not a 1-fold detection limit. Callers
# get the reason in `contrast_note`.
calculate_mdfc <- function(SE_x, SE_y, alpha = 0.05, power = 0.8,
                                    df = NULL, base = exp(1)) {
  se <- sqrt(SE_x^2 + SE_y^2)

  if (is.null(df) || length(df) != 1 || !is.finite(df) || df <= 0) {
    return(rep(NA_real_, length(se)))
  }

  MDFC <- base^((qt(1 - alpha / 2, df) + qt(power, df)) * se)
  MDFC[!is.finite(se) | se <= 0] <- NA_real_

  return(MDFC)
}

# Power to detect a change of a DECLARED size.
#
# The dual of calculate_mdfc(): that one fixes the power and reports
# the detectable effect, this one fixes the effect and reports the power. Both
# depend only on the standard errors, `alpha`, `df` and a declared quantity --
# never on the observed effect -- which is what makes either of them a fair
# answer to "were we powered to see a phenotype here?".
#
# `margin_log` is the effect size to be powered against, on the same log scale
# as the standard errors. Under the alternative that the true difference equals
# that margin, the t statistic follows a noncentral t with ncp = margin / SE,
# so two-sided power is the mass beyond +/- tcrit.
#
# Returns NA where SE or df make the question meaningless, matching
# calculate_mdfc().
calculate_power_at_margin <- function(SE_x, SE_y, alpha = 0.05, margin_log,
                                      df = NULL) {
  se <- sqrt(SE_x^2 + SE_y^2)

  if (is.null(df) || length(df) != 1 || !is.finite(df) || df <= 0) {
    return(rep(NA_real_, length(se)))
  }

  usable <- is.finite(se) & se > 0
  power <- rep(NA_real_, length(se))

  tcrit <- qt(1 - alpha / 2, df)
  ncp <- abs(margin_log) / se[usable]

  power[usable] <-
    (1 - pt(tcrit, df, ncp = ncp)) + pt(-tcrit, df, ncp = ncp)

  return(power)
}


get_current_base <- function(cond_a) {
  colnames <- colnames(cond_a)[grepl("log", colnames(cond_a))]
  base <- unlist(stringr::str_split(colnames[[1]], "_"))[1]
  return(base)
}

convert_base <- function(value, from_base, to_base) {
  if (from_base == "log2") {
    get(to_base)(2**value)
  } else if (from_base == "log10") {
    get(to_base)(10**value)
  } else {
    get(to_base)(exp(value))
  }
}

#' Compare two estimates of cell abundances from a Hooke model.
#'
#' @param ccm A cell_count_model.
#' @param cond_x tibble A cell type abundance estimate from estimate_abundances().
#' @param cond_y tibble A cell type abundance estimate from estimate_abundances().
#' @param by string The column name used to join the two estimates.
#' @param method string A method for correcting P-value multiple comparisons.
#'    This can be "BH" (Benjamini & Hochberg), "bonferroni" (Bonferroni),
#'    "hochberg" (Hochberg), "hommel", (Hommel), or "BYH" (Benjamini & Yekutieli).
#' @param alpha Desired significance level.
#' @param power Desired power level, used for `mdfc80`.
#' @param margin The effect size to be powered against, as a **fold change**
#'   (so `2` means a two-fold change), used for `power_at_margin`. Given as a
#'   fold change rather than on a log scale so that it means the same thing
#'   regardless of `log_scale`/`convert_scale`.
#' @param convert_scale Whether to convert to log2 scale.
#' @param adjust_q_values Whether to restrict multiple-testing correction to the
#'   rows with `power_at_margin >= power`. Changes what the FDR guarantee
#'   covers; see [adjust_q_values()]. Defaults to `FALSE`.
#' @param log_scale Log scale used for estimate_abundances.
#' @return tibble A table contrasting cond_x and cond_y (interpret as Y/X), with
#'   one row per `by` group. Columns of interest:
#'   \describe{
#'     \item{`delta_log_abund`, `delta_log_abund_se`}{Difference in abundance and
#'       its standard error, on the log scale given by `log_scale`/`convert_scale`.}
#'     \item{`delta_p_value`, `delta_q_value`}{Two-sided t test and its
#'       multiple-testing adjustment.}
#'     \item{`delta_log_abund_lo`, `delta_log_abund_hi`}{Wald confidence interval
#'       on `delta_log_abund` at level `1 - alpha`.}
#'     \item{`df_resid`}{Residual degrees of freedom of the fit, so a caller can
#'       rebuild the interval or an equivalence test from the table alone.}
#'     \item{`mdfc80`}{Minimum detectable fold change at the requested `power`:
#'       the smallest change the contrast could have caught. `NA` where no
#'       detection limit is defined.}
#'     \item{`power_at_margin`}{Power to detect a change of `margin`: the dual
#'       of `mdfc80`, fixing the effect and reporting the power rather than the
#'       reverse. `NA` on the same rows as `mdfc80`.}
#'     \item{`margin_fold_change`}{The `margin` used, so the table records what
#'       it was powered for.}
#'     \item{`contrast_note`}{Why `mdfc80` and `power_at_margin` are `NA`:
#'       `"degenerate_fit"` (zero or non-finite SE) or `"insufficient_df"`.
#'       `NA` when the row is fine.}
#'     \item{`observed_power`}{Post-hoc power, bit-identical to the column that
#'       was called `power` before 0.0.3. A deterministic restatement of
#'       `delta_p_value`, not a measure of precision. Do not filter on it.}
#'     \item{`power`}{**Deprecated** alias of `observed_power`, retained for one
#'       release and scheduled for removal.}
#'   }
#'
#'   `mdfc80` and `power_at_margin` are functions of `delta_log_abund_se`,
#'   `alpha`, `df_resid` and a declared constant only -- never of the observed
#'   effect. That effect-independence is what makes either of them a fair answer
#'   to "were we powered to see a phenotype here?", and what makes them safe to
#'   filter on.
#'
#' @details `power` and `mdfc` answered a different question than their names
#'   suggest, which is what 0.0.3 corrects. `power` substitutes
#'   the observed Wald statistic for the true effect, making it strictly
#'   increasing in `|Z|` and so carrying no information beyond `delta_p_value`
#'   -- `power >= 0.8` is exactly `p <= ~0.008`. It does not measure precision:
#'   a cell type with a huge standard error and a fluke estimate scores high,
#'   while a tightly measured genuine null scores the floor. `mdfc` inherits
#'   that contamination, because a `dplyr::mutate()` binds `power = power` to
#'   the column created on the line above rather than to the formal argument,
#'   making the stored value `exp((z_alpha + z_observed) * SE)` -- smallest
#'   exactly where the data are weakest, and `Inf` on the most significant rows.
#'
#'   As of 0.0.3, `power` survives only as a deprecated alias of
#'   `observed_power`, and the `mdfc` **column is gone**: `calculate_mdfc()` was
#'   fixed rather than frozen, and now emits `mdfc80`. `mdfc80` is therefore not
#'   a drop-in for `mdfc` -- it is a different quantity, it differs on
#'   essentially every row, and any threshold chosen against `mdfc` needs
#'   rechoosing. To state that a non-significant result excludes a meaningful
#'   change, compare `mdfc80` to your margin, read `power_at_margin`, or build
#'   an equivalence test from `delta_log_abund`, `delta_log_abund_se` and
#'   `df_resid`.
#' @importFrom dplyr full_join
#' @export
compare_abundances <- function(ccm,
                               cond_x,
                               cond_y,
                               by = "cell_group",
                               method = c("BH", "bonferroni", "hochberg", "hommel", "BY"),
                               alpha = 0.05,
                               power = 0.8,
                               margin = 2,
                               convert_scale = FALSE,
                               adjust_q_values = FALSE,
                               log_scale = c("log", "log10", "log2")) {
  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(tibble::is_tibble(cond_x))
  assertthat::assert_that(tibble::is_tibble(cond_y))
  assertthat::assert_that(is.character(by))

  assertthat::assert_that(
    tryCatch(
      expr = ifelse(match.arg(method) == "", TRUE, TRUE),
      error = function(e) FALSE
    ),
    msg = paste(
      'Argument method must be one of "size_factors",',
      '"BH", "bonferroni", "hochberg", "hommel", or "BY".'
    )
  )
  method <- match.arg(method)
  log_scale <- match.arg(log_scale)

  contrast_tbl <- dplyr::full_join(cond_x, cond_y, suffix = c("_x", "_y"), by = by)

  # num samples
  n <- nrow(model(ccm)$fitted)
  # num parameters
  k <- length(rownames(coef(ccm@best_full_model)))
  df.r <- n - k - 1

  # convert from log10 scale to log2
  if (convert_scale & log_scale == "log10") {
    contrast_tbl <- contrast_tbl %>%
      mutate(
        log_abund_x = log2(10^(log_abund_x)),
        log_abund_se_x = log2(10^(log_abund_se_x)),
        log_abund_sd_x = log2(10^(log_abund_sd_x))
      ) %>%
      mutate(
        log_abund_y = log2(10^(log_abund_y)),
        log_abund_se_y = log2(10^(log_abund_se_y)),
        log_abund_sd_y = log2(10^(log_abund_sd_y))
      )
  }
  # convert from natural log scale to log2
  else if (convert_scale & log_scale == "log") {
    contrast_tbl <- contrast_tbl %>%
      mutate(
        log_abund_x = log2(exp(log_abund_x)),
        log_abund_se_x = log2(exp(log_abund_se_x)),
        log_abund_sd_x = log2(exp(log_abund_sd_x))
      ) %>%
      mutate(
        log_abund_y = log2(exp(log_abund_y)),
        log_abund_se_y = log2(exp(log_abund_se_y)),
        log_abund_sd_y = log2(exp(log_abund_sd_y))
      )
  }

  # These scalars are bound OUTSIDE the mutate on purpose. dplyr::mutate()
  # evaluates sequentially against the data mask, so a line reading
  # `f(power = power)` binds to a `power` COLUMN created earlier in the same
  # mutate, not to the formal argument. That is precisely how `mdfc` came to be
  # a function of the observed effect instead of the requested power level.
  requested_power <- power
  requested_alpha <- alpha

  # Log scale the contrast columns are on once any conversion above has run,
  # and the linear base that matches it.
  result_log_scale <- if (convert_scale) "log2" else log_scale
  result_base <- log_base_value(result_log_scale)

  df_ok <- length(df.r) == 1 && is.finite(df.r) && df.r > 0
  tcrit <- if (df_ok) qt(1 - requested_alpha / 2, df.r) else NA_real_

  # `margin` is a FOLD CHANGE, so that it means the same thing whatever log
  # scale the estimates arrived on. Convert it into that scale here rather than
  # asking callers to remember whether a 2-fold change is 0.693 or 1.
  assertthat::assert_that(is.numeric(margin), length(margin) == 1, margin > 0)
  margin_log <- log(margin, base = result_base)

  contrast_tbl <- contrast_tbl %>%
    dplyr::mutate(
      delta_log_abund = log_abund_y - log_abund_x,
      delta_log_abund_se = sqrt(log_abund_se_y^2 + log_abund_se_x^2),
      tvalue = delta_log_abund / delta_log_abund_se,
      delta_p_value = 2 * pt(-abs(tvalue), df.r),
      delta_p_value = ifelse(is.nan(delta_p_value), 1, delta_p_value),
      # delta_p_value = pnorm(abs(delta_log_abund), sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2), lower.tail=FALSE),
      delta_q_value = p.adjust(delta_p_value, method = method),

      # Residual df, published so that a consumer can rebuild a confidence
      # interval or an equivalence test from the table alone. Previously
      # computed here and thrown away.
      df_resid = df.r,

      # Wald CI on the abundance difference, on the same log scale as
      # delta_log_abund.
      delta_log_abund_lo = delta_log_abund - tcrit * delta_log_abund_se,
      delta_log_abund_hi = delta_log_abund + tcrit * delta_log_abund_se,

      # Effect-INDEPENDENT minimum detectable fold change at the requested
      # `power`. This is the column to filter on, and the one that separates
      # "no phenotype" from "not powered".
      mdfc80 = calculate_mdfc(
        log_abund_se_x, log_abund_se_y,
        alpha = requested_alpha, power = requested_power,
        df = df.r, base = result_base
      ),

      # Why a row has no usable detection limit; NA means the row is fine.
      contrast_note = dplyr::case_when(
        rep(!df_ok, dplyr::n()) ~ "insufficient_df",
        !is.finite(delta_log_abund_se) | delta_log_abund_se <= 0 ~ "degenerate_fit",
        TRUE ~ NA_character_
      ),

      # Power to detect a change of `margin`, the dual of mdfc80. Also
      # effect-independent, so it is a fair answer to "were we powered here?"
      # and is safe to filter on.
      power_at_margin = calculate_power_at_margin(
        log_abund_se_x, log_abund_se_y,
        alpha = requested_alpha, margin_log = margin_log, df = df.r
      ),

      # The margin power_at_margin was computed against, so the table says what
      # it was powered for instead of leaving a reader to guess.
      margin_fold_change = margin,

      # Post-hoc power. A deterministic restatement of delta_p_value, not a
      # measure of precision; see the note above calculate_observed_power().
      # Kept because callers still want the number, under a name that says
      # what it is. Do not filter on it.
      observed_power = calculate_observed_power(
        log_abund_x, log_abund_se_x, log_abund_y, log_abund_se_y,
        alpha = requested_alpha
      ),

      # DEPRECATED alias for one release, bit-identical to hooke <= 0.0.2.
      # Scheduled for removal; use `observed_power`.
      power = observed_power
    ) %>%
    select(-tvalue)

  if (adjust_q_values) {
    warning(
      "adjust_q_values = TRUE restricts multiple-testing correction to rows ",
      "with power_at_margin >= ", requested_power, ". This changes what the ",
      "FDR guarantee covers: it applies to the adequately powered subset, not ",
      "to every row tested. Report it as such."
    )
    contrast_tbl <- adjust_q_values(contrast_tbl, power_threshold = requested_power)
  }

  return(contrast_tbl)
}

# Re-runs BH over only the rows whose `power_at_margin` clears a threshold,
# returning NA elsewhere.
#
# This is a defensible filter only because `power_at_margin` is
# effect-independent: it is a function of the standard errors, `alpha`, `df`
# and the declared margin, never of the observed effect. The pre-0.0.3 version
# of this function filtered on the old `power` column, which was a
# deterministic restatement of `delta_p_value` -- that selected rows by their
# p-value and then adjusted those same p-values, which is anti-conservative.
#
# It still changes the FDR guarantee: the correction now applies to the subset
# of rows that were adequately powered, not to all rows tested. That is a
# choice about what the guarantee covers, and callers should say so.
#
# Off by default and not enabled in production.
adjust_q_values <- function(contrast_tbl,
                            power_threshold = 0.8) {
  new_q_values <- contrast_tbl %>%
    mutate(rn = row_number()) %>%
    filter(power_at_margin >= power_threshold) %>%
    mutate(delta_q_value = p.adjust(delta_p_value, method = "BH")) %>%
    select(rn, delta_q_value)

  contrast_tbl <- contrast_tbl %>%
    mutate(rn = row_number()) %>%
    select(-delta_q_value) %>%
    left_join(new_q_values, by = "rn") %>%
    select(-rn)

  return(contrast_tbl)
}


#' @noRd
correlate_abundance_changes <- function(pln_model, cond_b_vs_a_tbl, edge_allowlist = NULL, edge_denylist = NULL) {
  cov_graph <- return_igraph(pln_model)
  cov_edges <- igraph::as_data_frame(cov_graph, what = "edges") %>% dplyr::filter(weight != 0.00)
  change_corr_tbl <- cov_edges %>%
    dplyr::select(from, to, weight) %>%
    dplyr::rename(pcor = weight)

  if (!is.null(edge_allowlist)) {
    edges_to_add <- edge_allowlist %>%
      dplyr::anti_join(change_corr_tbl, by = c("from", "to")) %>%
      dplyr::mutate(pcor = 0.001)
    
    change_corr_tbl <- dplyr::bind_rows(change_corr_tbl, edges_to_add)
  }
  if (!is.null(edge_denylist)) {
    change_corr_tbl <- change_corr_tbl %>%
      dplyr::anti_join(edge_denylist, by = c("from", "to"))
  }

  # corr_edge_coords_umap_delta_abund = corr_edge_coords_umap
  change_corr_tbl <- dplyr::left_join(change_corr_tbl, cond_b_vs_a_tbl %>% setNames(paste0("to_", names(.))), by = c("to" = "to_cell_group")) # %>%
  # dplyr::rename(log_abund_x,
  #              to_delta_log_abund = delta_log_abund)
  change_corr_tbl <- dplyr::left_join(change_corr_tbl, cond_b_vs_a_tbl %>% setNames(paste0("from_", names(.))), by = c("from" = "from_cell_group")) # %>%
  #  dplyr::rename(from_delta_log_abund = delta_log_abund)
  return(change_corr_tbl)
}

#' Helper function to plot kinetics
#' @param tp timepoint
#' @param perturbation_ccm a cell count model with a perturbation
#' @param interval_col column that matches the timepoint information
#' @param wt_pred_df control output from estimate_abundances_over_interval()
#' @param ko_pred_df perturbation output from estimate_abundances_over_interval()
#' @export
compare_ko_to_wt_at_timepoint <- function(tp, perturbation_ccm, wt_pred_df, ko_pred_df, interval_col) {
  cond_wt <- wt_pred_df %>% filter(!!sym(interval_col) == tp)
  cond_ko <- ko_pred_df %>% filter(!!sym(interval_col) == tp)
  return(compare_abundances(perturbation_ccm, cond_wt, cond_ko))
}
