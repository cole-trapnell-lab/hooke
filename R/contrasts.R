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
    responses = Yc,
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
  assertthat::assert_that(is.numeric(min_log_abund))
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

    below_thresh <- log_abund < min_log_abund
    log_abund[below_thresh] <- min_log_abund
    # log_abund_se[below_thresh] = 0

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
      pln_model = pln_model
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
#' @return A tibble of cell abundance predictions.
#' @importFrom tibble tibble
#' @export
estimate_abundances_over_interval <- function(ccm, interval_start, interval_stop, interval_col = "timepoint", interval_step = 2, min_log_abund = -5, newdata = tibble()) {
  assertthat::assert_that(is(ccm, "cell_count_model"))
  assertthat::assert_that(is.numeric(interval_start))
  assertthat::assert_that(is.numeric(interval_stop))
  assertthat::assert_that(interval_stop >= interval_start)
  assertthat::assert_that(is.numeric(interval_step))

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
      min_log_abund = min_log_abund
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


# Calculates the power to detect a significant change in a given contrast
# Power is the probability of avoiding a Type II error
calculate_power <- function(beta_x, SE_x, beta_y, SE_y, alpha = 0.05) {
  # Wald test statistic for comparing two groups
  Z <- (beta_x - beta_y) / sqrt(SE_x^2 + SE_y^2)

  # Critical value for a two-tailed test at alpha
  Z_alpha <- qnorm(1 - alpha / 2)

  # Compute power: 1 - Type II error probability
  power <- 1 - pnorm(Z_alpha - Z) + pnorm(-Z_alpha - Z)
  
  # return 0 if divide by 0 
  power = ifelse(SE_x == 0 & SE_y == 0, 0, power)
  
  return(power)

}

# Estimates the smallest fold change you can detect at a given power level
calculate_mdfc <- function(SE_x, SE_y, alpha = 0.05, power = 0.8) {
  # Z-scores for significance and power
  Z_alpha <- qnorm(1 - alpha / 2) # Two-tailed
  Z_power <- qnorm(power)

  # Compute minimum detectable effect size (beta difference)
  delta_beta <- (Z_alpha + Z_power) * sqrt(SE_x^2 + SE_y^2)

  # Convert to fold change
  MDFC <- exp(delta_beta)
  
  # if power is 0, no fold change is detectable
  MDFC <- ifelse(power == 0 & SE_x == 0 & SE_y == 0, Inf, MDFC )

  return(MDFC)
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
#' @param alpha Desired significance level
#' @param power Desired power level for calculating minimum detectable fold change
#' @param convert_scale Whether to convert to log2 scale.
#' @param log_scale Log scale used for estimate_abundances.
#' @return tibble A table contrasting cond_x and cond_y (interpret as Y/X).
#' @importFrom dplyr full_join
#' @export
compare_abundances <- function(ccm,
                               cond_x,
                               cond_y,
                               by = "cell_group",
                               method = c("BH", "bonferroni", "hochberg", "hommel", "BY"),
                               alpha = 0.05,
                               power = 0.8,
                               convert_scale = FALSE,
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
  
  contrast_tbl <- contrast_tbl %>%
    dplyr::mutate(
      delta_log_abund = log_abund_y - log_abund_x,
      delta_log_abund_se = sqrt(log_abund_se_y^2 + log_abund_se_x^2),
      tvalue = delta_log_abund / delta_log_abund_se,
      delta_p_value = 2 * pt(-abs(tvalue), df.r),
      delta_p_value = ifelse(is.nan(delta_p_value), 1, delta_p_value),
      # delta_p_value = pnorm(abs(delta_log_abund), sd = sqrt(log_abund_se_y^2 + log_abund_se_x^2), lower.tail=FALSE),
      delta_q_value = p.adjust(delta_p_value, method = method),
      power = calculate_power(log_abund_x, log_abund_se_x, log_abund_y, log_abund_se_y, alpha = alpha),
      mdfc = calculate_mdfc(log_abund_se_x, log_abund_se_y, power = power)
    ) %>%
    select(-tvalue)

  return(contrast_tbl)
}

# redo multiple hypothesis correction 
# by only including cell coutns where we have power
# otherwise return NA 
adjust_q_values <- function(contrast_tbl, 
                            power_threshold = 0.8) {
  
  
  new_q_values = contrast_tbl %>% 
    mutate(rn = row_number()) %>% 
    filter(power >= power_threshold) %>% 
    mutate(delta_q_value = p.adjust(delta_p_value, method = "BH")) %>% 
    select(rn, delta_q_value)
  
  contrast_tbl <- contrast_tbl %>% 
    mutate(rn = row_number()) %>% 
    select(-delta_q_value) %>% 
    left_join(new_q_values, by = "rn")  %>% 
    select(-rn)
  
  return(contrast_tbl)
  
}


#' @noRd
correlate_abundance_changes <- function(pln_model, cond_b_vs_a_tbl) {
  cov_graph <- return_igraph(pln_model)
  cov_edges <- igraph::as_data_frame(cov_graph, what = "edges") %>% dplyr::filter(weight != 0.00)
  change_corr_tbl <- cov_edges %>%
    dplyr::select(from, to, weight) %>%
    dplyr::rename(pcor = weight)

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
