#' CCN Bunching Correction Estimator
#'
#' Estimates the causal effect of a treatment variable censored at zero on an
#' outcome, correcting for bunching bias using the Caetano, Caetano, and
#' Nielsen (CCN) methodology.
#'
#' @param outcome Character string. Name of the outcome variable in `data`.
#' @param treatment Character string. Name of the treatment variable (censored
#'   at zero) in `data`.
#' @param instrument Character string. Name of the instrument variable in
#'   `data`.
#' @param data A data.frame containing all variables.
#' @param method Correction method. One of `"symmetric"` (default),
#'   `"uniform"`, `"tobit"` (pooled Tobit with instrument-bin dummies),
#'   `"het_tobit"` (per-bin Tobit), or `"naive"`.
#' @param n_bins Integer. Number of equal-frequency bins for discretizing the
#'   instrument. If `NULL`, uses the number of unique values (up to 50).
#' @param controls Character vector of control variable names, or a one-sided
#'   formula (e.g., `~ x1 + x2`).
#' @param het_beta Character vector of variable names to interact with the
#'   treatment for heterogeneous effects. If `NULL`, estimates a homogeneous
#'   effect.
#' @param het_delta Character vector of variable names to interact with the
#'   correction term. Only used when `method != "naive"`.
#' @param swap Fallback method when the symmetry estimator fails (>50\%
#'   censored in a bin). One of `"tobit"` (default) or `"het_tobit"`.
#' @param cluster_delta Integer or `NULL`. If provided, uses a separate
#'   discretization with this many bins for heterogeneous delta estimation
#'   (interacts correction term with instrument bin dummies).
#'
#' @return An object of class `"ccn"` containing:
#' \describe{
#'   \item{fit}{The fitted `lm` object from the corrected regression.}
#'   \item{beta}{Named vector of treatment effect estimates.}
#'   \item{delta}{Named vector of correction term coefficients (if applicable).}
#'   \item{method}{The correction method used.}
#'   \item{n_obs}{Number of observations.}
#'   \item{n_censored}{Number of censored (zero) observations.}
#'   \item{n_bins}{Number of instrument bins used.}
#'   \item{n_switched}{Number of bins where symmetry fallback was used.}
#'   \item{cens_exp_mean}{Mean censored expectation.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(42)
#' n <- 2000
#' z <- runif(n, 0, 10)
#' x_star <- 2 + 0.5 * z + rnorm(n)
#' x <- pmax(x_star, 0)
#' y <- 1 + 0.8 * x_star + rnorm(n)
#'
#' df <- data.frame(y = y, x = x, z = z)
#'
#' # Naive (biased)
#' ccn("y", "x", "z", df, method = "naive")
#'
#' # Corrected (symmetric)
#' ccn("y", "x", "z", df, method = "symmetric")
#'
#' # With heterogeneous effects
#' df$w <- rbinom(n, 1, 0.5)
#' ccn("y", "x", "z", df, method = "symmetric", het_beta = "w")
#' }
#'
#' @export
ccn <- function(outcome,
                treatment,
                instrument,
                data,
                method = c("symmetric", "uniform", "tobit", "het_tobit", "naive"),
                n_bins = NULL,
                controls = NULL,
                het_beta = NULL,
                het_delta = NULL,
                swap = c("tobit", "het_tobit"),
                cluster_delta = NULL) {

  cl <- match.call()
  method <- match.arg(method)
  swap   <- match.arg(swap)

  # --- validate inputs ---
  stopifnot(is.data.frame(data))
  for (v in c(outcome, treatment, instrument)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found in data.", v))
  }

  # parse controls if formula
  if (inherits(controls, "formula")) {
    controls <- all.vars(controls)
  }

  # --- extract vectors ---
  S <- data[[outcome]]
  I <- data[[treatment]]
  Z <- data[[instrument]]

  # censoring indicator
  cens <- as.integer(I == 0)

  # --- discretize instrument ---
  n_unique <- length(unique(Z))
  if (is.null(n_bins)) {
    n_bins <- min(n_unique, 50L)
  }
  if (n_unique <= n_bins) {
    z_bins <- as.integer(as.factor(Z))
  } else {
    z_bins <- as.integer(cut(Z, breaks = stats::quantile(Z, probs = seq(0, 1, length.out = n_bins + 1L)),
                             include.lowest = TRUE, labels = FALSE))
  }

  # --- build design matrix ---
  # instrument bin dummies (drop first for identification; skip if only 1 bin)
  if (length(unique(z_bins)) > 1) {
    bin_dummies <- stats::model.matrix(~ factor(z_bins))[, -1, drop = FALSE]
    colnames(bin_dummies) <- paste0("bin_", seq_len(ncol(bin_dummies)))
  } else {
    bin_dummies <- NULL
  }

  # treatment interactions for heterogeneous beta
  het_beta_mat <- NULL
  if (!is.null(het_beta)) {
    for (v in het_beta) {
      if (!v %in% names(data)) stop(sprintf("het_beta variable '%s' not found in data.", v))
    }
    het_beta_mat <- as.matrix(data[, het_beta, drop = FALSE]) * I
    colnames(het_beta_mat) <- paste0("I_", het_beta)
  }

  # control variables
  ctrl_mat <- NULL
  if (!is.null(controls)) {
    for (v in controls) {
      if (!v %in% names(data)) stop(sprintf("Control variable '%s' not found in data.", v))
    }
    ctrl_mat <- as.matrix(data[, controls, drop = FALSE])
  }

  # --- naive model ---
  if (method == "naive") {
    X <- cbind(I = I, het_beta_mat, bin_dummies, ctrl_mat)
    df_reg <- data.frame(S = S, X, check.names = FALSE)
    fit <- stats::lm(S ~ ., data = df_reg)

    beta <- coef(fit)["I"]
    names(beta) <- treatment
    if (!is.null(het_beta)) {
      het_coefs <- coef(fit)[paste0("I_", het_beta)]
      beta <- c(beta, het_coefs)
    }

    result <- list(
      fit         = fit,
      beta        = beta,
      delta       = NULL,
      method      = method,
      n_obs       = nrow(data),
      n_censored  = sum(cens),
      n_bins      = n_bins,
      n_switched  = 0L,
      cens_exp_mean = NA_real_,
      call        = cl,
      data        = data,
      outcome     = outcome,
      treatment   = treatment,
      instrument  = instrument
    )
    class(result) <- "ccn"
    return(result)
  }

  # --- compute censored expectations ---
  cens_exp <- switch(method,
    uniform   = cens_exp_uniform(I, z_bins, cens),
    tobit     = cens_exp_tobit_pooled(I, z_bins, cens),
    het_tobit = cens_exp_tobit(I, z_bins, cens),
    symmetric = cens_exp_symmetric(I, z_bins, cens, swap = swap)
  )

  # handle dropped observations (from tobit failures)
  dropped <- attr(cens_exp, "dropped")
  n_switched <- attr(cens_exp, "n_switched") %||% 0L

  keep <- !is.na(cens_exp)
  if (!is.null(dropped) && length(dropped) > 0) {
    keep[dropped] <- FALSE
  }

  S_k     <- S[keep]
  I_k     <- I[keep]
  cens_k  <- cens[keep]
  cens_exp_k <- cens_exp[keep]
  z_bins_k <- z_bins[keep]

  # correction term: E[I*|I=0]*1(I=0) + I
  reg_term <- cens_exp_k * cens_k + I_k

  cens_exp_mean <- mean(cens_exp_k, na.rm = TRUE)

  # bin dummies (recompute for kept obs; handle single-level case)
  if (length(unique(z_bins_k)) > 1) {
    bin_dummies_k <- stats::model.matrix(~ factor(z_bins_k))[, -1, drop = FALSE]
    colnames(bin_dummies_k) <- paste0("bin_", seq_len(ncol(bin_dummies_k)))
  } else {
    bin_dummies_k <- NULL
  }

  het_beta_mat_k <- NULL
  if (!is.null(het_beta)) {
    het_beta_mat_k <- as.matrix(data[keep, het_beta, drop = FALSE]) * I_k
    colnames(het_beta_mat_k) <- paste0("I_", het_beta)
  }

  ctrl_mat_k <- NULL
  if (!is.null(controls)) {
    ctrl_mat_k <- as.matrix(data[keep, controls, drop = FALSE])
  }

  # --- build correction term interactions ---
  delta_mat <- NULL
  if (!is.null(cluster_delta)) {
    # interact reg_term with instrument bin dummies for separate delta per bin
    Z_k <- Z[keep]
    n_delta_bins <- cluster_delta
    if (length(unique(Z_k)) <= n_delta_bins) {
      z_delta <- as.integer(as.factor(Z_k))
    } else {
      z_delta <- as.integer(cut(Z_k, breaks = stats::quantile(Z_k, probs = seq(0, 1, length.out = n_delta_bins + 1L)),
                                include.lowest = TRUE, labels = FALSE))
    }
    delta_dummies <- stats::model.matrix(~ factor(z_delta))[, , drop = FALSE]  # keep all levels
    delta_mat <- delta_dummies * reg_term
    colnames(delta_mat) <- paste0("delta_bin_", seq_len(ncol(delta_mat)))
  } else if (!is.null(het_delta)) {
    for (v in het_delta) {
      if (!v %in% names(data)) stop(sprintf("het_delta variable '%s' not found in data.", v))
    }
    delta_base <- cbind(reg_term = reg_term)
    delta_inter <- as.matrix(data[keep, het_delta, drop = FALSE]) * reg_term
    colnames(delta_inter) <- paste0("reg_term_", het_delta)
    delta_mat <- cbind(delta_base, delta_inter)
  } else {
    delta_mat <- cbind(reg_term = reg_term)
  }

  # --- run corrected regression (no constant, as in Stata) ---
  X <- cbind(I = I_k, het_beta_mat_k, bin_dummies_k, ctrl_mat_k, delta_mat)
  df_reg <- data.frame(S = S_k, X, check.names = FALSE)
  fit <- stats::lm(S ~ 0 + ., data = df_reg)

  # --- extract coefficients ---
  beta <- coef(fit)["I"]
  names(beta) <- treatment
  if (!is.null(het_beta)) {
    het_coefs <- coef(fit)[paste0("I_", het_beta)]
    beta <- c(beta, het_coefs)
  }

  delta_names <- colnames(delta_mat)
  delta <- coef(fit)[delta_names]

  result <- list(
    fit           = fit,
    beta          = beta,
    delta         = delta,
    method        = method,
    n_obs         = sum(keep),
    n_censored    = sum(cens_k),
    n_bins        = n_bins,
    n_switched    = n_switched,
    cens_exp_mean = cens_exp_mean,
    call          = cl,
    data          = data,
    outcome       = outcome,
    treatment     = treatment,
    instrument    = instrument,
    swap          = swap,
    het_beta_vars = het_beta,
    het_delta_vars = het_delta,
    cluster_delta = cluster_delta
  )
  class(result) <- "ccn"
  result
}

# base R does not have %||% before 4.4
`%||%` <- function(x, y) if (is.null(x)) y else x
