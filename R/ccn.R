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
#'   `"het_tobit"` (per-bin Tobit), `"naive"`, or `"cct"` (Caetano-
#'   Caetano-Tecchio 2025; identifies the local boundary slope
#'   \eqn{\beta(0^+)} via a nonparametric location-scale censored
#'   regression of \eqn{X} on \eqn{Z}). When `method = "cct"`, the
#'   `instrument` argument is interpreted as the covariate \eqn{Z}
#'   indexing conditional quantiles (Assumption 5), not as an IV-style
#'   instrument.
#' @param n_bins Integer. Number of equal-frequency bins for discretizing the
#'   instrument. If `NULL`, uses the number of unique values (up to 50).
#'   Ignored when `method = "cct"` (Z is used as-is).
#' @param controls Character vector of control variable names, or a one-sided
#'   formula (e.g., `~ x1 + x2`).
#' @param het_beta Character vector of variable names to interact with the
#'   treatment for heterogeneous effects. If `NULL`, estimates a homogeneous
#'   effect. Not used by `method = "cct"`.
#' @param het_delta Character vector of variable names to interact with the
#'   correction term. Only used when `method != "naive"`. Not used by
#'   `method = "cct"`.
#' @param swap Fallback method when the symmetry estimator fails (>50\%
#'   censored in a bin). One of `"tobit"` (default) or `"het_tobit"`.
#' @param cluster_delta Integer or `NULL`. If provided, uses a separate
#'   discretization with this many bins for heterogeneous delta estimation
#'   (interacts correction term with instrument bin dummies).
#' @param alpha_star Known quantile level at which
#'   \eqn{H^{-1}(\alpha^\ast) = 0} (Assumption 5). Used only when
#'   `method = "cct"`. Default `0.5` (median). Identification (Assumption
#'   4) requires `alpha_star > max_z P(X = 0 | Z = z)`: per-cell bunching
#'   must not exceed `alpha_star` or the CDK location-scale regression
#'   has no signal in over-bunched cells. Pick `alpha_star` above the
#'   heaviest per-cell bunching rate; the paper uses 0.85 for its
#'   Section 6.2 application where overall bunching is around 50%.
#' @param alpha_grid Numeric grid of quantile levels strictly above
#'   `alpha_star` used in the Chen-Dahl-Khan location-scale regression.
#'   Defaults to 8 equally spaced points in `(alpha_star, 0.99]`, i.e.
#'   `seq(alpha_star, 0.99, length.out = 9)[-1]`. The length-based rule
#'   adapts to heavy bunching: when `alpha_star` is high (e.g. 0.85), the
#'   grid contracts toward `alpha_star` instead of collapsing to 1-2 very
#'   extreme quantiles where most within-z cells have `q_hat = 0`.
#'   `method = "cct"`.
#' @param alpha_s Scale-anchor quantile in `(alpha_star, 1)`. Defaults to
#'   `max(alpha_grid)`. `method = "cct"`.
#' @param locscale One of `"auto"` (default), `"on"`, `"off"`. Controls
#'   the Remark 4.2 short-circuit: when every within-z `alpha_star`-
#'   quantile is strictly positive, the location-scale regression can be
#'   skipped in favor of `m_hat(z) = q_hat(z; alpha_star)`. `"auto"`
#'   picks this automatically; `"off"` forces it; `"on"` always runs the
#'   full CDK regression. `method = "cct"`.
#' @param zeta_0 Trimming constant for the CDK third step (drops grid
#'   quantiles below this). Defaults to `1e-3 * sd(X)`. `method = "cct"`.
#' @param zeta_1 Kink location of the Buchinsky-Hahn smooth weight.
#'   Defaults to `zeta_0`. `method = "cct"`.
#' @param boot_B Integer. Number of pairs-bootstrap replications for the
#'   standard error. Set to `0` to skip. Method-dependent default: `500`
#'   for `method = "cct"` (no closed-form SE is available) and `0` for the
#'   corrected CCN methods (`uniform`, `tobit`, `het_tobit`, `symmetric`)
#'   so that the legacy OLS SE in `$fit` is the default; opt-in to the
#'   pairs bootstrap by passing e.g. `boot_B = 500`. `method = "naive"`
#'   ignores `boot_B` and reports OLS SE only.
#' @param seed Optional integer seed used inside the bootstrap. Applies to
#'   both `method = "cct"` and the corrected CCN methods.
#'
#' @return An object of class `"ccn"` containing:
#' \describe{
#'   \item{fit}{The fitted `lm` object from the corrected regression.}
#'   \item{beta}{Named vector of treatment effect estimates. For
#'     `method = "cct"` this is a scalar labelled `"beta(0+)"`.}
#'   \item{delta}{Named vector of correction term coefficients (or `NULL`
#'     when `method` is `"naive"` or `"cct"`).}
#'   \item{method}{The correction method used.}
#'   \item{estimand}{Either `"beta_global"` (CCN-style pooled slope) or
#'     `"beta(0+)"` (CCT local boundary slope).}
#'   \item{n_obs}{Number of observations.}
#'   \item{n_censored}{Number of censored (zero) observations.}
#'   \item{n_bins}{Number of instrument bins used.}
#'   \item{n_switched}{Number of bins where symmetry fallback was used.}
#'   \item{cens_exp_mean}{Mean censored expectation.}
#'   \item{call}{The matched call.}
#' }
#' For `method = "cct"` the object additionally stores `m_hat`
#' (per-observation estimate of \eqn{E[X^\ast \mid Z_i]}), `pi_hat`,
#' `Delta_over_pi`, `boot_se`, `boot_ci`, and `locscale_used`.
#'
#' @references
#' Caetano, C., Caetano, G., and Tecchio, O. (2025). "Correcting
#' Endogeneity of Treatments with Bunching Using Censoring Strategies."
#' Working paper.
#'
#' Chen, S., Dahl, G. B., and Khan, S. (2005). "Nonparametric
#' Identification and Estimation of a Censored Location-Scale Regression
#' Model." \emph{Journal of the American Statistical Association}
#' 100(469): 212-221.
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
#' \donttest{
#' # CCT 2025: local boundary slope beta(0+) via Chen-Dahl-Khan m(Z)
#' set.seed(1)
#' n <- 5000
#' K <- 20
#' Z <- sample.int(K, n, replace = TRUE)
#' mu_z <- -0.5 + 0.15 * Z
#' eta <- rnorm(n)
#' x_star <- mu_z + eta
#' x <- pmax(x_star, 0)
#' y <- 2 + 0.8 * x_star + 1.5 * x * (x > 0) + rnorm(n)
#' df <- data.frame(y = y, x = x, z = Z)
#'
#' fit <- ccn("y", "x", "z", df, method = "cct", boot_B = 100, seed = 42)
#' print(fit)
#' }
#'
#' @export
ccn <- function(outcome,
                treatment,
                instrument,
                data,
                method = c("symmetric", "uniform", "tobit", "het_tobit",
                           "naive", "cct"),
                n_bins = NULL,
                controls = NULL,
                het_beta = NULL,
                het_delta = NULL,
                swap = c("tobit", "het_tobit"),
                cluster_delta = NULL,
                alpha_star = 0.5,
                alpha_grid = NULL,
                alpha_s = NULL,
                locscale = c("auto", "on", "off"),
                zeta_0 = NULL,
                zeta_1 = NULL,
                boot_B = NULL,
                seed = NULL) {

  cl <- match.call()
  method   <- match.arg(method)
  swap     <- match.arg(swap)
  locscale <- match.arg(locscale)

  # Default boot_B is method-dependent so we preserve back-compat:
  #   * CCT: default 500 (the CCT estimator has no closed-form SE)
  #   * CCN methods: default 0 (the OLS SE in $fit is the legacy default;
  #     opt-in to bootstrap by passing boot_B = e.g. 500)
  #   * naive: ignored entirely.
  if (is.null(boot_B)) {
    boot_B <- if (method == "cct") 500L else 0L
  }

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

  # --- CCT branch: Caetano-Caetano-Tecchio 2025 ---
  if (method == "cct") {
    # Build controls matrix if requested.
    controls_mat <- NULL
    if (!is.null(controls)) {
      for (v in controls) {
        if (!v %in% names(data)) stop(sprintf("Control variable '%s' not found in data.", v))
      }
      controls_mat <- as.matrix(data[, controls, drop = FALSE])
    }

    # Defaults for CCT tuning parameters.
    if (is.null(alpha_grid)) {
      # 8 equally spaced points in (alpha_star, 0.99]. Adapts to heavy
      # bunching: for high alpha_star (e.g. 0.85), the grid contracts
      # toward alpha_star rather than collapsing to 1-2 extreme quantiles
      # where most within-z cells have q_hat = 0.
      alpha_grid <- seq(alpha_star, 0.99, length.out = 9L)[-1L]
    }
    if (is.null(alpha_s)) {
      alpha_s <- max(alpha_grid)
    }
    if (is.null(zeta_0)) {
      zeta_0 <- 1e-3 * stats::sd(I)
      if (!is.finite(zeta_0) || zeta_0 <= 0) zeta_0 <- 1e-6
    }
    if (is.null(zeta_1)) zeta_1 <- zeta_0

    return(ccn_cct(
      Y            = S,
      X            = I,
      Z            = Z,
      controls_mat = controls_mat,
      alpha_star   = alpha_star,
      alpha_grid   = alpha_grid,
      alpha_s      = alpha_s,
      locscale     = locscale,
      zeta_0       = zeta_0,
      zeta_1       = zeta_1,
      boot_B       = boot_B,
      seed         = seed,
      call         = cl,
      data         = data,
      outcome      = outcome,
      treatment    = treatment,
      instrument   = instrument
    ))
  }

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
      estimand    = "beta_global",
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

  # --- corrected non-naive branch: delegate to ccn_internal_fit ---
  het_beta_mat_full <- NULL
  if (!is.null(het_beta)) {
    for (v in het_beta) {
      if (!v %in% names(data)) stop(sprintf("het_beta variable '%s' not found in data.", v))
    }
    # Pass the raw het-beta covariate matrix (unmultiplied by I); the
    # internal fit recomputes the I_var interaction from it.
    het_beta_mat_full <- as.matrix(data[, het_beta, drop = FALSE])
    colnames(het_beta_mat_full) <- het_beta
  }

  het_delta_mat_full <- NULL
  if (!is.null(het_delta)) {
    for (v in het_delta) {
      if (!v %in% names(data)) stop(sprintf("het_delta variable '%s' not found in data.", v))
    }
    het_delta_mat_full <- as.matrix(data[, het_delta, drop = FALSE])
    colnames(het_delta_mat_full) <- het_delta
  }

  fit_obj <- ccn_internal_fit(
    S             = S,
    I             = I,
    Z             = Z,
    method        = method,
    n_bins        = n_bins,
    controls_mat  = ctrl_mat,
    het_beta_mat  = het_beta_mat_full,
    het_delta_mat = het_delta_mat_full,
    cluster_delta = cluster_delta,
    swap          = swap,
    treatment     = treatment
  )

  if (is.null(fit_obj)) {
    stop("ccn(): corrected fit failed (no usable observations after dropping ",
         "failed z-cells).")
  }

  # --- bootstrap inference (optional) ---
  boot <- NULL
  if (!is.null(boot_B) && boot_B > 0L) {
    boot <- ccn_bootstrap(
      S             = S,
      I             = I,
      Z             = Z,
      method        = method,
      n_bins        = n_bins,
      controls_mat  = ctrl_mat,
      het_beta_mat  = het_beta_mat_full,
      het_delta_mat = het_delta_mat_full,
      cluster_delta = cluster_delta,
      swap          = swap,
      treatment     = treatment,
      beta_names    = names(fit_obj$beta),
      delta_names   = names(fit_obj$delta),
      B             = as.integer(boot_B),
      seed          = seed
    )
  }

  result <- list(
    fit           = fit_obj$fit,
    beta          = fit_obj$beta,
    delta         = fit_obj$delta,
    method        = method,
    estimand      = "beta_global",
    n_obs         = fit_obj$n_obs,
    n_censored    = fit_obj$n_censored,
    n_bins        = n_bins,
    n_switched    = fit_obj$n_switched,
    cens_exp_mean = fit_obj$cens_exp_mean,
    call          = cl,
    data          = data,
    outcome       = outcome,
    treatment     = treatment,
    instrument    = instrument,
    swap          = swap,
    het_beta_vars = het_beta,
    het_delta_vars = het_delta,
    cluster_delta = cluster_delta,
    boot_B        = if (is.null(boot)) 0L else boot$B_ok,
    boot_se_beta  = if (is.null(boot)) NULL else boot$se_beta,
    boot_se_delta = if (is.null(boot)) NULL else boot$se_delta,
    boot_ci_beta  = if (is.null(boot)) NULL else boot$ci_beta,
    boot_ci_delta = if (is.null(boot)) NULL else boot$ci_delta,
    boot_beta     = if (is.null(boot)) NULL else boot$beta_b,
    boot_delta    = if (is.null(boot)) NULL else boot$delta_b
  )
  class(result) <- "ccn"
  result
}


# ---------------------------------------------------------------------------- #
# Internal fit for non-naive CCN methods.
#
# Encapsulates the full pipeline (censored expectation -> correction term ->
# bin dummies + interactions -> corrected lm). Used both by ccn() for the
# main fit and by ccn_bootstrap() for each resample.
#
# Returns NULL on degenerate input (e.g., no observations remain after
# dropping failed z-cells, or singular design); callers handle the NULL.
# ---------------------------------------------------------------------------- #
#' @keywords internal
ccn_internal_fit <- function(S, I, Z,
                             method,
                             n_bins,
                             controls_mat,
                             het_beta_mat,
                             het_delta_mat,
                             cluster_delta,
                             swap,
                             treatment) {

  cens <- as.integer(I == 0)

  # --- discretize instrument (re-binning per call so bootstrap resamples
  #     use a fresh quantile cut on the resampled Z) ---
  n_unique <- length(unique(Z))
  nb <- if (is.null(n_bins)) min(n_unique, 50L) else n_bins
  if (n_unique <= nb) {
    z_bins <- as.integer(as.factor(Z))
  } else {
    brks <- stats::quantile(Z,
                            probs = seq(0, 1, length.out = nb + 1L),
                            names = FALSE)
    # Guard against ties producing fewer than nb+1 unique breaks.
    brks <- unique(brks)
    if (length(brks) < 2L) {
      z_bins <- rep(1L, length(Z))
    } else {
      z_bins <- as.integer(cut(Z, breaks = brks,
                               include.lowest = TRUE, labels = FALSE))
    }
  }

  # --- compute censored expectations ---
  cens_exp <- switch(method,
    uniform   = cens_exp_uniform(I, z_bins, cens),
    tobit     = cens_exp_tobit_pooled(I, z_bins, cens),
    het_tobit = cens_exp_tobit(I, z_bins, cens),
    symmetric = cens_exp_symmetric(I, z_bins, cens, swap = swap),
    stop(sprintf("ccn_internal_fit: unsupported method '%s'", method))
  )

  dropped <- attr(cens_exp, "dropped")
  n_switched <- attr(cens_exp, "n_switched") %||% 0L

  keep <- !is.na(cens_exp)
  if (!is.null(dropped) && length(dropped) > 0L) {
    keep[dropped] <- FALSE
  }
  if (sum(keep) < 2L) return(NULL)

  S_k     <- S[keep]
  I_k     <- I[keep]
  cens_k  <- cens[keep]
  cens_exp_k <- cens_exp[keep]
  z_bins_k <- z_bins[keep]

  # correction term: E[I*|I=0]*1(I=0) + I
  reg_term <- cens_exp_k * cens_k + I_k
  cens_exp_mean <- mean(cens_exp_k, na.rm = TRUE)

  # bin dummies (recompute for kept obs; use ALL K levels for no-constant model)
  z_f_k <- factor(z_bins_k)
  if (nlevels(z_f_k) > 1L) {
    bin_dummies_k <- stats::model.matrix(~ 0 + z_f_k)
    colnames(bin_dummies_k) <- paste0("bin_", seq_len(ncol(bin_dummies_k)))
  } else {
    bin_dummies_k <- matrix(1, nrow = length(z_bins_k), ncol = 1L)
    colnames(bin_dummies_k) <- "bin_1"
  }

  het_beta_mat_k <- NULL
  if (!is.null(het_beta_mat)) {
    het_beta_mat_k <- het_beta_mat[keep, , drop = FALSE] * I_k
    colnames(het_beta_mat_k) <- paste0("I_", colnames(het_beta_mat))
  }

  ctrl_mat_k <- NULL
  if (!is.null(controls_mat)) {
    ctrl_mat_k <- controls_mat[keep, , drop = FALSE]
  }

  # --- build correction term interactions ---
  if (!is.null(cluster_delta)) {
    Z_k <- Z[keep]
    n_delta_bins <- cluster_delta
    if (length(unique(Z_k)) <= n_delta_bins) {
      z_delta <- as.integer(as.factor(Z_k))
    } else {
      brks_d <- stats::quantile(Z_k,
                                probs = seq(0, 1, length.out = n_delta_bins + 1L),
                                names = FALSE)
      brks_d <- unique(brks_d)
      if (length(brks_d) < 2L) {
        z_delta <- rep(1L, length(Z_k))
      } else {
        z_delta <- as.integer(cut(Z_k, breaks = brks_d,
                                  include.lowest = TRUE, labels = FALSE))
      }
    }
    delta_dummies <- stats::model.matrix(~ factor(z_delta))[, , drop = FALSE]
    delta_mat <- delta_dummies * reg_term
    colnames(delta_mat) <- paste0("delta_bin_", seq_len(ncol(delta_mat)))
  } else if (!is.null(het_delta_mat)) {
    delta_base <- cbind(reg_term = reg_term)
    delta_inter <- het_delta_mat[keep, , drop = FALSE] * reg_term
    colnames(delta_inter) <- paste0("reg_term_", colnames(het_delta_mat))
    delta_mat <- cbind(delta_base, delta_inter)
  } else {
    delta_mat <- cbind(reg_term = reg_term)
  }

  # --- run corrected regression (no constant) ---
  X <- cbind(I = I_k, het_beta_mat_k, bin_dummies_k, ctrl_mat_k, delta_mat)
  df_reg <- data.frame(S = S_k, X, check.names = FALSE)
  fit <- tryCatch(stats::lm(S ~ 0 + ., data = df_reg),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  co <- stats::coef(fit)
  beta <- co["I"]
  if (is.na(beta)) return(NULL)
  names(beta) <- treatment
  if (!is.null(het_beta_mat)) {
    het_coefs <- co[paste0("I_", colnames(het_beta_mat))]
    beta <- c(beta, het_coefs)
  }

  delta_names <- colnames(delta_mat)
  delta <- co[delta_names]

  list(
    fit           = fit,
    beta          = beta,
    delta         = delta,
    n_obs         = sum(keep),
    n_censored    = sum(cens_k),
    n_switched    = n_switched,
    cens_exp_mean = cens_exp_mean
  )
}

# base R does not have %||% before 4.4
`%||%` <- function(x, y) if (is.null(x)) y else x
