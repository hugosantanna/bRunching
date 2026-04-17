# ---------------------------------------------------------------------------- #
# CCT-2025: Caetano, Caetano, Tecchio bunching-as-censoring estimator.
#
# Target parameter: beta(0+) := lim_{x|0} E[(Y(x) - Y(0)) / x | X = x].
#
# Algorithm A (Chen-Dahl-Khan): nonparametric location-scale estimation of
# m(Z) := E[X* | Z] from within-z sample quantiles of X | Z.
#
# Algorithm B (Proposition 5.1): under Assumption 7 (linear E[Y | X > 0]),
# beta(0+) is the coefficient on X in
#     Y ~ 1 + X + (X + 1{X = 0} * pi_hat)  (+ controls),
# where pi_hat = (E[m(Z)] - E[X]) / P(X = 0).
#
# This module is self-contained and does not reuse `R/censored_exp.R`.
# ---------------------------------------------------------------------------- #


#' Buchinsky-Hahn smooth approximation to the indicator of q > 0
#'
#' C^1 weight used in the CDK third-step weighted OLS. Vanishes below
#' `zeta_1`, equals 1 above `3 * zeta_1`, and is smooth in between.
#'
#' @param q Numeric vector of points at which to evaluate the weight.
#' @param zeta_1 Positive scalar kink location.
#' @return Numeric vector of weights the same length as `q`.
#' @keywords internal
cct_weight <- function(q, zeta_1) {
  # zeta_1 must be positive
  z1 <- zeta_1
  out <- numeric(length(q))

  # region 1: q <= zeta_1 -> weight is 0 (leave as initialized)
  # region 3: q >= 3 * zeta_1 -> weight is 1
  upper <- q >= 3 * z1
  out[upper] <- 1

  # region 2: zeta_1 < q < 3 * zeta_1 -> smooth transition
  mid <- q > z1 & q < 3 * z1
  if (any(mid)) {
    qm <- q[mid]
    num1 <- exp(qm - 2 * z1) / (1 + exp(qm - 2 * z1)) -
      exp(-z1) / (1 + exp(-z1))
    scale <- (2 + exp(z1) + exp(-z1)) / (exp(z1) - exp(-z1))
    out[mid] <- num1 * scale
  }

  out
}


#' Within-Z sample quantiles of X | Z = z for a grid of alpha levels
#'
#' Returns a matrix of within-cell sample quantiles. Because Z is assumed
#' to have finite support (Assumption 9), the conditional quantile
#' minimizer of the check loss reduces to the unconditional sample
#' quantile within each z-cell.
#'
#' @param X Numeric vector (treatment, censored at zero).
#' @param Z Vector (categorical / integer) of covariate values.
#' @param alphas Numeric vector of quantile levels in (0, 1).
#' @return A list with `z_levels` (unique Z values, in order) and
#'   `q` (matrix with one row per z-level and one column per alpha).
#' @keywords internal
cct_quantiles <- function(X, Z, alphas) {
  z_levels <- sort(unique(Z))
  K <- length(z_levels)
  L <- length(alphas)

  q <- matrix(NA_real_, nrow = K, ncol = L,
              dimnames = list(as.character(z_levels), as.character(alphas)))

  for (k in seq_len(K)) {
    zk <- z_levels[k]
    xk <- X[Z == zk]
    if (length(xk) == 0L) next
    # type = 1 is the check-loss minimizer (inverse-ecdf / empirical quantile)
    q[k, ] <- stats::quantile(xk, probs = alphas, type = 1, names = FALSE)
  }

  list(z_levels = z_levels, q = q)
}


#' Estimate m(Z) via the Chen-Dahl-Khan (2005) three-step procedure
#'
#' Implements Algorithm A from the CCT-2025 strategy memo, including the
#' Remark 4.2 short-circuit: when every alpha_star-quantile is strictly
#' positive, the location-scale regression is skipped and
#' `m_hat(z) = q_hat(z; alpha_star)` is returned.
#'
#' @param X Numeric vector (treatment).
#' @param Z Categorical / integer vector indexing the conditioning cells.
#' @param alpha_star Known quantile with `H^{-1}(alpha_star) = 0`.
#' @param alpha_grid Numeric grid of quantile levels strictly above
#'   `alpha_star` used in the location-scale step.
#' @param alpha_s Scale-anchor quantile in `(alpha_star, 1)`.
#' @param zeta_0 Trimming constant (drop alpha_ell quantiles below this).
#' @param zeta_1 Kink location of the smooth weight.
#' @param locscale One of `"auto"`, `"on"`, `"off"`.
#' @return A list with `m_hat` (per-obs vector of length `length(X)`),
#'   `m_hat_by_z` (named numeric), `z_levels`, and `locscale_used`
#'   (`"on"` or `"off"`).
#' @keywords internal
cct_m_hat <- function(X, Z,
                      alpha_star,
                      alpha_grid,
                      alpha_s,
                      zeta_0,
                      zeta_1,
                      locscale = c("auto", "on", "off")) {
  locscale <- match.arg(locscale)

  # --- Step 1: conditional quantile estimation ---
  all_alphas <- sort(unique(c(alpha_star, alpha_grid, alpha_s)))
  qres <- cct_quantiles(X, Z, all_alphas)
  z_levels <- qres$z_levels
  Q <- qres$q

  col_star <- match(alpha_star, all_alphas)
  col_s    <- match(alpha_s,    all_alphas)
  cols_ell <- match(alpha_grid, all_alphas)

  q_star <- Q[, col_star]
  q_s    <- Q[, col_s]
  Q_ell  <- Q[, cols_ell, drop = FALSE]

  # Identification check (Assumption 4). CDK's location-scale estimator
  # needs q_hat(z; alpha_star) > 0 in the z-cells that contribute — else
  # m_hat(z) collapses to 0 in those cells and Eq. (9) is biased.
  n_zero_star <- sum(q_star <= zeta_0)
  if (n_zero_star > 0L) {
    frac <- n_zero_star / length(q_star)
    if (frac > 0.2) {
      warning(sprintf(
        paste0("%d of %d z-cells (%.0f%%) have q_hat(z; alpha_star=%.2f) <= zeta_0. ",
               "Per-cell bunching exceeds alpha_star in these cells, violating ",
               "Assumption 4 of Caetano-Caetano-Tecchio (2025). Consider raising ",
               "alpha_star (typical fix: alpha_star = 1 - max per-cell bunching)."),
        n_zero_star, length(q_star), 100 * frac, alpha_star))
    }
  }

  # --- Remark 4.2 short-circuit: m_hat(z) = q_hat(z; alpha_star) ---
  use_locscale <- switch(
    locscale,
    "off"  = FALSE,
    "on"   = TRUE,
    "auto" = !all(q_star > 0)
  )

  if (!use_locscale) {
    m_by_z <- q_star
    names(m_by_z) <- as.character(z_levels)
    m_hat <- m_by_z[match(Z, z_levels)]
    names(m_hat) <- NULL
    return(list(
      m_hat         = m_hat,
      m_hat_by_z    = m_by_z,
      z_levels      = z_levels,
      locscale_used = "off"
    ))
  }

  # --- Step 2: compute C_ell for each grid quantile ---
  # Map each observation to its cell's q_star and q_s via z index.
  z_idx <- match(Z, z_levels)
  q_star_i <- q_star[z_idx]
  q_s_i    <- q_s[z_idx]
  w_i      <- cct_weight(q_star_i, zeta_1)

  denom <- q_s_i - q_star_i
  # Guard against zero scale at obs-level (when q_s == q_star); these obs
  # get zero weight anyway if q_star is near zero, but be defensive.
  safe <- denom > 0 & w_i > 0
  w_sum <- sum(w_i[safe])

  L <- length(alpha_grid)
  C_ell <- numeric(L)
  if (w_sum > 0) {
    for (ell in seq_len(L)) {
      q_ell_i <- Q_ell[z_idx, ell]
      ratio <- (q_ell_i - q_star_i) / denom
      C_ell[ell] <- sum(w_i[safe] * ratio[safe]) / w_sum
    }
  } else {
    # Degenerate: no informative cells. Fall back to short-circuit.
    m_by_z <- q_star
    names(m_by_z) <- as.character(z_levels)
    m_hat <- m_by_z[match(Z, z_levels)]
    names(m_hat) <- NULL
    return(list(
      m_hat         = m_hat,
      m_hat_by_z    = m_by_z,
      z_levels      = z_levels,
      locscale_used = "off"
    ))
  }

  # --- Step 3: weighted OLS per z-cell ---
  # For each z, stack rows ell with q_hat(z; alpha_ell) >= zeta_0.
  # Design row: (1, C_ell); response: q_hat(z; alpha_ell).
  # m_hat(z) = first coefficient.
  K <- length(z_levels)
  m_by_z <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    q_row <- Q_ell[k, ]
    keep <- q_row >= zeta_0
    if (sum(keep) < 2L) {
      # Not enough valid grid points to identify (1, C_ell); fall back
      # to q_star if it is informative, otherwise NA.
      m_by_z[k] <- if (q_star[k] > 0) q_star[k] else NA_real_
      next
    }
    D <- cbind(1, C_ell[keep])
    y <- q_row[keep]
    # Solve normal equations; fall back to q_star on singularity.
    sol <- tryCatch(
      qr.solve(crossprod(D), crossprod(D, y)),
      error = function(e) NULL
    )
    if (is.null(sol)) {
      m_by_z[k] <- if (q_star[k] > 0) q_star[k] else NA_real_
    } else {
      m_by_z[k] <- sol[1L]
    }
  }

  names(m_by_z) <- as.character(z_levels)
  m_hat <- m_by_z[match(Z, z_levels)]
  names(m_hat) <- NULL

  list(
    m_hat         = m_hat,
    m_hat_by_z    = m_by_z,
    z_levels      = z_levels,
    locscale_used = "on"
  )
}


#' Algorithm B (linear form): compute beta(0+) via OLS on a generated regressor
#'
#' Under Assumption 7 (Proposition 5.1), beta(0+) is the coefficient on X
#' in `Y ~ 1 + X + Xtilde (+ controls)` where
#' `Xtilde = X + 1{X = 0} * pi_hat` and
#' `pi_hat = sum(m_hat - X) / sum(X == 0)`.
#'
#' @param Y Outcome vector.
#' @param X Treatment vector.
#' @param m_hat Per-observation estimate of \eqn{E[X^\ast \mid Z_i]}.
#' @param controls Optional matrix of control variables.
#' @return A list with `beta` (scalar), `pi_hat`, `Delta_over_pi`,
#'   `gamma0`, and `fit` (the underlying `lm`).
#' @keywords internal
cct_beta <- function(Y, X, m_hat, controls = NULL) {
  stopifnot(length(Y) == length(X), length(X) == length(m_hat))

  cens <- as.integer(X == 0)
  n_cens <- sum(cens)
  if (n_cens == 0L) {
    stop("cct_beta: no censored observations (sum(X == 0) == 0); ",
         "the CCT estimand is not identified without bunching.")
  }

  # drop observations where m_hat is NA (z-cell where Algorithm A failed)
  keep <- !is.na(m_hat) & !is.na(Y) & !is.na(X)
  if (sum(keep) < 3L) {
    stop("cct_beta: fewer than 3 usable observations after dropping NAs.")
  }

  Yk <- Y[keep]; Xk <- X[keep]; mk <- m_hat[keep]
  cens_k <- as.integer(Xk == 0)
  n_cens_k <- sum(cens_k)
  if (n_cens_k == 0L) {
    stop("cct_beta: after dropping NAs, no censored observations remain.")
  }

  pi_hat <- sum(mk - Xk) / n_cens_k
  Xtilde <- Xk + cens_k * pi_hat

  dat <- data.frame(Y = Yk, X = Xk, Xtilde = Xtilde, check.names = FALSE)
  if (!is.null(controls)) {
    if (nrow(controls) != length(Y)) {
      stop("cct_beta: controls must have as many rows as Y.")
    }
    Ck <- controls[keep, , drop = FALSE]
    # Guard against non-matrix input
    if (!is.matrix(Ck)) Ck <- as.matrix(Ck)
    # Preserve column names; rename blanks to avoid lm collisions.
    cn <- colnames(Ck)
    if (is.null(cn)) cn <- paste0("W", seq_len(ncol(Ck)))
    cn[cn == ""] <- paste0("W", which(cn == ""))
    colnames(Ck) <- make.names(cn, unique = TRUE)
    dat <- cbind(dat, as.data.frame(Ck, check.names = FALSE))
  }

  fit <- stats::lm(Y ~ ., data = dat)
  co <- stats::coef(fit)

  beta_hat <- unname(co["X"])
  delta_over_pi <- unname(co["Xtilde"])
  gamma0 <- unname(co["(Intercept)"])

  list(
    beta           = beta_hat,
    pi_hat         = pi_hat,
    Delta_over_pi  = delta_over_pi,
    gamma0         = gamma0,
    fit            = fit,
    n_used         = sum(keep)
  )
}


#' Pairs bootstrap for beta(0+)
#'
#' Re-runs Algorithm A and Algorithm B on each resample. Iterations where
#' CDK / OLS degenerate (rank-deficient design, all-zero cells, etc.) are
#' silently skipped and reported via `B_ok`.
#'
#' @param Y,X,Z As in `ccn_cct`.
#' @param controls Matrix of controls (or NULL).
#' @param alpha_star,alpha_grid,alpha_s,zeta_0,zeta_1,locscale Pass-through
#'   parameters for `cct_m_hat`.
#' @param B Number of bootstrap replications.
#' @param seed Optional integer seed.
#' @return A list with `beta_b`, `se`, `ci_95`, `B_ok`.
#' @keywords internal
cct_bootstrap <- function(Y, X, Z,
                          controls,
                          alpha_star, alpha_grid, alpha_s,
                          zeta_0, zeta_1, locscale,
                          B = 500L,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(Y)
  beta_b <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Yb <- Y[idx]; Xb <- X[idx]; Zb <- Z[idx]
    Cb <- if (!is.null(controls)) controls[idx, , drop = FALSE] else NULL

    # Algorithm A must have at least 2 unique z-cells to identify m(Z).
    if (length(unique(Zb)) < 2L) next

    m_res <- tryCatch(
      cct_m_hat(Xb, Zb,
                alpha_star = alpha_star,
                alpha_grid = alpha_grid,
                alpha_s    = alpha_s,
                zeta_0     = zeta_0,
                zeta_1     = zeta_1,
                locscale   = locscale),
      error = function(e) NULL
    )
    if (is.null(m_res)) next

    fit_b <- tryCatch(
      cct_beta(Yb, Xb, m_res$m_hat, controls = Cb),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next
    if (is.na(fit_b$beta)) next

    beta_b[b] <- fit_b$beta
  }

  ok <- !is.na(beta_b)
  B_ok <- sum(ok)
  se <- if (B_ok >= 2L) stats::sd(beta_b[ok]) else NA_real_
  ci <- if (B_ok >= 10L) {
    stats::quantile(beta_b[ok], probs = c(0.025, 0.975),
                    names = FALSE, type = 7)
  } else c(NA_real_, NA_real_)

  list(beta_b = beta_b, se = se, ci_95 = ci, B_ok = B_ok)
}


#' CCT-2025 top-level driver (called from `ccn(method = "cct")`)
#'
#' Not exported. See `?ccn` for the user-facing API.
#'
#' @keywords internal
ccn_cct <- function(Y, X, Z, controls_mat,
                    alpha_star, alpha_grid, alpha_s,
                    locscale,
                    zeta_0, zeta_1,
                    boot_B, seed,
                    call,
                    data, outcome, treatment, instrument) {

  # --- input validation ---
  if (is.null(Z)) {
    stop("method = 'cct' requires a non-NULL instrument (interpreted as Z ",
         "in Assumption 5, not as an IV-style instrument).")
  }
  n_unique_z <- length(unique(Z))
  if (n_unique_z < 3L) {
    stop(sprintf(
      "method = 'cct' requires at least 3 unique values of Z (got %d); ",
      n_unique_z
    ), "Assumption 4 needs cross-cell variation.")
  }
  if (!(alpha_star > 0 && alpha_star < 1)) {
    stop("alpha_star must be strictly between 0 and 1.")
  }
  if (any(alpha_grid <= alpha_star) || any(alpha_grid >= 1)) {
    stop("alpha_grid must lie strictly in (alpha_star, 1).")
  }
  if (!(alpha_s > alpha_star && alpha_s < 1)) {
    stop("alpha_s must lie strictly in (alpha_star, 1).")
  }
  if (zeta_0 <= 0 || zeta_1 <= 0) {
    stop("zeta_0 and zeta_1 must be strictly positive.")
  }

  n <- length(X)
  n_cens <- sum(X == 0)
  bunch_frac <- n_cens / n
  if (bunch_frac < 0.01) {
    warning(sprintf(
      "Bunching fraction is %.2f%% (< 1%%); the CCT estimator may be unstable.",
      100 * bunch_frac
    ))
  }
  if (bunch_frac > 0.80) {
    warning(sprintf(
      "Bunching fraction is %.2f%% (> 80%%); the CCT estimator may be unstable.",
      100 * bunch_frac
    ))
  }

  # --- Algorithm A: m_hat(Z) ---
  m_res <- cct_m_hat(X, Z,
                     alpha_star = alpha_star,
                     alpha_grid = alpha_grid,
                     alpha_s    = alpha_s,
                     zeta_0     = zeta_0,
                     zeta_1     = zeta_1,
                     locscale   = locscale)

  # --- Algorithm B: beta(0+) ---
  fit_b <- cct_beta(Y, X, m_res$m_hat, controls = controls_mat)

  # --- bootstrap inference (optional) ---
  boot <- NULL
  if (!is.null(boot_B) && boot_B > 0L) {
    boot <- cct_bootstrap(Y, X, Z,
                          controls    = controls_mat,
                          alpha_star  = alpha_star,
                          alpha_grid  = alpha_grid,
                          alpha_s     = alpha_s,
                          zeta_0      = zeta_0,
                          zeta_1      = zeta_1,
                          locscale    = locscale,
                          B           = boot_B,
                          seed        = seed)
  }

  beta_vec <- fit_b$beta
  names(beta_vec) <- "beta(0+)"

  result <- list(
    fit           = fit_b$fit,
    beta          = beta_vec,
    delta         = NULL,
    method        = "cct",
    estimand      = "beta(0+)",
    n_obs         = n,
    n_censored    = n_cens,
    n_bins        = n_unique_z,
    n_switched    = 0L,
    cens_exp_mean = NA_real_,
    m_hat         = m_res$m_hat,
    m_hat_by_z    = m_res$m_hat_by_z,
    pi_hat        = fit_b$pi_hat,
    Delta_over_pi = fit_b$Delta_over_pi,
    gamma0        = fit_b$gamma0,
    locscale_used = m_res$locscale_used,
    alpha_star    = alpha_star,
    alpha_grid    = alpha_grid,
    alpha_s       = alpha_s,
    zeta_0        = zeta_0,
    zeta_1        = zeta_1,
    boot_B        = if (is.null(boot)) 0L else boot$B_ok,
    boot_se       = if (is.null(boot)) NA_real_ else boot$se,
    boot_ci       = if (is.null(boot)) c(NA_real_, NA_real_) else boot$ci_95,
    boot_beta     = if (is.null(boot)) NULL else boot$beta_b,
    call          = call,
    data          = data,
    outcome       = outcome,
    treatment     = treatment,
    instrument    = instrument
  )
  class(result) <- "ccn"
  result
}
