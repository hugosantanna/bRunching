# ---------------------------------------------------------------------------- #
# Pairs bootstrap for the corrected CCN methods (uniform, tobit, het_tobit,
# symmetric).
#
# Each iteration resamples observations with replacement and re-runs the
# FULL pipeline:
#   (i)  re-bin Z on the resample (so quantile cuts adapt),
#   (ii) re-compute the censored expectation under the chosen method,
#   (iii) re-fit the corrected lm and harvest beta and delta.
#
# Iterations that fail (rank-deficient design, Tobit non-convergence, all
# observations dropped, etc.) are silently skipped and reported via B_ok.
# Mirrors the structure of cct_bootstrap() in R/cct.R.
# ---------------------------------------------------------------------------- #

#' Pairs bootstrap for the corrected CCN methods
#'
#' Re-runs `ccn_internal_fit()` on each resample. Returns SE and percentile
#' CIs for `beta` and `delta`.
#'
#' @param S,I,Z Outcome, treatment, instrument vectors (full sample).
#' @param method One of `"uniform"`, `"tobit"`, `"het_tobit"`, `"symmetric"`.
#' @param n_bins Number of instrument bins (already resolved upstream).
#' @param controls_mat,het_beta_mat,het_delta_mat Numeric matrices (or
#'   `NULL`) over the full sample.
#' @param cluster_delta,swap Pass-through arguments for `ccn_internal_fit`.
#' @param treatment Treatment variable name (for labelling beta).
#' @param beta_names,delta_names Coefficient names from the main fit; used
#'   to size and label the bootstrap result matrices so even iterations
#'   where the resample's design has different bin counts are aligned by
#'   common position (extras are dropped, gaps filled with NA).
#' @param B Number of bootstrap replications.
#' @param seed Optional integer seed.
#' @return A list with `beta_b`, `delta_b`, `se_beta`, `se_delta`,
#'   `ci_beta`, `ci_delta`, `B_ok`.
#' @keywords internal
ccn_bootstrap <- function(S, I, Z,
                          method,
                          n_bins,
                          controls_mat,
                          het_beta_mat,
                          het_delta_mat,
                          cluster_delta,
                          swap,
                          treatment,
                          beta_names,
                          delta_names,
                          B = 500L,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- length(S)
  p_beta  <- length(beta_names)
  p_delta <- length(delta_names)

  beta_b  <- matrix(NA_real_, nrow = B, ncol = p_beta,
                    dimnames = list(NULL, beta_names))
  delta_b <- matrix(NA_real_, nrow = B, ncol = p_delta,
                    dimnames = list(NULL, delta_names))

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)

    Sb <- S[idx]
    Ib <- I[idx]
    Zb <- Z[idx]
    Cb <- if (!is.null(controls_mat)) controls_mat[idx, , drop = FALSE] else NULL
    Hb <- if (!is.null(het_beta_mat))  het_beta_mat[idx, , drop = FALSE]  else NULL
    Db <- if (!is.null(het_delta_mat)) het_delta_mat[idx, , drop = FALSE] else NULL

    fit_b <- tryCatch(
      ccn_internal_fit(
        S             = Sb,
        I             = Ib,
        Z             = Zb,
        method        = method,
        n_bins        = n_bins,
        controls_mat  = Cb,
        het_beta_mat  = Hb,
        het_delta_mat = Db,
        cluster_delta = cluster_delta,
        swap          = swap,
        treatment     = treatment
      ),
      error   = function(e) NULL,
      warning = function(w) NULL
    )

    if (is.null(fit_b)) next
    if (is.null(fit_b$beta) || any(is.na(fit_b$beta))) next

    # Store beta by position. Names should match the main fit (treatment +
    # I_<het_beta_vars>) so direct positional assignment is safe.
    nb_use <- min(length(fit_b$beta), p_beta)
    beta_b[b, seq_len(nb_use)] <- fit_b$beta[seq_len(nb_use)]

    # Delta column count CAN vary across iterations (cluster_delta with
    # quantile cuts can yield different #bins on the resample). Match by
    # name when possible; otherwise fall back to positional alignment.
    if (!is.null(fit_b$delta) && p_delta > 0L) {
      d_names <- names(fit_b$delta)
      common <- intersect(delta_names, d_names)
      if (length(common) > 0L) {
        delta_b[b, common] <- fit_b$delta[common]
      } else {
        nd_use <- min(length(fit_b$delta), p_delta)
        delta_b[b, seq_len(nd_use)] <- fit_b$delta[seq_len(nd_use)]
      }
    }
  }

  ok <- !is.na(beta_b[, 1L])
  B_ok <- sum(ok)

  se_beta <- apply(beta_b, 2L, function(v) {
    vv <- v[!is.na(v)]
    if (length(vv) >= 2L) stats::sd(vv) else NA_real_
  })
  names(se_beta) <- beta_names

  if (p_delta > 0L) {
    se_delta <- apply(delta_b, 2L, function(v) {
      vv <- v[!is.na(v)]
      if (length(vv) >= 2L) stats::sd(vv) else NA_real_
    })
    names(se_delta) <- delta_names
  } else {
    se_delta <- NULL
  }

  ci_beta <- t(apply(beta_b, 2L, function(v) {
    vv <- v[!is.na(v)]
    if (length(vv) >= 10L) {
      stats::quantile(vv, probs = c(0.025, 0.975), names = FALSE, type = 7)
    } else c(NA_real_, NA_real_)
  }))
  rownames(ci_beta) <- beta_names
  colnames(ci_beta) <- c("2.5%", "97.5%")

  if (p_delta > 0L) {
    ci_delta <- t(apply(delta_b, 2L, function(v) {
      vv <- v[!is.na(v)]
      if (length(vv) >= 10L) {
        stats::quantile(vv, probs = c(0.025, 0.975), names = FALSE, type = 7)
      } else c(NA_real_, NA_real_)
    }))
    rownames(ci_delta) <- delta_names
    colnames(ci_delta) <- c("2.5%", "97.5%")
  } else {
    ci_delta <- NULL
  }

  list(
    beta_b   = beta_b,
    delta_b  = delta_b,
    se_beta  = se_beta,
    se_delta = se_delta,
    ci_beta  = ci_beta,
    ci_delta = ci_delta,
    B_ok     = B_ok
  )
}
