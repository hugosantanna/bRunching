# ---------------------------------------------------------------------------- #
# Internal utilities
# ---------------------------------------------------------------------------- #

#' Kernel-weighted local linear regression (lpoly port)
#'
#' Evaluates a local linear (degree-1) kernel regression at specified points,
#' using an Epanechnikov kernel. This matches Stata's `lpoly` with `degree(1)`.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param eval_at Numeric vector of points at which to evaluate the fit.
#' @param h Bandwidth (scalar).
#' @return A list with `fit` (predicted values) and `se` (standard errors),
#'   each the same length as `eval_at`.
#' @keywords internal
lpoly_fit <- function(x, y, eval_at, h) {
  n_eval <- length(eval_at)
  fit <- rep(NA_real_, n_eval)
  se  <- rep(NA_real_, n_eval)

  for (j in seq_len(n_eval)) {
    x0 <- eval_at[j]
    u  <- (x - x0) / h

    # Epanechnikov kernel: K(u) = 0.75*(1 - u^2) for |u| <= 1
    w <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)

    idx <- which(w > 0)
    if (length(idx) < 3) next  # need at least 3 obs for lm + SE

    xw <- x[idx] - x0  # center at evaluation point
    yw <- y[idx]
    ww <- w[idx]

    # weighted local linear: y = a + b*(x - x0)
    wlm <- stats::lm.wfit(
      x = cbind(1, xw),
      y = yw,
      w = ww
    )

    fit[j] <- wlm$coefficients[1]  # intercept = fitted value at x0

    # SE: sqrt(diag(solve(X'WX) * sigma^2_w))
    # where sigma^2_w = sum(w * resid^2) / (sum(w) - 2)
    X <- cbind(1, xw)
    W <- diag(ww)
    XWX <- crossprod(X * ww, X)
    XWX_inv <- tryCatch(solve(XWX), error = function(e) NULL)
    if (is.null(XWX_inv)) next

    resid_w <- yw - X %*% wlm$coefficients
    sigma2 <- sum(ww * resid_w^2) / (sum(ww) - 2)
    var_beta <- XWX_inv * sigma2
    se[j] <- sqrt(max(var_beta[1, 1], 0))
  }

  list(fit = fit, se = se)
}
