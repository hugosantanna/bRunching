# ---------------------------------------------------------------------------- #
# S3 methods for class "ccn"
# ---------------------------------------------------------------------------- #

#' @export
print.ccn <- function(x, ...) {
  cat("\nCCN Bunching Correction Estimator\n")
  cat("Method:", x$method, "\n")
  cat("Obs:", x$n_obs, " | Censored:", x$n_censored,
      sprintf(" (%.1f%%)\n", 100 * x$n_censored / x$n_obs))

  cat("\nTreatment effect (beta):\n")
  beta_se <- sqrt(diag(vcov(x$fit)))[names(x$beta)]

  # handle name mismatch: in the lm the var is "I", but we label with treatment name
  beta_names <- names(x$beta)
  fit_coefs <- coef(x$fit)
  fit_vcov  <- vcov(x$fit)

  for (i in seq_along(x$beta)) {
    nm <- beta_names[i]
    # find matching coef in the fit
    fit_nm <- if (nm == x$treatment) "I" else nm
    se_val <- if (fit_nm %in% names(fit_coefs)) sqrt(fit_vcov[fit_nm, fit_nm]) else NA
    cat(sprintf("  %-20s %10.5f  (SE: %.5f)\n", nm, x$beta[i], se_val))
  }

  if (!is.null(x$delta) && x$method != "naive") {
    cat("\nCorrection term (delta):\n")
    for (i in seq_along(x$delta)) {
      nm <- names(x$delta)[i]
      se_val <- if (nm %in% rownames(fit_vcov)) sqrt(fit_vcov[nm, nm]) else NA
      cat(sprintf("  %-20s %10.5f  (SE: %.5f)\n", nm, x$delta[i], se_val))
    }
  }

  if (x$method != "naive") {
    cat(sprintf("\nInstrument bins: %d", x$n_bins))
    if (x$n_switched > 0) {
      cat(sprintf(" | Fallback bins: %d", x$n_switched))
    }
    cat(sprintf("\nMean censored expectation: %.4f\n", x$cens_exp_mean))
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.ccn <- function(object, ...) {
  cat("\nCCN Bunching Correction Estimator\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("Call:\n")
  print(object$call)
  cat("\nMethod:", object$method, "\n")
  cat("Observations:", object$n_obs, "\n")
  cat("Censored:", object$n_censored,
      sprintf("(%.1f%%)\n", 100 * object$n_censored / object$n_obs))

  if (object$method != "naive") {
    cat("Instrument bins:", object$n_bins, "\n")
    if (object$n_switched > 0) {
      cat("Bins with fallback:", object$n_switched, "\n")
    }
    cat(sprintf("Mean censored expectation: %.4f\n", object$cens_exp_mean))
  }

  cat("\n")
  cat("--- Underlying regression ---\n")
  print(summary(object$fit))
  invisible(object)
}


#' @export
coef.ccn <- function(object, ...) {
  c(object$beta, object$delta)
}


#' @export
plot.ccn <- function(x, ...) {
  discontinuity_plot(
    outcome    = x$outcome,
    treatment  = x$treatment,
    data       = x$data,
    instrument = x$instrument,
    ...
  )
}
